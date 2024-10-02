import gc

gc.collect()

import patsy
import multiprocessing
import pandas as pd
import numpy as np
import doubleml as dml
from sklearn.linear_model import LinearRegression, LogisticRegression, ElasticNet
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from scipy.stats import norm


def set_model(model, param, discrete_outcome):
    model_lower = model.lower()
    if model_lower in {
        "default",
        "linear",
        "logistic",
        "l",
        "d",
    }:
        model_set = LogisticRegression() if discrete_outcome else LinearRegression()

    elif model_lower in {"regularization", "r"}:
        model_set = ElasticNet(**param)

    elif model_lower in {
        "randomforest",
        "random forest",
        "random_forest",
        "rf",
        "forest",
    }:
        model_set = (
            RandomForestClassifier(**param)
            if discrete_outcome
            else RandomForestRegressor(**param)
        )

    elif model_lower in {
        "boosting",
        "gradient_boosting",
        "gradient boosting",
        "hist_gradient_boosting",
        "hist gradient boosting",
        "boost",
        "gradient_boost",
        "gradient boost",
        "hist_gradient_boost",
        "hist gradient boost",
        "b",
        "gb",
        "hgb",
    }:
        model_set = (
            HistGradientBoostingClassifier(**param)
            if discrete_outcome
            else HistGradientBoostingRegressor(**param)
        )

    elif model_lower in {"network", "neural_network", "neural network", "nn"}:
        if "max_iter" not in param:
            param["max_iter"] = 1000

        model_set = (
            MLPClassifier(**param) if discrete_outcome else MLPRegressor(**param)
        )

    else:
        raise Exception(
            "'model_y' and 'model_t' should be one of 'linear', 'regularization', 'rf', 'hgb', and 'nn'."
        )

    return model_set


def marginal_effect_for_treatment(
    df,
    Y,
    D,
    X,
    Z,
    model_y="rf",
    param_y={},
    param_grid_y={},
    scoring_y=None,
    model_t="rf",
    param_t={},
    param_grid_t={},
    scoring_t=None,
    CV=False,
    n_folds=10,
    n_jobs=-1,
    cf_n_folds=5,
    cf_n_rep=1,
    gate=False,
):
    ### Initialization
    n_folds = int(n_folds)
    n_jobs = int(n_jobs)
    cf_n_folds = int(cf_n_folds)
    cf_n_rep = int(cf_n_rep)
    if n_jobs == -1:
        n_jobs = multiprocessing.cpu_count()

    if sorted(list(df[Y].unique())) == [0, 1]:
        discrete_outcome = True
    else:
        discrete_outcome = False

    if len(df[D].unique()) > 5:
        discrete_treatment = False
    else:
        discrete_treatment = True

    if type(Z) is str:
        Z = [Z]
    elif Z is None:
        Z = []
    covariates = Z.copy()
    covariates.append(X)

    df[D] = df[D].astype("float")

    model_y_set = set_model(
        model_y,
        param_y,
        discrete_outcome,
    )

    model_t_set = set_model(
        model_t,
        param_t,
        discrete_treatment,
    )

    data_dml_base = dml.DoubleMLData(
        df[[X, Y, D, *Z]],
        y_col=Y,
        d_cols=D,
        x_cols=covariates,
    )

    if discrete_treatment:
        dml_model = dml.DoubleMLIRM(
            data_dml_base,
            ml_g=model_y_set,
            ml_m=model_t_set,
            n_folds=cf_n_folds,
            n_rep=cf_n_rep,
        )
        model_y_key = "ml_g"

    else:
        dml_model = dml.DoubleMLPLR(
            data_dml_base,
            ml_l=model_y_set,
            ml_m=model_t_set,
            n_folds=cf_n_folds,
            n_rep=cf_n_rep,
        )
        model_y_key = "ml_l"

    params = {"model.y": None, "model.t": None}
    if CV:
        dml_model.tune(
            param_grids={model_y_key: param_grid_y, "ml_m": param_grid_t},
            n_folds_tune=n_folds,
            n_jobs_cv=n_jobs,
            scoring_methods={model_y_key: scoring_y, "ml_m": scoring_t},
            search_mode="grid_search",
            return_tune_res=True,
        )
        if "ml_g" in dml_model.params:
            params["model.y"] = dml_model.params["ml_g"][D]
        elif "ml_g1" in dml_model.params:
            params["model.y"] = dml_model.params["ml_g1"][D]
        params["model.t"] = dml_model.params["ml_m"][D]

    dml_model.fit()

    df_gate = pd.DataFrame()
    if gate:
        groups = pd.get_dummies(df[X])
        gate = dml_model.gate(groups=groups)
        df_gate = gate.confint(level=0.95, joint=True, n_rep_boot=2000)

        np_basis = pd.DataFrame(
            np.diag(v=np.full((gate._basis.shape[1]), True))
        ).to_numpy()
        df_gate["X"] = df_gate.index
        df_gate["ME"] = df_gate["effect"]
        df_gate["sd"] = np.sqrt(
            (np.dot(np_basis, gate._blp_omega) * np_basis).sum(axis=1)
        )
        df_gate["lower CI(95%)"] = df_gate["ME"] + norm.ppf(q=0.05 / 2) * df_gate["sd"]
        df_gate["upper CI(95%)"] = (
            df_gate["ME"] + norm.ppf(q=1 - 0.05 / 2) * df_gate["sd"]
        )
        df_gate["lower uniform CI(95%)"] = df_gate["2.5 %"]
        df_gate["upper uniform CI(95%)"] = df_gate["97.5 %"]
        df_gate = df_gate[
            [
                "X",
                "ME",
                "sd",
                "lower CI(95%)",
                "upper CI(95%)",
                "lower uniform CI(95%)",
                "upper uniform CI(95%)",
            ]
        ].reset_index(drop=True)

    design_matrix = patsy.dmatrix("bs(x, df=5, degree=2)", {"x": df[X]})
    spline_basis = pd.DataFrame(design_matrix)
    spline_basis.to_csv("/Users/tianzhuqin/Documents/DML-Interaction/1.csv")
    cate = dml_model.cate(spline_basis)
    spline_basis.to_csv("/Users/tianzhuqin/Documents/DML-Interaction/2.csv")
    new_data = {"x": np.linspace(df[X].min(), df[X].max(), 50)}
    spline_grid = pd.DataFrame(
        patsy.build_design_matrices([design_matrix.design_info], new_data)[0]
    )
    spline_grid_np = spline_grid.to_numpy(df[X])

    df_cate = cate.confint(spline_grid, level=0.95, joint=True, n_rep_boot=2000)
    df_cate["X"] = new_data["x"]
    df_cate["ME"] = df_cate["effect"]
    df_cate["sd"] = np.sqrt(
        (np.dot(spline_grid_np, cate._blp_omega) * spline_grid_np).sum(axis=1)
    )
    df_cate["lower CI(95%)"] = df_cate["ME"] + norm.ppf(q=0.05 / 2) * df_cate["sd"]
    df_cate["upper CI(95%)"] = df_cate["ME"] + norm.ppf(q=1 - 0.05 / 2) * df_cate["sd"]
    df_cate["lower uniform CI(95%)"] = df_cate["2.5 %"]
    df_cate["upper uniform CI(95%)"] = df_cate["97.5 %"]
    df_cate = df_cate[
        [
            "X",
            "ME",
            "sd",
            "lower CI(95%)",
            "upper CI(95%)",
            "lower uniform CI(95%)",
            "upper uniform CI(95%)",
        ]
    ]
    return (
        df_cate.to_dict("list"),
        df_gate.to_dict("list"),
        params,
    )
