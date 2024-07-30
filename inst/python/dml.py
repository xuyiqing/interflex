## estimating marginal effects: continuous treatment + continuous outcome
import warnings

warnings.filterwarnings("ignore")

import json
import numpy as np
from sklearn.preprocessing import PolynomialFeatures, SplineTransformer
from sklearn.linear_model import LinearRegression, ElasticNet
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import (
    RandomForestClassifier,
    RandomForestRegressor,
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from econml.dml import NonParamDML
from econml.inference import BootstrapInference
from econml.grf import RegressionForest
from econml.sklearn_extensions.linear_model import (
    DebiasedLasso,
    StatsModelsLinearRegression,
)


def set_model(model, param, CV, discrete_outcome, param_grid, n_folds, scoring, n_jobs):
    model_lower = model.lower()
    if model_lower in {
        "default",
        "linear",
        "l",
        "d",
    }:
        model = LinearRegression()

    elif model_lower in {"regularization", "r"}:
        if CV:
            model = GridSearchCV(
                estimator=ElasticNet(),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = ElasticNet(**param)
            else:
                model = ElasticNet()

    elif model_lower in {
        "randomforest",
        "random forest",
        "random_forest",
        "rf",
        "forest",
    }:
        if CV:
            model = GridSearchCV(
                estimator=(
                    RandomForestClassifier()
                    if discrete_outcome
                    else RandomForestRegressor()
                ),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = (
                    RandomForestClassifier(**param)
                    if discrete_outcome
                    else RandomForestRegressor(**param)
                )
            else:
                model = (
                    RandomForestClassifier()
                    if discrete_outcome
                    else RandomForestRegressor()
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
        if CV:
            model = GridSearchCV(
                estimator=(
                    HistGradientBoostingClassifier()
                    if discrete_outcome
                    else HistGradientBoostingRegressor()
                ),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = (
                    HistGradientBoostingClassifier(**param)
                    if discrete_outcome
                    else HistGradientBoostingRegressor(**param)
                )
            else:
                model = (
                    HistGradientBoostingClassifier()
                    if discrete_outcome
                    else HistGradientBoostingRegressor()
                )

    elif model_lower in {"network", "neural_network", "neural network", "nn"}:
        if CV:
            model = GridSearchCV(
                estimator=MLPClassifier() if discrete_outcome else MLPRegressor(),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = (
                    MLPClassifier(**param)
                    if discrete_outcome
                    else MLPRegressor(**param)
                )
            else:
                model = MLPClassifier() if discrete_outcome else MLPRegressor()

    else:
        raise Exception(
            "'model_y' and 'model_t' should be one of 'linear', 'regularization', 'rf', 'hgb', and 'nn'."
        )

    return model


def set_model_final(model, param, CV, param_grid, n_folds, scoring, n_jobs):
    model_lower = model.lower()
    if model_lower in {
        "default",
        "linear",
        "l",
        "d",
    }:
        model = StatsModelsLinearRegression()

    elif model_lower in {"regularization", "r"}:
        if CV:
            model = GridSearchCV(
                estimator=DebiasedLasso(),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = DebiasedLasso(**param)
            else:
                model = DebiasedLasso()

    elif model_lower in {
        "randomforest",
        "random forest",
        "random_forest",
        "rf",
        "forest",
    }:
        if CV:
            model = GridSearchCV(
                estimator=RegressionForest(),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = RegressionForest(**param)
            else:
                model = RegressionForest()

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
        if CV:
            model = GridSearchCV(
                estimator=HistGradientBoostingRegressor(),
                param_grid=param_grid,
                cv=n_folds,
                n_jobs=n_jobs,
                scoring=scoring,
            )
        else:
            if param != {}:
                model = HistGradientBoostingRegressor(**param)
            else:
                model = HistGradientBoostingRegressor()

    # elif model_lower in {"network", "neural_network", "neural network", "nn"}:
    #     if CV:
    #         model = GridSearchCV(
    #             estimator=MLPRegressor(),
    #             param_grid=param_grid,
    #             cv=n_folds,
    #             n_jobs=n_jobs,
    #             scoring=scoring,
    #         )
    #     else:
    #         if param != {}:
    #             model = MLPRegressor(**param)
    #         else:
    #             model = MLPRegressor()

    else:
        raise Exception(
            "'model_final' should be one of 'linear', 'regularization', 'rf', and 'hgb'."
        )

    return model


def marginal_effect_for_treatment(
    # data
    df,
    Y,
    D,
    X,
    Z,
    d_ref=1,
    # model
    FSCF_n_folds=2,
    model_y="rf",
    param_y={},
    CV_y=False,
    param_grid_y={},
    n_folds_y=10,
    scoring_y="neg_mean_squared_error",
    model_t="rf",
    param_t={},
    CV_t=False,
    param_grid_t={},
    n_folds_t=10,
    scoring_t="neg_mean_squared_error",
    model_final="linear",
    param_final={},
    CV_final=False,
    param_grid_final={},
    n_folds_final=10,
    scoring_final="neg_mean_squared_error",
    featurizer_model_final=None,
    featurizer_param_final={},
    n_jobs=-1,
):

    df[D] = df[D].astype(float)

    ### Initialization
    FSCF_n_folds = int(FSCF_n_folds)
    n_folds_y = int(n_folds_y)
    n_folds_t = int(n_folds_t)
    n_folds_final = int(n_folds_final)
    n_jobs = int(n_jobs)

    if type(Z) is str:
        Z = [Z]

    if len(df[Y].unique()) > 5:
        discrete_outcome = False
    else:
        discrete_outcome = True

    if len(df[D].unique()) > 5:
        discrete_treatment = False
    else:
        discrete_treatment = True

    featurizer_model_final_lower = featurizer_model_final.lower()
    if featurizer_model_final_lower in {"p", "poly", "polynomial"}:
        if featurizer_param_final == {}:
            featurizer_final = PolynomialFeatures()
        else:
            featurizer_final = PolynomialFeatures(**featurizer_param_final)
    elif featurizer_model_final_lower in {"s", "spline"}:
        if featurizer_param_final == {}:
            featurizer_final = SplineTransformer()
        else:
            featurizer_final = SplineTransformer(**featurizer_param_final)
    else:
        raise Exception(
            "'featurizer_model_final' should be one of 'poly', and 'spline'."
        )

    model_y_set = set_model(
        model_y,
        param_y,
        CV_y,
        discrete_outcome,
        param_grid_y,
        n_folds_y,
        scoring_y,
        n_jobs,
    )
    model_t_set = set_model(
        model_t,
        param_t,
        CV_t,
        discrete_treatment,
        param_grid_t,
        n_folds_t,
        scoring_t,
        n_jobs,
    )
    model_final_set = set_model_final(
        model_final,
        param_final,
        CV_final,
        param_grid_final,
        n_folds_final,
        scoring_final,
        n_jobs,
    )

    est = NonParamDML(
        model_y=model_y_set,
        model_t=model_t_set,
        model_final=model_final_set,
        discrete_outcome=discrete_outcome,
        discrete_treatment=discrete_treatment,
        cv=FSCF_n_folds,
        featurizer=featurizer_final,
    )

    est.fit(df[[Y]].values.ravel(), df[[D]].values.ravel(), X=df[[X]], W=df[[*Z]])

    if CV_final and (model_final not in {"default", "linear", "l", "d"}):
        model_final_set_after_cv = est.model_cate.best_estimator_
    else:
        model_final_set_after_cv = est.model_cate

    est = NonParamDML(
        model_y=est.models_y[0][0],
        model_t=est.models_t[0][0],
        model_final=model_final_set_after_cv,
        discrete_outcome=discrete_outcome,
        discrete_treatment=discrete_treatment,
        cv=FSCF_n_folds,
        featurizer=featurizer_final,
    )

    if model_final in {
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
        "network",
        "neural_network",
        "neural network",
        "nn",
    }:
        est.fit(
            df[[Y]].values.ravel(),
            df[[D]].values.ravel(),
            X=df[[X]],
            W=df[[*Z]],
            inference=BootstrapInference(n_bootstrap_samples=100, n_jobs=n_jobs),
        )
    else:
        est.fit(df[[Y]].values.ravel(), df[[D]].values.ravel(), X=df[[X]], W=df[[*Z]])

    length = 50
    new_data = {"x": np.linspace(df[X].min(), df[X].max(), length)}
    df_cate = est.effect_inference(
        new_data["x"].reshape(length, 1), T0=0, T1=d_ref
    ).summary_frame(alpha=0.05, value=0, decimals=99)

    df_cate["X"] = new_data["x"]
    df_cate["ME"] = df_cate["point_estimate"]
    df_cate["sd"] = df_cate["stderr"]
    df_cate["lower CI(95%)"] = df_cate["ci_lower"]
    df_cate["upper CI(95%)"] = df_cate["ci_upper"]
    df_cate = df_cate[["X", "ME", "sd", "lower CI(95%)", "upper CI(95%)"]]

    return (
        df_cate.to_dict("list"),
        est.models_y[0][0].get_params(),
        est.models_t[0][0].get_params(),
        est.model_cate.get_params(),
    )
