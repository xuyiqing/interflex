## estimating marginal effects: continuous treatment + continuous outcome
import warnings
warnings.filterwarnings("ignore")

import json
import numpy as np
from econml.dml import DML, LinearDML, SparseLinearDML, CausalForestDML
from sklearn.linear_model import Lasso, LassoCV
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.model_selection import GridSearchCV

def marginal_effect_for_treatment(df, ml_method, Y, D, X, Z, d_ref = 1,
                                  n_estimators=500, solver='adam', max_iter=10000, alpha=1e-5, hidden_layer_sizes=(5, 3, 2), random_state=1,
                                  dml_method='default',
                                  poly_degree=3, lasso_alpha=0.0001,
                                  casual_forest_criterion="mse",
                                  casual_forest_n_estimators=1000,
                                  casual_forest_in_impurity_decrease=0.001, 
                                  CV_y=False, param_grid_y=None, n_folds_y=10, scoring_y='neg_mean_squared_error',
                                  CV_t=False, param_grid_t=None, n_folds_t=10, scoring_t='neg_mean_squared_error',
                                  CV_f=False, param_grid_f=None, n_folds_f=10, scoring_f='neg_mean_squared_error',
                                  n_jobs=-1):
    if param_grid_y is not None:
        try:
            param_grid_y = json.loads(param_grid_y)
        except json.JSONDecodeError:
            raise Exception("'param_grid_y' should be a dictionary-like string.")

    if param_grid_t is not None:
        try:
            param_grid_t = json.loads(param_grid_t)
        except json.JSONDecodeError:
            raise Exception("'param_grid_t' should be a dictionary-like string.")
        
    if param_grid_f is not None:
        try:
            param_grid_f = json.loads(param_grid_f)
        except json.JSONDecodeError:
            raise Exception("'param_grid_f' should be a dictionary-like string.")

    ### Initialization
    n_folds_y = int(n_folds_y)
    n_folds_t = int(n_folds_t)
    n_folds_f = int(n_folds_f)
    n_jobs = int(n_jobs)
    n_estimators = int(n_estimators)
    max_iter = int(max_iter)
    hidden_layer_sizes = tuple(int(1) for x in hidden_layer_sizes)
    random_state = int(random_state)
    poly_degree = int(poly_degree)
    casual_forest_n_estimators = int(casual_forest_n_estimators)
    if type(Z) is str:
        Z = [Z]
    covariates = Z.copy()
    covariates.append(X)

    if len(df[Y].unique()) > 5:
        discrete_outcome = False
    else:
        discrete_outcome = True

    if len(df[D].unique()) > 5:
        discrete_treatment = False
    else:
        discrete_treatment = True
    
    ml_method_lower = ml_method.lower()
    dml_method_lower = dml_method.lower()

    ### Choose ML method 
    if ml_method_lower in {"randomforest", "random forest", "random_forest", "rf", "forest"}:
        if CV_y:
            stage_model_y = lambda: GridSearchCV(estimator=RandomForestClassifier() if discrete_outcome else RandomForestRegressor(), 
                                                 param_grid=param_grid_y, cv=n_folds_y, n_jobs=n_jobs, scoring=scoring_y)
            model_y = stage_model_y().fit(df[[X]], df[[Y]].values.ravel()).best_estimator_
        else:
            model_y = RandomForestClassifier(n_estimators=n_estimators) if discrete_outcome else RandomForestRegressor(n_estimators=n_estimators)
        
        if CV_t:
            stage_model_t = lambda: GridSearchCV(estimator=RandomForestClassifier() if discrete_treatment else RandomForestRegressor(),
                                                 param_grid=param_grid_t, cv=n_folds_t, n_jobs=n_jobs, scoring=scoring_t)
            model_t = stage_model_t().fit(df[[X]], df[[D]].values.ravel()).best_estimator_
        else:
            model_t = RandomForestClassifier(n_estimators=n_estimators) if discrete_treatment else RandomForestRegressor()

    elif ml_method_lower in {"boosting", "gradient_boosting", "gradient boosting",
                        "hist_gradient_boosting", "hist gradient boosting",
                        "boost", "gradient_boost", "gradient boost", 
                        "hist_gradient_boost", "hist gradient boost",
                        "b", "gb", "hgb"}:
        if CV_y:
            stage_model_y = lambda: GridSearchCV(estimator=HistGradientBoostingClassifier() if discrete_outcome else HistGradientBoostingRegressor(),
                                                 param_grid=param_grid_y, cv=n_folds_y, n_jobs=n_jobs, scoring=scoring_y)
            model_y = stage_model_y().fit(df[[X]], df[[Y]].values.ravel()).best_estimator_
        else:
            model_y = HistGradientBoostingClassifier() if discrete_outcome else HistGradientBoostingRegressor()
        
        if CV_t:
            stage_model_t = lambda: GridSearchCV(estimator=HistGradientBoostingClassifier() if discrete_treatment else HistGradientBoostingRegressor(),
                                                 param_grid=param_grid_t, cv=n_folds_t, n_jobs=n_jobs, scoring=scoring_t)
            model_t = stage_model_t().fit(df[[X]], df[[D]].values.ravel()).best_estimator_
        else:
            model_t = HistGradientBoostingClassifier() if discrete_treatment else HistGradientBoostingRegressor()

    elif ml_method_lower in {"network", "neural_network", "neural network", "nn"}:
        if CV_y:
            stage_model_y = lambda: GridSearchCV(estimator=MLPClassifier() if discrete_outcome else MLPRegressor(),
                                                 param_grid=param_grid_y, cv=n_folds_y, n_jobs=n_jobs, scoring=scoring_y)
        else:
            model_y = MLPClassifier(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state) if discrete_outcome else MLPRegressor(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        
        if CV_t:
            stage_model_t = lambda: GridSearchCV(estimator=MLPClassifier() if discrete_treatment else MLPRegressor(),
                                                 param_grid=param_grid_t, cv=n_folds_t, n_jobs=n_jobs, scoring=scoring_t)
        else:
            model_t = MLPClassifier(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state) if discrete_treatment else MLPRegressor(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)

    else:
        raise Exception("'ml_method' should be one of 'rf', 'hgb', and 'nn'.")

    ### Choose DML method
    if dml_method_lower == "default"[:len(dml_method_lower)]:
        est = LinearDML(model_y=model_y, model_t=model_t, random_state=random_state, discrete_outcome=discrete_outcome, discrete_treatment=discrete_treatment)

    elif dml_method_lower == "polynomial"[:len(dml_method_lower)]:
        est = SparseLinearDML(model_y=model_y, model_t=model_t,
                            featurizer=PolynomialFeatures(degree=poly_degree),
                            random_state=random_state, discrete_outcome=discrete_outcome, discrete_treatment=discrete_treatment)
        
    elif dml_method_lower == "regularization"[:len(dml_method_lower)]:
        if CV_f:
            est = DML(model_y=model_y, model_t=model_t,
                    model_final=LassoCV(fit_intercept=False, cv=n_folds_f, n_jobs=n_jobs, alphas=param_grid_f["lasso_alpha"]),
                    featurizer=PolynomialFeatures(degree=poly_degree),
                    random_state=random_state, discrete_outcome=discrete_outcome, discrete_treatment=discrete_treatment)
        else:
            est = DML(model_y=model_y, model_t=model_t,
                    model_final=Lasso(alpha=lasso_alpha, fit_intercept=False),
                    featurizer=PolynomialFeatures(degree=poly_degree),
                    random_state=random_state, discrete_outcome=discrete_outcome, discrete_treatment=discrete_treatment)
        
    elif dml_method_lower == "non-parametric"[:len(dml_method_lower)]:
        est = CausalForestDML(model_y=model_y, model_t=model_t,
                            criterion=casual_forest_criterion, 
                            n_estimators=casual_forest_n_estimators,
                            min_impurity_decrease=casual_forest_in_impurity_decrease,
                            random_state=random_state, discrete_outcome=discrete_outcome, discrete_treatment=discrete_treatment)
        
        if CV_f:
            est.tune(df[[Y]], df[[D]], X=df[[X]], W=df[[*Z]], params=param_grid_f)

    else:
        raise Exception("'dml_method' should be one of 'default', 'polynomial', 'regularization', and 'non-parametric'.")
        
    est.fit(df[[Y]], df[[D]], X=df[[X]], W=df[[*Z]])
    
    ### Get results
    length = 50
    new_data = {"x": np.linspace(df[X].min(), df[X].max(), length)}
    df_cate = est.effect_inference(new_data["x"].reshape(length, 1), T0=0, T1=d_ref).summary_frame(alpha=0.05, value=0, decimals=99)

    df_cate["X"] = new_data["x"]
    df_cate["ME"] = df_cate["point_estimate"]
    if dml_method_lower == "regularization":
        df_cate["sd"] = 0
        df_cate["lower CI(95%)"] = df_cate["point_estimate"]
        df_cate["upper CI(95%)"] = df_cate["point_estimate"]
    else:
        df_cate["sd"] = df_cate["stderr"]
        df_cate["lower CI(95%)"] = df_cate["ci_lower"]
        df_cate["upper CI(95%)"] = df_cate["ci_upper"]
    df_cate = df_cate[["X", "ME", "sd", "lower CI(95%)", "upper CI(95%)"]]
    return df_cate.to_dict("list")