## estimating marginal effects: continuous treatment + continuous outcome
import warnings
warnings.filterwarnings("ignore")

import numpy as np
from econml.dml import DML, LinearDML, SparseLinearDML, CausalForestDML
from sklearn.linear_model import Lasso
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor

def marginal_effect_for_continuous_treatment(df, ml_method, Y, D, X, Z, d_ref = 1,
                          n_estimators=500, solver='adam', max_iter=10000, alpha=1e-5, hidden_layer_sizes=(5, 3, 2), random_state=1,
                          dml_method='default',
                          poly_degree=3, lasso_alpha=0.1,
                          casual_forest_criterion="mse",
                          casual_forest_n_estimators=1000,
                          casual_forest_in_impurity_decrease=0.001):

    ### Initialization
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

    ## Below are methods for binary treatment 
    # df_base_index = df[D] == "Group.1"
    # df.loc[df_base_index, D] = int(0)
    # df.loc[~df_base_index, D] = int(1)
    # df[D] = df[D].astype('int')

    ### Choose ML method 
    ml_method_lower = ml_method.lower()
    if ml_method_lower in {"randomforest", "random forest", "random_forest", "rf", "forest"}:
        model_y = RandomForestClassifier(n_estimators=n_estimators) if discrete_outcome else RandomForestRegressor(n_estimators=n_estimators)
        model_t = RandomForestRegressor(n_estimators=n_estimators)

    elif ml_method_lower in {"boosting", "gradient_boosting", "gradient boosting",
                          "hist_gradient_boosting", "hist gradient boosting",
                          "boost", "gradient_boost", "gradient boost", 
                          "hist_gradient_boost", "hist gradient boost",
                          "b", "gb", "hgb"}:                  
        model_y = HistGradientBoostingClassifier() if discrete_outcome else HistGradientBoostingRegressor()
        model_t = HistGradientBoostingRegressor()

    elif ml_method_lower in {"network", "neural_network", "neural network", "nn"}:
        model_y = MLPClassifier(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state) if discrete_outcome else MLPRegressor(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        model_t = MLPRegressor(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)

    else:
        print("'ml_method' should be one of 'rf', 'hgb', and 'nn'.")
        return None

    ### Choose DML method
    dml_method_lower = dml_method.lower()
    if dml_method_lower == "default"[:len(dml_method_lower)]:
        est = LinearDML(model_y=model_y, model_t=model_t, random_state=random_state, discrete_outcome=discrete_outcome)

    elif dml_method_lower == "polynomial"[:len(dml_method_lower)]:
        est = SparseLinearDML(model_y=model_y, model_t=model_t,
                            featurizer=PolynomialFeatures(degree=poly_degree),
                            random_state=random_state, discrete_outcome=discrete_outcome)
        
    elif dml_method_lower == "regularization"[:len(dml_method_lower)]:
        est = DML(model_y=model_y, model_t=model_t,
                  model_final=Lasso(alpha=lasso_alpha, fit_intercept=False),
                  featurizer=PolynomialFeatures(degree=poly_degree),
                  random_state=random_state, discrete_outcome=discrete_outcome)
        
    elif dml_method_lower == "non-parametric"[:len(dml_method_lower)]:
        est = CausalForestDML(model_y=model_y, model_t=model_t,
                              criterion=casual_forest_criterion, 
                              n_estimators=casual_forest_n_estimators,
                              min_impurity_decrease=casual_forest_in_impurity_decrease,
                              random_state=random_state, discrete_outcome=discrete_outcome)
        est.tune(df[[Y]], df[[D]], X=df[[X]], W=df[[*Z]])

    else:
        print("'dml_method' should be one of 'default', 'polynomial', 'regularization', and 'non-parametric'.")
        return None

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