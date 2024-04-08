## estimating marginal effects: binary treatment + continuous outcome

import patsy
import pandas as pd
import numpy as np
import doubleml as dml
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.neural_network import MLPClassifier, MLPRegressor

def marginal_effect_for_R(df, ml_method, Y, D, X, Z, 
                          trimming_threshold=0.01, n_folds=5,
                          n_estimators=500, 
                          solver='adam', max_iter=10000, alpha=1e-5, hidden_layer_sizes=(5, 3, 2), random_state=1):

    ### Initialization
    n_folds = int(n_folds)
    n_estimators = int(n_estimators)
    max_iter = int(max_iter)
    hidden_layer_sizes = tuple(int(1) for x in hidden_layer_sizes)
    random_state = int(random_state)

    # get labels of covariates
    if type(Z) is str:
        Z = [Z]
    covariates = Z.copy()
    covariates.append(X)

    df_base_index = df[D] == "Group.1"
    df.loc[df_base_index, D] = int(0)
    df.loc[~df_base_index, D] = int(1)
    df[D] = df[D].astype('int')
        
    ml_method_lower = ml_method.lower()
    ### First stage: estimate nuisance parameters
    # create random forest estimators for nuisance parameters g and m
    if ml_method_lower in {"randomforest", "random forest", "random_forest", "rf", "forest"}:
        ml_reg = RandomForestRegressor(n_estimators=n_estimators)
        ml_class = RandomForestClassifier(n_estimators=n_estimators)

    elif ml_method_lower in {"boosting", "gradient_boosting", "gradient boosting",
                          "hist_gradient_boosting", "hist gradient boosting",
                          "boost", "gradient_boost", "gradient boost", 
                          "hist_gradient_boost", "hist gradient boost",
                          "b", "gb", "hgb"}:                  
        ml_reg = HistGradientBoostingRegressor()
        ml_class = HistGradientBoostingClassifier()

    elif ml_method_lower in {"network", "neural_network", "neural network", "nn"}:
        ml_reg = MLPRegressor(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)
        ml_class = MLPClassifier(solver=solver, max_iter=max_iter, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state)

    else:
        print("'ml_method' should be one of 'rf', 'hgb', and 'nn'.")
        return None

    # create dml dataset object
    data_dml_base = dml.DoubleMLData(df[[X, Y, D, *Z]], y_col=Y, d_cols=D, x_cols=covariates)

    # create dml IRM object and start training
    dml_irm = dml.DoubleMLIRM(data_dml_base, ml_g=ml_reg, ml_m=ml_class, trimming_threshold=trimming_threshold, n_folds=n_folds)
    dml_irm.fit()

    ### Second stage: estimate marginal effects using orthogonal estimator
    # create b-spline basis
    design_matrix = patsy.dmatrix("bs(x, df=5, degree=2)", {"x":df[X]})
    spline_basis = pd.DataFrame(design_matrix)
    
    # estimate marginal effects
    cate = dml_irm.cate(spline_basis)
    
    # get confidence intervals
    new_data = {"x": np.linspace(df[X].min(), df[X].max(), 50)}
    spline_grid = pd.DataFrame(patsy.build_design_matrices([design_matrix.design_info], new_data)[0])
    df_cate = cate.confint(spline_grid, level=0.95, joint=True, n_rep_boot=2000)
    df_cate["X"] = new_data["x"]
    df_cate["ME"] = df_cate["effect"]
    df_cate["lower CI(95%)"] = df_cate["2.5 %"]
    df_cate["upper CI(95%)"] = df_cate["97.5 %"]
    df_cate = df_cate[["X", "ME", "lower CI(95%)", "upper CI(95%)"]]
    return df_cate.to_dict("list")