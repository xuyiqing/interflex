model.kernel <- interflex(estimator = "kernel",data = s2,
                          Y = 'Y', D = 'D', X = 'X', Z = 'Z1', bw = 1.8354,
                          full.moderate = T)

param_grid <- '{
  "max_depth": [3],
  "n_estimators": [10, 30, 50, 100, 200, 400, 600, 800, 1000],
  "max_features": [1,2,3]
}'

s7.DML.rf <- interflex(estimator='DML', data = s7, ml_method="rf", CV=TRUE, param_grid=param_grid, dml_method = "regularization",
                       Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), treat.type = "continuous")
plot(s7.DML.rf)

s7.DML.nn <- interflex(estimator="DML", data = s7, ml_method="rf",
                       Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
                       treat.type = "continuous")
plot(s7.DML.nn)
