# interflex 1.3.1-1.3.2

1. Introduce the DML estimators

2. Fix bugs

# interflex 1.3.0

1. Add support for DML estimators.

2. Replace `fastplm` with `fixest` for estimating fixed effects models

3. Add analytical standard errors for the kernel estimator. 

# interflex 1.2.1

1. Allows users to conduct linear, binning and kernel estimation using glm functions including logit, probit, poisson and negative binomial model. 

2. Incorporate more flexible methods for uncertainty estimation, including simulation, delta method, and bootstrap.

3. Because the estimated marginal effects/treatment effects in glm rely on the value of covariates **Z**, we now add an option `Z.ref`, by which users can pass in the values of covariates as what they want. The mean of each covariate in the data will be set as the default value of `Z.ref`.

4. Incorporate the estimation of average treatment effects(if **D** is discrete) or average marginal effects(if **D** is continuous) of all observations.

5. When the outcome takes binary values (0 and 1), the program supports four criteria in cross-validation for the kernel method, including Mean Squared Error(MSE), Mean Absolute Error(MAE), Cross Entropy and Area Under the Curve(AUC). 

6. For the binning estimator, the program allows likelihood ratio test(lr-test) and wald test, providing a double insurance for verifying the linear interaction effect (LIE) assumption.

7. For all three estimators, the program allows a more flexible fully moderated model. When the option `full.moderate` is set to **TRUE**, all covariaes will interact with the moderater variable **X** in the same way as the treatment variable **D** does.


# interflex 1.1.1

1. Add a function **inter.test** to test the difference in treatment effects at two or three specific values in the moderator after the estimation.

# interflex 1.1.0

1. Accommodate discrete treatments variables (>2 treatment arms) interacting with continuous moderators; users can decide whether to draw the multiple marginal effects in one plot or in several by specifying the `pool` option.

2. Incorporate adaptive bandwidth search in cross validation; optimize cross-validation when  fixed effects are present.

3. Add an umbrella function **interflex** (S3 method) to incorporate the functionalities of both **inter.binning** and **inter.kernel** while **inter.binning** and **inter.kernel** are deprecated (for usage, see [old User's Guide](https://yiqingxu.org/packages/interflex/RGuide_1.0.9.html)).

4. Add a **predict** function, which allows users estimate and plot the expected value of Y given the values of the treatment, moderators, and controls. 

5. Allow users to test the difference in treatment effects at three different values in the moderate using the `diff.values` option

---

### Note:

Starting from v.1.1.0, we incorporate the interation terms between the moderator `X` and all covariates `Z` in our models to reduce biases induced by their correlations ([Blackwell and Olson 2019](https://www.mattblackwell.org/research/lasso-inters/)). Though not recommended, users can turn off this option by setting `full.moderate = FALSE`. In [22 published papers](https://yiqingxu.org/packages/interflex/SI_update.pdf) that we look into, we find the differences between the original model and the fully moderated model to be small. 

*Reference:* Matthew Blackwell and Michael Olson (2019). "Reducing Model Misspecification and Bias in the Estimation of Interactions." Mimeo, Harvard University.
