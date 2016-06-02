# Low-dim survival data (train1):
# - plot Kaplan-Meier curve
# - fit cox ph model & interpret
# - use pec() to plot prediction error curve for the cox ph fit. Do this on the training data and using 10fold Bootstrap cross validation.
# 
# High-dim survival data (train2):
# - plot Kaplan-Meier curve
# - fit lasso model via cross validation (variables X1 and X2 are mandatory variables, should not be penalised)
# - make sure your result is reproducible
# - interpret
# 
# competing risks data (train3) :
# - plot Kaplan-Meier curve for each cause
# - plot the cumulative cause-specific hazards (Nelson-Aalen estimator) for the competing risks model (can use mvna() here)
# - plot the cumulative incidence functions / empirical transition matrix (Aalen-Johansen estimator) of the transition probability matrix of the competing risks model (can use etm() here)
# - (Low-dim): fit cause-specific cox ph models using the first 10 variables. Interpret.
# - (High-dim:) fit cause-specific lasso models using all variables (again via cross validation, variables X1 and X2 are mandatory variables, should not be penalised)

setwd("C:\\Users\\Matthew\\Desktop\\DKFZ\\Competing_Risks")
load("simdat1.Rdata")

require("survival")

# For reproducability
set.seed(248)

# ** Low dim survival analysis - train1 **
#   Time = train1[,11]
#   Status = train[,12] - 0 = censored, 1 = event

# Fit Kaplan Meier curve and plot it.
#   + = point where an individual became censored
#   Vertical drops = point where an event occured
#   Dotted lines = confidence intervals
simpKM <- survfit(Surv(time, status) ~ 1, type="kaplan-meier", data = train1)
plot(simpKM, main="Kaplan Meier Plot", xlab="t", ylab="Survival")

summary(simpKM)

# Fit a cox PH model
#   . = all variables not already mentioned in formula
coxPHModel = coxph(Surv(time, status) ~ ., data = train1)
# Interpretation:
#   exp(coef) = Hazard ratio: Hazard in treatment arm / Hazard in control arm
#     - A hazard ratio greater than one indicates that as the parameter associated with the coefficient increases,
#       the event hazard increases (length of survival decreases)
#     - A hazard ratio less than one is negatively associated with event probability and positively associated with survival
#   The only predictors found to be signifcantly correlated with survival (p < 0.05) are X1 and X2.
#   A patient with higher X1 is more likely to have an event earlier on
#   A patient with higher X2 is likely to survive longer (or have an event later/no event)

# Package for creeating prediction error curves
#install.packages("pec")
require("pec")

# Prediction error curves: Brier score plotted over time
# Brier score: Weighted average of the squared distance between observed survival status and predicted survival probability
# The weights used while calculating the Brier score correspond to probability of not being censored.
# Boostrapping can be automatically applied. 
# Reference line is kaplan meier (it does not use information about covariates to model prediction)
# We would like to see less error in our model than in the reference - shows it is better than more basic KM model
# The formula argument allows specification of variables that effect censoring.
#   Right now we have it set so it assumes all covariates are independent of censoring (~1)
# Leaving out bootstrapping would test model on same data it trained on. This would give a artificially good looking model.
predictionError = pec(object = coxPHModel, 
                      formula = Surv(time, status) ~ 1, data = train1,
                      splitMethod = "BootCv", B = 10)
plot(predictionError)
summary(predictionError)

# ** High-dimensional Data - train2 **
highDimKM <- survfit(Surv(time, status) ~ 1, type="kaplan-meier", data = train2)
plot(highDimKM, main="Kaplan Meier Plot", xlab="t", ylab="Survival")

# Use penalized package for lasso - See glmnet
#install.packages("penalized")
require("penalized")
vignette("penalized")

# First argument is response variable (Surv object)
# Second argument are the covariates to be penalized (can be specified as formula or matrix)
# Third argument are the covariates which you do not want to have penalized (also as formuala/matrix)
# This function uses lambda1 for lasso, and lambda2 for L2 penalty 
exLasso = penalized(Surv(time, status), ~.-X1-X2-time-status, ~X1+X2, data=train2, lambda1=1, model="cox")
# Same as:
exLasso = penalized(Surv(time, status), train2[,3:300], train2[,1:2], data=train2, lambda1=1, model="cox")

# See Coefficients
coefficients(exLasso, "nonzero")
coefficients(exLasso, "all")
coefficients(exLasso, "unpenalized")

# Use cross-validation to choose lambda
lassoFit = optL1(Surv(time, status), train2[,3:300], train2[,1:2], data=train2, model="cox")
# Print lambda
lassoFit$lambda
# Print Coefficients
coefficients(lassoFit$fullfit)
# 6 Coefficients are non-zero: X1, X2, X4, X5, X115, X203
# The coefficients for X5 and X115 are very small and with a slightly higher value of lambda, 
# they would drop out of the model. 
 
# Create a plot showing likelihood for different values of lambda
lassoProfile = profL1(Surv(time, status), train2[,3:300], train2[,1:2], data=train2, model="cox", fold=10)
plot(lassoProfile$lambda, lassoProfile$cvl, type="l", log="x")
# The plotpath function can again be used to visualize the effect of the tuning parameter on the regression coefficients
plotpath(lassoProfile$fullfit, log="x")

# How do you make a prediction error curve with a lasso model???

# ** Competing Risks **
# 1) Cause-specific Kaplan-meier curves
cause1_km = survfit(Surv(time, status) ~ 1, type="kaplan-meier", data = subset(train3, (status == 1 | status == 0)))
plot(cause1_km, main="Cause 1 Only, KM Plot", xlab="t", ylab="Survival")

cause2_km = survfit(Surv(time, status == 2) ~ 1, type="kaplan-meier", data = subset(train3, (status == 2 | status == 0)))
plot(cause2_km, main="Cause 2 Only, KM Plot", xlab="t", ylab="Survival")

# 2) Nelson-Aalen estimator for the competing risks model
# Package for Nelson-Aalen estimator
#install.packages("mvna")
require("mvna")

# Create a transition matrix
tmat = matrix( c(F,F, F, T, F, F, T, F, F), nrow=3, ncol=3) 

# Call mvna function 
NA_est = mvna(train3, c("0", "1", "2"), tmat, "cens")
summary(NA_est)

# plot
plot(NA_est)

# 3) Cumulative incidence functions / empirical transition matrix (Aalen-Johansen estimator)
#install.packages("etm")
require("etm")

AJ_est = etm(train3, c("0", "1", "2"), tmat, "cens", 0)
summary(AJ_est)

# Plot
require("lattice")
xyplot(AJ_est, tr.choice=c("0 0", "0 1", "0 2"), layout=c(1, 3), strip=strip.custom(bg="white", 
      factor.levels= c("0 to 0", "0 to 1", "0 to 2")))

# Can also be done with the prodlim package
#install.packages("prodlim")
require("prodlim")

crFit <- prodlim(Hist(time,status)~1,data=train3)
summary(crFit)
plot(crFit, cause=1)
plot(crFit,cause=2)

# The prodlim return object is implemented in pec
predictionError = pec(object = crFit, 
                      formula = Hist(time, status) ~ 1, data = train3,
                      splitMethod = "BootCv", B = 10)
plot(predictionError)
# Q: Should we use Hist for left side of formula?

# 4) Fit cause-specific cox ph models using the first 10 variables. Interpret.
#install.packages("mstate")
require("mstate")
#install.packages("riskRegression")
require("riskRegression")

# Use Hist rather than Surv: Hist provides the suitable extensions for dealing with right censored and 
#   interval censored data from competing risks and other multi state models.
# CSC: Cause-specific Cox proportional hazard regression
#   First argument is a list of formulae, one for each cause, each specifying a cause-specific Cox regression model
#   "survtype" is either 'survival' or 'hazard'
#     If "hazard" fit cause-specific Cox regression models for all causes.
#     If "survival" fit one cause-specific Cox regression model for the cause of interest and also a 
#       Cox regression model for event-free survival.
f10Form <- as.formula(paste("Hist(time, status) ~", paste(colnames(train3[,6:15]), collapse=" + "))) # Build formula
cscFit1 <- CSC(list(f10Form, f10Form), data=train3, cause="1", survtype="hazard")
cscFit1

# Compare to regular coxph function w/ competing cause being treated as censored
f10FormReg <- as.formula(paste("Surv(time, status==1) ~", paste(colnames(train3[,6:15]), collapse=" + "))) 
cphFitCH = coxph(f10FormReg, data = train3)
cphFitCH

 # 5) fit cause-specific lasso models using all variables (again via cross validation, 
#    variables X1 and X2 are mandatory variables, should not be penalised).

# First argument is response variable (Surv object)
# Second argument are the covariates to be penalized (can be specified as formula or matrix)
# Third argument are the covariates which you do not want to have penalized (also as formuala/matrix)
# This function uses lambda1 for lasso, and lambda2 for L2 penalty 

# X1 and X2 do not exist, treating c1 and c2 as manditory variables instead
lassoFit = optL1(Surv(time, status==1), train3[,8:411], train3[,6:7], data=train3, model="cox")
# Print lambda
lassoFit$lambda
# Print Coefficients
coefficients(lassoFit$fullfit)

lassoFit2 = optL1(Surv(time, status==2), train3[,8:411], train3[,6:7], data=train3, model="cox")
# Print lambda
lassoFit2$lambda
# Print Coefficients
coefficients(lassoFit2$fullfit)


