# Matthew Lueder
# M-state tutarial

#install.packages("mstate")
require("mstate")

data(aidssi)

si <- aidssi # Just a shorter name
head(si)

# 0 = censored, 1 = AIDS, 2 = SI appearance
table(si$status)


# The function trans.comprisk prepares a transition matrix for competing risks models.
# The first argument is the number of causes of failure.
# The names argument is a character vector of length three (the total number of states in the multi-state model 
#   including the failure-free state) may be given.
# The transition matrix has three states with state 1 being the failure-free state and the subsequent 
#   states representing the different causes of failure.
tmat <- trans.comprisk(2, names = c("event-free", "AIDS", "SI"))

# print out the transistion matrix
tmat

# Add new column that treats AIDS (1) as event of interest, and all other states as censored
si$stat1 <- as.numeric(si$status == 1)
# Add new column that treats SI (2) as event of interest, and all other states as censored
si$stat2 <- as.numeric(si$status == 2)

# To prepare data in long format, use msprep. (Converts from wide format to long format)
# Wide format: one subject per line, multiple columns indicating time and status for different states
# Long format: one line for each transition for which a subject is at risk
silong <- msprep(time = c(NA, "time", "time"), status = c(NA, "stat1", "stat2"), data = si, keep = "ccr5", trans = tmat)

# events: this function counts the number of observed transitions in the multi-state model and gives their percentages
# check whether the number of events from original data (si) corresponds with long data
events(silong)

# Add cause-specific covariates
# ccr5 is genotype (WW = wildtype, WM = one mutant allele)
silong <- expand.covs(silong, "ccr5")

# Show that niave KM estimators are biased
# Fit KM model to both events
# AIDS is shown as survival curve, SI shown as incidence curve
c1 <- coxph(Surv(time, status) ~ 1, data = silong, subset = (trans == 1), method = "breslow")
c2 <- coxph(Surv(time, status) ~ 1, data = silong, subset = (trans == 2), method = "breslow")
h1 <- survfit(c1)
h1 <- data.frame(time = h1$time, surv = h1$surv)
h2 <- survfit(c2)
h2 <- data.frame(time = h2$time, surv = h2$surv)
idx1 <- (h1$time<13) # this restricts the plot to the first 13 years
plot(c(0, h1$time[idx1], 13), c(1, h1$surv[idx1], min(h1$surv[idx1])), type="s", 
     xlim=c(0,13), ylim=c(0, 1), xlab="Years from HIV infection", ylab="Probability", lwd=2)
idx2 <- (h2$time<13)
lines(c(0, h2$time[idx2], 13), c(0, 1-h2$surv[idx2], max(1-h2$surv[idx2])), type="s", lwd=2)
text(8, 0.71, adj=0, "AIDS")
text(8, 0.32, adj=0, "SI")

# Now use mstate
# This function computes nonparametric cumulative incidence functions and 
# associated standard errors for each value of a group variable (two ways do the same thing shown)
ci <- Cuminc(time = si$time, status = si$status)
ci <- Cuminc(time = "time", status = "status", data = aidssi)
# ci = data frame containing the failure-free probabilities (Surv) and the cumulative incidence 
# functions with their standard error.
idx0 <- (ci$time < 13)
plot(c(0, ci$time[idx0], 13), c(1, 1 - ci$CI.1[idx0], min(1 - ci$CI.1[idx0])), type = "s", 
     xlim = c(0, 13), ylim = c(0, 1), xlab = "Years from HIV infection", ylab = "Probability", lwd = 2)
idx1 <- (h1$time < 13)
lines(c(0, h1$time[idx1], 13), c(1, h1$surv[idx1], min(h1$surv[idx1])), type = "s", lwd = 2, col = 8)
lines(c(0, ci$time[idx0], 13), c(0, ci$CI.2[idx0], max(ci$CI.2[idx0])), type = "s", lwd = 2)
idx2 <- (h2$time < 13)
lines(c(0, h2$time[idx2], 13), c(0, 1 - h2$surv[idx2], max(1 - h2$surv[idx2])), type = "s", lwd = 2, col = 8)
text(8, 0.77, adj = 0, "AIDS")
text(8, 0.275, adj = 0, "SI")

# Stacked plot: The cumulative incidence functions are stacked; the distances between two curves 
# represent the probabilities of the different events
idx0 <- (ci$time < 13)
plot(c(0, ci$time[idx0]), c(0, ci$CI.1[idx0]), type = "s", xlim = c(0, 13), ylim = c(0, 1), 
     xlab = "Years from HIV infection", ylab = "Probability", lwd = 2)
lines(c(0, ci$time[idx0]), c(0, ci$CI.1[idx0] + ci$CI.2[idx0]), type = "s", lwd = 2)
text(13, 0.5 * max(ci$CI.1[idx0]), adj = 1, "AIDS")
text(13, max(ci$CI.1[idx0]) + 0.5 * max(ci$CI.2[idx0]), adj = 1, "SI")
text(13, 0.5 + 0.5 * max(ci$CI.1[idx0]) + 0.5 * max(ci$CI.2[idx0]), adj = 1, "Event-free")

# * Now we use cox regression *
#   Using the original dataset, we can apply ordinary Cox regression for cause 1 (AIDS),
#   taking only the AIDS cases as events. This is done by specifying status==1 below (observations
#   with status=0 (true censorings) and status=2 (SI) are treated as censorings).
coxph(Surv(time, status == 1) ~ ccr5, data = si) # Treats SI events as censored
coxph(Surv(time, status == 2) ~ ccr5, data = si) # Treats AIDS events as censored




