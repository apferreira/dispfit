logsech.function <- function (data, chi.res.hist, ks.res.hist) {
log.dist.logsech <- function (par, r) {
  a <- par[1] ## location, median
  b <- par[2] ## scale
  flogsech <- (1 / ((pi^2) * b * (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))
  -sum(log(flogsech)) ##
}
dist.logsech <- function (r, a, b) {
  flogsech <- 2*pi*r * (1 / ((pi^2) * b * (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))
}
# initial values estimation # from Halley and Inchausti 2002
location <- median(log(data))
scale <- IQR(log(data))/(2*-log(tan(pi/8)))
# optimization procedure
dist.logsech.opt <- optim (par = c(location, scale), ## valor inicial para o "a"
                           fn = log.dist.logsech, ## função a minimizar
                           r = data, ## dados
                           control = list(maxit = 10000),
                           method = "Nelder-Mead")
#lower = c(0.00001, 0.00001),

# output values
# AIC
aic.logsech <- 2*length(dist.logsech.opt$par) + 2 * dist.logsech.opt$value
# AICc
aicc.logsech <- aic.logsech + (2 * length(dist.logsech.opt$par)^2 + 2 * length(dist.logsech.opt$par))/(length(data) - length(dist.logsech.opt$par) - 1 )
# BIC
bic.logsech <-  2 * dist.logsech.opt$value + length(dist.logsech.opt$par)*log(length(data))
# Chi-squared
chi.expected.values.logsech <- dist.logsech(chi.res.hist$mids, dist.logsech.opt$par[1], dist.logsech.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
chi.squared.statistic.logsech <- sum((chi.res.hist$counts - chi.expected.values.logsech)^2 / chi.expected.values.logsech)
chi.squared.pvalue.logsech <- 1-pchisq(chi.squared.statistic.logsech, length(chi.res.hist$counts)-3)
# Kolmogorov-Smirnov
ks.expected.values.logsech <- dist.logsech(ks.res.hist$mids, dist.logsech.opt$par[1], dist.logsech.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
simul.logsech <- c()
for (i in seq_along(ks.res.hist$mids)) {
  simul.logsech <- c(simul.logsech, rep(ks.res.hist$mids[i], round(ks.expected.values.logsech[i], 0)))
}
ks.logsech <- ks.test(data, simul.logsech)
g.max.logsech <- as.numeric(ks.logsech$statistic)
KS.logsech <- as.numeric(ks.logsech$p.value)

# cumulative.expected.values.logsech <- c(expected.values.logsech[1])
# for (i in 1+seq_along(expected.values.logsech)) {
#   cumulative.expected.values.logsech[i] <- cumulative.expected.values.logsech[i-1] + expected.values.logsech[i]
# }
# cumulative.expected.values.logsech <- cumulative.expected.values.logsech/sum(expected.values.logsech)
# cumulative.expected.values.logsech <- cumulative.expected.values.logsech[!is.na(cumulative.expected.values.logsech)]
# g.max.logsech <- max(abs(cumulative.data - cumulative.expected.values.logsech))
# if (g.max.logsech < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
#   KS.logsech <- "Accept"
# } else {KS.logsech <- "Reject"}
# parameter estimate
par.1.logsech <- dist.logsech.opt$par[1]
par.2.logsech <- dist.logsech.opt$par[2]
# parameter estimate standard error
par.1.se.logsech <- sqrt(diag(solve(numDeriv::hessian(log.dist.logsech, x=dist.logsech.opt$par, r=data))))[1]
par.2.se.logsech <- sqrt(diag(solve(numDeriv::hessian(log.dist.logsech, x=dist.logsech.opt$par, r=data))))[2]
# mean dispersal distance
if (dist.logsech.opt$par[2] < 1) {
  mean.logsech <- 2*pi* (integrate(function(r) (r^2 / ((pi^2) * dist.logsech.opt$par[2] * (r^2))) / (((r / dist.logsech.opt$par[1])^(1 / dist.logsech.opt$par[2])) + ((r / dist.logsech.opt$par[1]) ^ -(1 / dist.logsech.opt$par[2]))),
                                   lower = 0,
                                   upper = Inf))$value
} else {
  mean.logsech <- "Infinite Value"
}
if (dist.logsech.opt$par[2] < 1) {
  mean.stderr.logsech <- "in progress"
} else {
  mean.stderr.logsech <- "Infinite Value"
}
# variance
variance.logsech <- "in progress"
variance.stderr.logsech <-"in progress"
# skewness
skewness.logsech <- "in progress"
skewness.stderr.logsech <- "in progress"
# kurtosis
kurtosis.logsech <- "in progress"
kurtosis.stderr.logsech <- "in progress"
# output
res <- data.frame(aic.logsech, aicc.logsech, bic.logsech,
                             chi.squared.statistic.logsech, chi.squared.pvalue.logsech,g.max.logsech, KS.logsech,
                             par.1.logsech, par.1.se.logsech, par.2.logsech, par.2.se.logsech,
                             mean.logsech, mean.stderr.logsech, variance.logsech, variance.stderr.logsech,
                             skewness.logsech, skewness.stderr.logsech, kurtosis.logsech, kurtosis.stderr.logsech)
logsech.values <- list("opt" = dist.logsech.opt, "res" = res)
}
