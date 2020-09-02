cauchy.function <- function (data, chi.res.hist, ks.res.hist) {
log.dist.cauchy <- function (par, r) {
  a <- par[1]
  b <- par[2]
  fcauchy <- (1/(2*pi*r))*(1 / (pi*b)) / (1 + ((r-a)/b)^2)
  -sum(log(fcauchy)) ##
}
dist.cauchy <- function (r, a, b) {
  fcauchy <- 2*pi*r * (1/(2*pi*r))*(1 / (pi*b)) / (1 + ((r-a)/b)^2)
}
# initial values estimation # from Halley and Inchausti 2002
location <- median(data)
scale <- IQR(data)/2
# optimization procedure
dist.cauchy.opt <- optim (par = c(location, scale), ## valor inicial para o "a"
                          fn = log.dist.cauchy, ## função a minimizar
                          r = data, ## dados
                          #control = list(maxit = 10000), ## limite superior do parametro
                          method = "Nelder-Mead"
                          # lower = c(0.00001, 0.00001)
)
kernel.fit$cauchy <- dist.cauchy.opt
# output values
# AIC
aic.cauchy <- 2*length(dist.cauchy.opt$par) + 2 * dist.cauchy.opt$value
# AICc
aicc.cauchy <- aic.cauchy + (2 * length(dist.cauchy.opt$par)^2 + 2 * length(dist.cauchy.opt$par))/(length(data) - length(dist.cauchy.opt$par) - 1 )
# BIC
bic.cauchy <-  2 * dist.cauchy.opt$value + length(dist.cauchy.opt$par)*log(length(data))
# Chi-squared
chi.expected.values.cauchy <- dist.cauchy(chi.res.hist$mids, dist.cauchy.opt$par[1], dist.cauchy.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
chi.squared.statistic.cauchy <- sum((chi.res.hist$counts - chi.expected.values.cauchy)^2 / chi.expected.values.cauchy)
chi.squared.pvalue.cauchy <- 1-pchisq(chi.squared.statistic.cauchy, length(chi.res.hist$counts)-3)
# Kolmogorov-Smirnov
ks.expected.values.cauchy <- dist.cauchy(ks.res.hist$mids, dist.cauchy.opt$par[1], dist.cauchy.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
simul.cauchy <- c()
for (i in seq_along(ks.res.hist$mids)) {
  simul.cauchy <- c(simul.cauchy, rep(ks.res.hist$mids[i], round(ks.expected.values.cauchy[i], 0)))
}
ks.cauchy <- ks.test(data, simul.cauchy)
g.max.cauchy <- as.numeric(ks.cauchy$statistic)
KS.cauchy <- as.numeric(ks.cauchy$p.value)

# cumulative.expected.values.cauchy <- c(expected.values.cauchy[1])
# for (i in 1+seq_along(expected.values.cauchy)) {
#   cumulative.expected.values.cauchy[i] <- cumulative.expected.values.cauchy[i-1] + expected.values.cauchy[i]
# }
# cumulative.expected.values.cauchy <- cumulative.expected.values.cauchy/sum(expected.values.cauchy)
# cumulative.expected.values.cauchy <- cumulative.expected.values.cauchy[!is.na(cumulative.expected.values.cauchy)]
# g.max.cauchy <- max(abs(cumulative.data - cumulative.expected.values.cauchy))
# if (g.max.cauchy < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
#   KS.cauchy <- "Accept"
# } else {KS.cauchy <- "Reject"}
# parameter estimate
par.1.cauchy <- dist.cauchy.opt$par[1]
par.2.cauchy <- dist.cauchy.opt$par[2]
# parameter estimate standard error
par.1.se.cauchy <- sqrt(diag(solve(numDeriv::hessian(log.dist.cauchy, x=dist.cauchy.opt$par, r=data))))[1]
par.2.se.cauchy <- sqrt(diag(solve(numDeriv::hessian(log.dist.cauchy, x=dist.cauchy.opt$par, r=data))))[2]
# mean dispersal distance # the mean is undefined, so we used the median in the case of Cauchy
mean.cauchy <- par.1.cauchy
mean.stderr.cauchy <- par.1.se.cauchy
# variance
variance.cauchy <- "undefined"
variance.stderr.cauchy <- "undefined"
# skewness
skewness.cauchy <- "undefined"
skewness.stderr.cauchy <- "undefined"
# kurtosis
kurtosis.cauchy <- "undefined"
kurtosis.stderr.cauchy <- "undefined"
# output
res <- data.frame(aic.cauchy, aicc.cauchy, bic.cauchy,
                            chi.squared.statistic.cauchy, chi.squared.pvalue.cauchy,g.max.cauchy, KS.cauchy,
                            par.1.cauchy, par.1.se.cauchy, par.2.cauchy, par.2.se.cauchy,
                            mean.cauchy, mean.stderr.cauchy, variance.cauchy, variance.stderr.cauchy,
                            skewness.cauchy, skewness.stderr.cauchy, kurtosis.cauchy, kurtosis.stderr.cauchy)
cauchy.values <- list("opt" = dist.cauchy.opt, "res" = res)
}
