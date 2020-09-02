weibull.function <- function (data, chi.res.hist, ks.res.hist) {
log.dist.weibull <- function (r, par) {
  a <- par[1] ## scale
  b <- par[2] ## shape
  fw <- (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
  -sum(log(fw)) ##
}
dist.weibull <- function (r, a, b) {
  fw <- 2*pi*r * (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
}
# initial values estimation
m <- mean(log(data))
v <- var(log(data))
shape <- 1.2/sqrt(v)
scale <- exp(m + 0.572/shape)
# optimization procedure
dist.weibull.opt <- optim (par = c(scale, shape), ## valor inicial para o "a"
                           fn = log.dist.weibull, ## função a minimizar
                           r = data, ## dados
                           method = "Nelder-Mead",
                           # lower = c(0.00001, 0.00001),
                           control = list(maxit = 10000))
# output values
# AIC
aic.weibull <- 2*length(dist.weibull.opt$par) + 2 * dist.weibull.opt$value
# AICc
aicc.weibull <- aic.weibull + (2 * length(dist.weibull.opt$par)^2 + 2 * length(dist.weibull.opt$par))/(length(data) - length(dist.weibull.opt$par) - 1 )
# BIC
bic.weibull <-  2 * dist.weibull.opt$value + length(dist.weibull.opt$par)*log(length(data))
# Chi-squared
chi.expected.values.weibull <- dist.weibull(chi.res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
chi.squared.statistic.weibull <- sum((chi.res.hist$counts - chi.expected.values.weibull)^2 / chi.expected.values.weibull)
chi.squared.pvalue.weibull <- 1-pchisq(chi.squared.statistic.weibull, length(chi.res.hist$counts)-3)
# Kolmogorov-Smirnov
ks.expected.values.weibull <- dist.weibull(ks.res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
simul.weibull <- c()
for (i in seq_along(ks.res.hist$mids)) {
  simul.weibull <- c(simul.weibull, rep(ks.res.hist$mids[i], round(ks.expected.values.weibull[i], 0)))
}
ks.weibull <- ks.test(data, simul.weibull)
g.max.weibull <- as.numeric(ks.weibull$statistic)
KS.weibull <- as.numeric(ks.weibull$p.value)

# cumulative.expected.values.weibull <- c(expected.values.weibull[1])
# for (i in 1+seq_along(expected.values.weibull)) {
#   cumulative.expected.values.weibull[i] <- cumulative.expected.values.weibull[i-1] + expected.values.weibull[i]
# }
# cumulative.expected.values.weibull <- cumulative.expected.values.weibull/sum(expected.values.weibull)
# cumulative.expected.values.weibull <- cumulative.expected.values.weibull[!is.na(cumulative.expected.values.weibull)]
# g.max.weibull <- max(abs(cumulative.data - cumulative.expected.values.weibull))
# if (g.max.weibull < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
#   KS.weibull <- "Accept"
# } else {KS.weibull <- "Reject"}
# parameter estimate
par.1.weibull <- dist.weibull.opt$par[1]
par.2.weibull <- dist.weibull.opt$par[2]
# parameter estimate standard error
par.1.se.weibull <- sqrt(diag(solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data))))[1]
par.2.se.weibull <- sqrt(diag(solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data))))[2]
# mean dispersal distance ## from Austerlitz 2004
mean.weibull <- dist.weibull.opt$par[1] * (gamma(1 + 1/dist.weibull.opt$par[2]))
mean.stderr.weibull <- msm::deltamethod(~ x1 * (gamma(1 + 1/x2)), mean = dist.weibull.opt$par, cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
# variance
variance.weibull <- dist.weibull.opt$par[1]^2*(gamma(1+2/dist.weibull.opt$par[2])-gamma(1+1/dist.weibull.opt$par[2])^2)
variance.stderr.weibull <- msm::deltamethod(~ x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2), mean = dist.weibull.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
# skewness
skewness.weibull <- (dist.weibull.opt$par[1]^3*gamma(1+3/dist.weibull.opt$par[2])-3*mean.weibull*variance.weibull-mean.weibull^3)/variance.weibull^(3/2)
skewness.stderr.weibull <- msm::deltamethod(~ (x1^3*gamma(1+3/x2)-3*(x1 * (gamma(1 + 1/x2)))*(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))-(x1 * (gamma(1 + 1/x2)))^3)/(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2),
                                            mean = dist.weibull.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
# kurtosis
kurtosis.weibull <- (dist.weibull.opt$par[1]^4*gamma(1+4/dist.weibull.opt$par[2])-4*skewness.weibull*variance.weibull^(3/2)*mean.weibull-6*variance.weibull*mean.weibull^2-mean.weibull^4)/variance.weibull^2-3
kurtosis.stderr.weibull <- msm::deltamethod(~ (x1^4*gamma(1+4/x2) - 4 * ((x1^3*gamma(1+3/x2)-3*(x1 * (gamma(1 + 1/x2)))*(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))-(x1 * (gamma(1 + 1/x2)))^3)/(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2)) * (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2) * (x1 * (gamma(1 + 1/x2))) - 6 * (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2)) * (x1 * (gamma(1 + 1/x2)))^2 - (x1 * (gamma(1 + 1/x2)))^4) / (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^2-3,
                                            mean = dist.weibull.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
# output
res <- data.frame(aic.weibull, aicc.weibull, bic.weibull,
                             chi.squared.statistic.weibull, chi.squared.pvalue.weibull,g.max.weibull, KS.weibull,
                             par.1.weibull, par.1.se.weibull, par.2.weibull, par.2.se.weibull,
                             mean.weibull, mean.stderr.weibull, variance.weibull, variance.stderr.weibull,
                             skewness.weibull, skewness.stderr.weibull, kurtosis.weibull, kurtosis.stderr.weibull)
weibull.values <- list("opt" = dist.weibull.opt, "res" = res)
}
