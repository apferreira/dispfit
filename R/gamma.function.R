gamma.function <- function (data, chi.res.hist, ks.res.hist) {
log.dist.gamma <- function (par, r) {
  a <- par[1] ## scale
  b <- par[2] ## shape
  fgamma <- (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
  -sum(log(fgamma)) ## Negative Log Likelihood
}
dist.gamma <- function (r, a, b) {
  fgamma <- 2*pi*r * (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
}
# initial values estimation
shape <- (mean(data)/sd(data))
scale <- sd(data)^2/mean(data)
# optimization procedure
dist.gamma.opt <- optim (par = c(scale, shape), ##
                         fn = log.dist.gamma, ## função a minimizar
                         r = data,
                         method = "L-BFGS-B",
                         lower = c(0.00001, 0.00001)
)
# output values
# AIC
aic.gamma <- 2 * length(dist.gamma.opt$par) + 2 * dist.gamma.opt$value
# AICc
aicc.gamma <- aic.gamma + (2 * length(dist.gamma.opt$par)^2 + 2 * length(dist.gamma.opt$par))/(length(data) - length(dist.gamma.opt$par) - 1 )
# BIC
bic.gamma <-  2 * dist.gamma.opt$value + length(dist.gamma.opt$par)*log(length(data))
# Chi-squared
chi.expected.values.gamma <- dist.gamma(chi.res.hist$mids, dist.gamma.opt$par[1],  dist.gamma.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
chi.squared.statistic.gamma <- sum((chi.res.hist$counts - chi.expected.values.gamma)^2 / chi.expected.values.gamma)
chi.squared.pvalue.gamma <- 1-pchisq(chi.squared.statistic.gamma, length(chi.res.hist$counts)-3)
# Kolmogorov-Smirnov
ks.expected.values.gamma <- dist.gamma(ks.res.hist$mids, dist.gamma.opt$par[1],  dist.gamma.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
simul.gamma <- c()
for (i in seq_along(ks.res.hist$mids)) {
  simul.gamma <- c(simul.gamma, rep(ks.res.hist$mids[i], round(ks.expected.values.gamma[i], 0)))
}
ks.gamma <- ks.test(data, simul.gamma)
g.max.gamma <- as.numeric(ks.gamma$statistic)
KS.gamma <- as.numeric(ks.gamma$p.value)

# cumulative.expected.values.gamma <- c(expected.values.gamma[1])
# for (i in 1+seq_along(expected.values.gamma)) {
#   cumulative.expected.values.gamma[i] <- cumulative.expected.values.gamma[i-1] + expected.values.gamma[i]
# }
# cumulative.expected.values.gamma <- cumulative.expected.values.gamma/sum(expected.values.gamma)
# cumulative.expected.values.gamma <- cumulative.expected.values.gamma[!is.na(cumulative.expected.values.gamma)]
# g.max.gamma <- max(abs(cumulative.data - cumulative.expected.values.gamma))
# if (g.max.gamma < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
#   KS.gamma <- "Accept"
# } else {KS.gamma <- "Reject"}
# parameter estimate
par.1.gamma <- dist.gamma.opt$par[1]
par.2.gamma <- dist.gamma.opt$par[2]
# parameter estimate standard error
par.1.se.gamma <- sqrt(diag(solve(numDeriv    ::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[1]
par.2.se.gamma <- sqrt(diag(solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[2]
# mean dispersal distance
mean.gamma <- dist.gamma.opt$par[1] * dist.gamma.opt$par[2]
mean.stderr.gamma <- msm::deltamethod(~ x1 * x2,
                                      mean = dist.gamma.opt$par,
                                      cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
# variance
variance.gamma <- dist.gamma.opt$par[1]^2 * dist.gamma.opt$par[2]
variance.stderr.gamma <- msm::deltamethod(~ x2*x1^2,
                                          mean = dist.gamma.opt$par,
                                          cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
# skewness
skewness.gamma <- 2/sqrt(dist.gamma.opt$par[2])
skewness.stderr.gamma <- msm::deltamethod(~ 2/sqrt(x2),
                                          mean = dist.gamma.opt$par,
                                          cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
# kurtosis
kurtosis.gamma <- 6/dist.gamma.opt$par[2]
kurtosis.stderr.gamma <- msm::deltamethod(~ 6/x2,
                                          mean = dist.gamma.opt$par,
                                          cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
# output
res <- data.frame(aic.gamma, aicc.gamma, bic.gamma,
                           chi.squared.statistic.gamma, chi.squared.pvalue.gamma,g.max.gamma, KS.gamma,
                           par.1.gamma, par.1.se.gamma, par.2.gamma, par.2.se.gamma,
                           mean.gamma, mean.stderr.gamma, variance.gamma, variance.stderr.gamma,
                           skewness.gamma, skewness.stderr.gamma, kurtosis.gamma, kurtosis.stderr.gamma)
gamma.values <- list("opt" = dist.gamma.opt, "res" = res)
}
