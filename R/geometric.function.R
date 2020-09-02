geometric.function <- function (data, chi.res.hist, ks.res.hist) {
log.dist.geometric <- function (r, par) {
  a <- par[1]
  b <- par[2]
  fgeometric <- (((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
  -sum(log(fgeometric)) ##
}
dist.geometric <- function (r, a, b) {
  fgeometric <- 2*pi*r*(((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
}
# initial values estimation
while (TRUE) {
  SANN.geometric.opt <- optim (par = c(1, 2.000001), ## valor inicial para o "a"
                               fn = log.dist.geometric, ## função a minimizar
                               r = data,
                               method = "SANN",
                               # lower = c(0, 0),
                               control = list(maxit = 10000))
  try.geometric <- try(
    dist.geometric.opt.try <- optim (par = c(SANN.geometric.opt$par[1], SANN.geometric.opt$par[2]), ## valor inicial para o "a"
                                     fn = log.dist.geometric, ## função a minimizar
                                     r = data,
                                     method = "L-BFGS-B",
                                     lower = c(0.000001, 2.000001),
                                     upper = c(Inf, Inf),
                                     control = list(maxit = 10000)),
    silent=T)

  if (class(try.geometric) != "try-error") {
    dist.geometric.opt.try
    break
  }
}
# optimization procedure
dist.geometric.opt <- dist.geometric.opt.try
# dist.geometric.opt <- optim (par = c(0.00001, 3.00001), ## valor inicial para o "a"
#                       fn = log.dist.geometric, ## função a minimizar
#                       r = data, ## dados
#                       lower = c(0.00001, 3.00001),
#                       control = list(maxit = 10000), ## limite superior do parametro
#                       hessian = T)

# output values
# AIC
aic.geometric <- 2*length(dist.geometric.opt$par) + 2 * dist.geometric.opt$value
# AICc
aicc.geometric <- aic.geometric + (2 * length(dist.geometric.opt$par)^2 + 2 * length(dist.geometric.opt$par))/(length(data) - length(dist.geometric.opt$par) - 1 )
# BIC
bic.geometric <-  2 * dist.geometric.opt$value + length(dist.geometric.opt$par)*log(length(data))
# Chi-squared
chi.expected.values.geometric <- dist.geometric(chi.res.hist$mids, dist.geometric.opt$par[1], dist.geometric.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
chi.squared.statistic.geometric <- sum((chi.res.hist$counts - chi.expected.values.geometric)^2 / chi.expected.values.geometric)
chi.squared.pvalue.geometric <- 1-pchisq(chi.squared.statistic.geometric, length(chi.res.hist$counts)-3)
# Kolmogorov-Smirnov
ks.expected.values.geometric <- dist.geometric(ks.res.hist$mids, dist.geometric.opt$par[1], dist.geometric.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
simul.geometric <- c()
for (i in seq_along(ks.res.hist$mids)) {
  simul.geometric <- c(simul.geometric, rep(ks.res.hist$mids[i], round(ks.expected.values.geometric[i], 0)))
}
ks.geometric <- ks.test(data, simul.geometric)
g.max.geometric <- as.numeric(ks.geometric$statistic)
KS.geometric <- as.numeric(ks.geometric$p.value)

# cumulative.expected.values.geometric <- c(expected.values.geometric[1])
# for (i in 1+seq_along(expected.values.geometric)) {
#   cumulative.expected.values.geometric[i] <- cumulative.expected.values.geometric[i-1] + expected.values.geometric[i]
# }
# cumulative.expected.values.geometric <- cumulative.expected.values.geometric/sum(expected.values.geometric)
# cumulative.expected.values.geometric <- cumulative.expected.values.geometric[!is.na(cumulative.expected.values.geometric)]
# g.max.geometric <- max(abs(cumulative.data - cumulative.expected.values.geometric))
# if (g.max.geometric < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
#   KS.geometric <- "Accept"
# } else {KS.geometric <- "Reject"}
# parameter estimate
par.1.geometric <- dist.geometric.opt$par[1]
par.2.geometric <- dist.geometric.opt$par[2]
# parameter estimate standard error
par.1.se.geometric <- sqrt(diag(solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data))))[1]
par.2.se.geometric <- sqrt(diag(solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data))))[2]
# mean dispersal distance
if (dist.geometric.opt$par[2] >= 3) {
  mean.geometric <- (2 * dist.geometric.opt$par[1]) / (dist.geometric.opt$par[2]-3)
} else {
  mean.geometric <- "Infinite Value"
}
if (dist.geometric.opt$par[2] >= 3) {
  mean.stderr.geometric <- msm::deltamethod(~ (2 * x1) / (x2-3),
                                            mean = dist.geometric.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
} else {
  mean.stderr.geometric <- "Infinite Value"
}
# variance
if (dist.geometric.opt$par[2] >= 4) {
  variance.geometric <- sqrt(6) * dist.geometric.opt$par[1] * sqrt(1/((dist.geometric.opt$par[2]-4) * (dist.geometric.opt$par[2]-3)))
} else {
  variance.geometric <- "Infinite Value"
}
if (dist.geometric.opt$par[2] >= 4) {
  variance.stderr.geometric <- msm::deltamethod(~ sqrt(6) * x1 * sqrt(1/((x2-4) * (x2-3))),
                                                mean = dist.geometric.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
} else {
  variance.stderr.geometric <- "Infinite Value"
}
# standard deviation
if (dist.geometric.opt$par[2] >= 4) {
  stdev.geometric <- sqrt(sqrt(6) * dist.geometric.opt$par[1] * sqrt(1/((dist.geometric.opt$par[2]-4) * (dist.geometric.opt$par[2]-3))))
} else {
  stdev.geometric <- "Infinite Value"
}
if (dist.geometric.opt$par[2] >= 4) {
  stdev.stderr.geometric <- msm::deltamethod(~ sqrt(sqrt(6) * x1 * sqrt(1/((x2-4) * (x2-3)))),
                                             mean = dist.geometric.opt$par,
                                             cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
} else {
  stdev.stderr.geometric <- "Infinite Value"
}
# skewness
if (dist.geometric.opt$par[2] >= 5) {
  skewness.geometric <- 2 * 3^(1/3) * dist.geometric.opt$par[1] * (1/((dist.geometric.opt$par[2] - 5) * (dist.geometric.opt$par[2] - 4) * (dist.geometric.opt$par[2] - 3)))^(1/3)
} else {
  skewness.geometric <- "Infinite Value"
}
if (dist.geometric.opt$par[2] >= 5) {
  skewness.stderr.geometric <- msm::deltamethod(~ 2 * 3^(1/3) * x1 * (1/((x2 - 5) * (x2 - 4) * (x2 - 3)))^(1/3),
                                                mean = dist.geometric.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
} else {
  skewness.stderr.geometric <- "Infinite Value"
}
# kurtosis
if (dist.geometric.opt$par[2] >= 6) {
  kurtosis.geometric <- 2^(3/4) * 15^(1/4) * dist.geometric.opt$par[1] * (1/((dist.geometric.opt$par[2] - 6) * (dist.geometric.opt$par[2] - 5) * (dist.geometric.opt$par[2] - 4) * (dist.geometric.opt$par[2] - 3)))^(1/4)
} else {
  kurtosis.geometric <- "Infinite Value"
}
if (dist.geometric.opt$par[2] >= 6) {
  kurtosis.stderr.geometric <- msm::deltamethod(~ 2^(3/4) * 15^(1/4) * x2 * (1/((x2 - 6) * (x2 - 5) * (x2 - 4) * (x2 - 3)))^(1/4),
                                                mean = dist.geometric.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
} else {
  kurtosis.stderr.geometric <- "Infinite Value"
}
# output
res <- data.frame(aic.geometric, aicc.geometric, bic.geometric,
                               chi.squared.statistic.geometric, chi.squared.pvalue.geometric,g.max.geometric, KS.geometric,
                               par.1.geometric, par.1.se.geometric, par.2.geometric, par.2.se.geometric,
                               mean.geometric, mean.stderr.geometric, variance.geometric, variance.stderr.geometric,
                               skewness.geometric, skewness.stderr.geometric, kurtosis.geometric, kurtosis.stderr.geometric)
geometric.values <- list("opt" = dist.geometric.opt, "res" = res)
}
