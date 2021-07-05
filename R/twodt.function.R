twodt.function <- function (data, chi.res.hist, ks.res.hist) {
  log.dist.2dt <- function (r, par) {
    a <- par[1] ## scale parameter
    b <- par[2] ## shape parameter
    f2dt <- ((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
    -sum(log(f2dt)) ##
  }
  dist.2dt <- function (r, a, b) {
    f2dt <- 2*pi*r*((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
  }
  # initial values estimation
  while (TRUE) {
    SANN.2dt.opt <- optim (par = c(1, 1.000001), ## valor inicial para o "a"
                           fn = log.dist.2dt, ## função a minimizar
                           r = data,
                           method = "SANN",
                           # lower = c(0, 0),
                           control = list(maxit = 10000))
    try.2dt <- try(
      dist.2dt.opt.try <- optim (par = c(SANN.2dt.opt$par[1], SANN.2dt.opt$par[2]), ## valor inicial para o "a"
                                 fn = log.dist.2dt, ## função a minimizar
                                 r = data,
                                 method = "L-BFGS-B",
                                 lower = c(0.000001, 1.000001),
                                 upper = c(Inf, Inf),
                                 control = list(maxit = 10000)),
      silent=T)

    if (class(try.2dt) != "try-error") {
      dist.2dt.opt.try
      break
    }
  }
  # optimization procedure
  dist.2dt.opt <- dist.2dt.opt.try
  # dist.2dt.opt <- optim (par = c(0.001, 1.001), ## valor inicial para o "a"
  #                        fn = log.dist.2dt, ## função a minimizar
  #                        r = data, ## dados
  #                        control = list(maxit = 10000), ## limite superior do parametro
  #                        hessian = T)
  # output values
  # AIC
  aic.2dt <- 2*length(dist.2dt.opt$par) + 2 * dist.2dt.opt$value
  # AICc
  aicc.2dt <- aic.2dt + (2 * length(dist.2dt.opt$par)^2 + 2 * length(dist.2dt.opt$par))/(length(data) - length(dist.2dt.opt$par) - 1 )
  # BIC
  bic.2dt <-  2 * dist.2dt.opt$value + length(dist.2dt.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.2dt <- dist.2dt(chi.res.hist$mids, dist.2dt.opt$par[1], dist.2dt.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.2dt <- sum((chi.res.hist$counts - chi.expected.values.2dt)^2 / chi.expected.values.2dt)
  chi.squared.pvalue.2dt <- 1-pchisq(chi.squared.statistic.2dt, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.2dt <- dist.2dt(ks.res.hist$mids, dist.2dt.opt$par[1], dist.2dt.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.2dt <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.2dt <- c(simul.2dt, rep(ks.res.hist$mids[i], round(ks.expected.values.2dt[i], 0)))
  }
  ks.2dt <- ks.test(data, simul.2dt)
  g.max.2dt <- as.numeric(ks.2dt$statistic)
  KS.2dt <- as.numeric(ks.2dt$p.value)

  # cumulative.expected.values.2dt <- c(expected.values.2dt[1])
  # for (i in 1+seq_along(expected.values.2dt)) {
  #   cumulative.expected.values.2dt[i] <- cumulative.expected.values.2dt[i-1] + expected.values.2dt[i]
  # }
  # cumulative.expected.values.2dt <- cumulative.expected.values.2dt/sum(expected.values.2dt)
  # cumulative.expected.values.2dt <- cumulative.expected.values.2dt[!is.na(cumulative.expected.values.2dt)]
  # g.max.2dt <- max(abs(cumulative.data - cumulative.expected.values.2dt))
  # if (g.max.2dt < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
  #   KS.2dt <- "Accept"
  # } else {KS.2dt <- "Reject"}
  # parameter estimate
  par.1.2dt <- dist.2dt.opt$par[1]
  par.2.2dt <- dist.2dt.opt$par[2]
  # parameter estimate standard error
  par.1.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[1]
  par.2.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[2]
  # mean dispersal distance
  if (dist.2dt.opt$par[2] >= 3/2) {
    mean.2dt <- dist.2dt.opt$par[1] * (sqrt(pi)/2) * (exp(lgamma(dist.2dt.opt$par[2]-(3/2))-lgamma(dist.2dt.opt$par[2]-1)))
  } else {
    mean.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] >= 1.5) {
    mean.stderr.2dt <- msm::deltamethod(~ x1 * (sqrt(pi)/2)*(exp(lgamma(x2-(3/2))-lgamma(x2-1))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    mean.stderr.2dt <- "Infinite Value"
  }
  # variance
  if (dist.2dt.opt$par[2] > 2) {
    variance.2dt <- dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2))
  } else {
    variance.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] > 2) {
    variance.stderr.2dt <- msm::deltamethod(~ x1 * sqrt(1/(x2-2)), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    variance.stderr.2dt <- "Infinite Value"
  }
  # standard deviation
  if (dist.2dt.opt$par[2] > 2) {
    stdev.2dt <- sqrt(dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2)))
  } else {
    stdev.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] > 2) {
    stdev.stderr.2dt <- msm::deltamethod(~ sqrt(x1 * sqrt(1/(x2-2))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    stdev.stderr.2dt <- "Infinite Value"
  }
  # skewness
  if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
    skewness.2dt <- (dist.2dt.opt$par[1]*(gamma(2.5) * exp(lgamma(dist.2dt.opt$par[2]-2.5) - lgamma(dist.2dt.opt$par[2]-1)))^(1/3))
  } else {
    skewness.2dt <- "Infinite Value"
  }
  if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
    skewness.stderr.2dt <- msm::deltamethod(~ x1 * (gamma(2.5) * exp(lgamma(x2-2.5) - lgamma(x2-1)))^(1/3), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    skewness.stderr.2dt <- "Infinite Value"
  }
  # kurtosis
  if (2*(dist.2dt.opt$par[2]-1) > 3) {
    kurtosis.2dt <- 2^(1/4) * dist.2dt.opt$par[1] * (1/((dist.2dt.opt$par[2] - 3) * (dist.2dt.opt$par[2] - 2)))^(1/4)
  } else {
    kurtosis.2dt <- "Infinite Value"
  }
  if (2*(dist.2dt.opt$par[2]-1) > 3) {
    kurtosis.stderr.2dt <- msm::deltamethod(~ 2^(1/4) * x1 * (1/((x2 - 3) * (x2 - 2)))^(1/4),
                                            mean = dist.2dt.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    kurtosis.stderr.2dt <- "Infinite Value"
  }
  # output
  res <- data.frame(aic.2dt, aicc.2dt, bic.2dt,
                             chi.squared.statistic.2dt, chi.squared.pvalue.2dt,g.max.2dt, KS.2dt,
                             par.1.2dt, par.1.se.2dt, par.2.2dt, par.2.se.2dt,
                             mean.2dt, mean.stderr.2dt, stdev.2dt, stdev.stderr.2dt,
                    skewness.2dt, skewness.stderr.2dt, kurtosis.2dt, kurtosis.stderr.2dt)
twodt.values <- list("opt" = dist.2dt.opt, "res" = res)
}
