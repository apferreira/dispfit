gamma.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.gamma <- function (par, r) {
    a <- par[1] ## scale
    b <- par[2] ## shape
    if(a < 0 || b < 0) return(Inf)
    
    fgamma <- (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
    fgamma[fgamma == 0] <- .Machine$double.xmin
    -sum(log(fgamma)) ## Negative Log Likelihood
  }
  dist.gamma <- function (r, a, b) {
    fgamma <- 2*pi*r * (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
  }
  # initial values estimation
  shape <- (mean(data)/sd(data))
  scale <- sd(data)^2/mean(data)
  # optimization procedure
  dist.opt <- optim (par = c(scale, shape), ##
                           fn = log.dist.gamma, ## função a minimizar
                           r = data,
                           method = "L-BFGS-B",
                           lower = c(0.00001, 0.00001)
  )
  # output values
  # AIC
  aic.gamma <- 2 * length(dist.opt$par) + 2 * dist.opt$value
  # AICc
  aicc.gamma <- aic.gamma + (2 * length(dist.opt$par)^2 + 2 * length(dist.opt$par))/(length(data) - length(dist.opt$par) - 1 )
  # BIC
  bic.gamma <-  2 * dist.opt$value + length(dist.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.gamma <- dist.gamma(chi.res.hist$mids, dist.opt$par[1],  dist.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.gamma <- sum((chi.res.hist$counts - chi.expected.values.gamma)^2 / chi.expected.values.gamma)
  chi.squared.pvalue.gamma <- 1-pchisq(chi.squared.statistic.gamma, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.gamma <- dist.gamma(ks.res.hist$mids, dist.opt$par[1],  dist.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
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

  par2.upper.limit <- function(pars, data) {
    # 169 is near the maximum value that can be used in gamma(x) without returning Inf
    min(169, (2 * log(max(data) / pars[1]) + log(.Machine$double.xmax)) / log(max(data) / pars[1]))
  }
  CI <- confint.dispfit(dist.opt, log.dist.gamma, data=data, lower=c(0, 0), upper=list(100000, par2.upper.limit), confidence.level=confidence.level)
  
  # mean dispersal distance
  mean.gamma <- dist.opt$par[1] * dist.opt$par[2]
  mean.stderr.gamma <- msm::deltamethod(~ x1 * x2,
                                        mean = dist.opt$par,
                                        cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.opt$par, r=data)) )
  # variance
  variance.gamma <- dist.opt$par[1]^2 * dist.opt$par[2]
  variance.stderr.gamma <- msm::deltamethod(~ x2*x1^2,
                                            mean = dist.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.opt$par, r=data)) )
  # standard deviation
  stdev.gamma <- sqrt(dist.opt$par[1]^2 * dist.opt$par[2])
  stdev.stderr.gamma <- msm::deltamethod(~ sqrt(x2*x1^2),
                                         mean = dist.opt$par,
                                         cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.opt$par, r=data)) )
  # skewness
  skewness.gamma <- 2/sqrt(dist.opt$par[2])
  skewness.stderr.gamma <- msm::deltamethod(~ 2/sqrt(x2),
                                            mean = dist.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.opt$par, r=data)) )
  # kurtosis
  kurtosis.gamma <- 6/dist.opt$par[2]
  kurtosis.stderr.gamma <- msm::deltamethod(~ 6/x2,
                                            mean = dist.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.opt$par, r=data)) )
  # output
  res <- data.frame(aic.gamma, aicc.gamma, bic.gamma,
                    chi.squared.statistic.gamma, chi.squared.pvalue.gamma,g.max.gamma, KS.gamma,
                    dist.opt$par[1], CI["par1.CIlow"], CI["par1.CIupp"],
                    dist.opt$par[2], CI["par2.CIlow"], CI["par2.CIupp"],
                    mean.gamma, mean.stderr.gamma, stdev.gamma, stdev.stderr.gamma,
                    skewness.gamma, skewness.stderr.gamma, kurtosis.gamma, kurtosis.stderr.gamma)
  gamma.values <- list("opt" = dist.opt, "res" = res)
}
