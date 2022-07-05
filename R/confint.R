compute.par.limits <- function(lower.limits, upper.limits, pars, data, parameter) {
	if(inherits(lower.limits, "list")) {
		if(inherits(lower.limits[[parameter]], "function"))
			par2.lower <- lower.limits[[parameter]](pars, data)
		else
			par2.lower <- lower.limits[[parameter]]
	} else if(inherits(lower.limits, "numeric")) {
		par2.lower <- lower.limits[parameter]
	}

	if(inherits(upper.limits, "list")) {
		if(inherits(upper.limits[[parameter]], "function"))
			par2.upper <- upper.limits[[parameter]](pars, data)
		else
			par2.upper <- upper.limits[[parameter]]
	} else if(inherits(upper.limits, "numeric")) {
		par2.upper <- upper.limits[parameter]
	}
	return(c(par2.lower, par2.upper))
}

confint.dispfit <- function(init.pars, logdistfun, data, lower.limits=c(0, 0), upper.limits=c(100000, 100000), confidence.level, debug=TRUE) {
	twopars <- length(init.pars$par) > 1
	
	par1 <- init.pars$par[1]
	# parameter estimate standard error
	par.1.se <- sqrt(diag(solve(numDeriv::hessian(logdistfun, x=init.pars$par, r=data))))[1]

	if(twopars) {
		par2 <- init.pars$par[2]
		par.2.se <- sqrt(diag(solve(numDeriv::hessian(logdistfun, x=init.pars$par, r=data))))[2]
	}
	
	if(par1 > 100000 || (twopars && par2 > 100000)) {
		warning("Distribution parameters likely diverged")
		return(c(
			par1.CIlow=NA,
			par1.CIupp=NA,
			par2.CIlow=NA,
			par2.CIupp=NA
		))
	}
	
	# parameter estimate confidence intervals
	log.dist.ci <- function (a, b=NA, r) {
		return(logdistfun(par=c(a, b), r))
	}
	n.se <- 30
	len <- 1000

	par.1.ini <- par1 - n.se * par.1.se
	par.1.fin <- par1 + n.se * par.1.se
	step <- (par.1.fin - par.1.ini) / len
	threshold <- init.pars$value + qchisq(confidence.level, 1) / 2
	max.trials <- 10000

	# Par 1 CI LOWER
	par1.test <- par1
	if(debug) par.1.prof <- matrix(ncol=2, nrow=0)
	count <- 0
	last.value <- NA
	last.par <- NA

	repeat {
		if(twopars) {
			limits <- compute.par.limits(lower.limits, upper.limits, c(par1.test, NA), data, parameter=2)
			if(limits[2] > limits[1] && is.finite(log.dist.ci(par1.test, par2, data))) {
	#				tryCatch({
				prev.value <- last.value
				last.value <- optim(log.dist.ci, par = par2, a = par1.test, r = data, lower=limits[1], upper=limits[2], method = "Brent")$value
				if(debug) par.1.prof <- rbind(par.1.prof, c(par1.test, last.value))
	#				}, warning=function(w) browser())
			}
		} else {
			prev.value <- last.value
			last.value <- log.dist.ci(a = par1.test, r = data)
			if(debug) par.1.prof <- rbind(par.1.prof, c(par1.test, last.value))
		}
		prev.par <- last.par
		last.par <- par1.test
		par1.test <- par1.test - step

		if(par1.test <= lower.limits[1]) {	# ups, we crossed the lower limit
			par1.test <- par1.test + step	# undo step
			repeat {
				step <- step / 10	# increase resolution
				if(debug) message("Reduced step")
				if(par1.test - step > lower.limits[1]) break
			}
			par1.test <- par1.test - step
		}

		count <- count + 1
		if(last.value > threshold || count > max.trials || step < 1.e-12) break
	}
	if(debug) {
		par(mfrow=c(2, 2))
		plot(par.1.prof, pch=19, cex=0.7)
		abline(v=approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y, h=threshold, lwd=2)
	}

	if(count > max.trials || step < 1.e-12) {
		par.1.CIlow <- lower.limits[1]
		warning("Lower CI for 'a' is not accurate, I've given up after ", max.trials, " trials.")
	} else
		par.1.CIlow <- approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y

	# Par 1 CI UPPER
	step <- (par.1.fin - par.1.ini) / len
	par1.test <- par1
	if(debug) par.1.prof <- matrix(ncol=2, nrow=0)
	count <- 0
	last.value <- NA
	last.par <- NA

	repeat {
		if(twopars) {
			limits <- compute.par.limits(lower.limits, upper.limits, c(par1.test, NA), data, parameter=2)

			if(limits[2] > limits[1] && is.finite(log.dist.ci(par1.test, par2, data))) {
	#				tryCatch({
				prev.value <- last.value
				last.value <- optim(log.dist.ci, par = par2, a = par1.test, r = data, lower=limits[1], upper=limits[2], method = "Brent")$value
				if(debug) par.1.prof <- rbind(par.1.prof, c(par1.test, last.value))
	#				}, warning=function(w) browser())
			}
		} else {
			prev.value <- last.value
			last.value <- log.dist.ci(a = par1.test, r = data)
			if(debug) par.1.prof <- rbind(par.1.prof, c(par1.test, last.value))
		}
		prev.par <- last.par
		last.par <- par1.test
		par1.test <- par1.test + step
		count <- count + 1
		if(last.value > threshold || count > max.trials) break;
	}

	if(debug) {
		plot(par.1.prof, pch=19, cex=0.7)
		abline(v=approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y, h=threshold, lwd=2)
	}
	if(count > max.trials) {
		par.1.CIupp <- Inf
		warning("Upper CI for 'a' is not accurate, I've given up after ", max.trials, " trials.")
	} else
		par.1.CIupp <- approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y
		

###################

	if(twopars) {
		par.2.ini <- par2 - n.se * par.2.se
		par.2.fin <- par2 + n.se * par.2.se
		step <- (par.2.fin - par.2.ini) / len
		threshold <- init.pars$value + qchisq(confidence.level, 1) / 2
		max.trials <- 10000

		# Par 2 CI LOWER
		par2.test <- par2
		if(debug) par.2.prof <- matrix(ncol=2, nrow=0)
		count <- 0
		last.value <- NA
		last.par <- NA

		repeat {
			limits <- compute.par.limits(lower.limits, upper.limits, c(NA, par2.test), data, parameter=1)

			if(limits[2] > limits[1] && is.finite(log.dist.ci(par1, par2.test, data))) {
	#				tryCatch({
				prev.value <- last.value
				last.value <- optim(log.dist.ci, par = par1, b = par2.test, r = data, lower=limits[1], upper=limits[2], method = "Brent")$value
				if(debug) par.2.prof <- rbind(par.2.prof, c(par2.test, last.value))
	#				}, warning=function(w) browser())
			}
			prev.par <- last.par
			last.par <- par2.test
			par2.test <- par2.test - step

			if(par2.test <= lower.limits[2]) {	# ups, we crossed the lower limit
				par2.test <- par2.test + step	# undo step
				repeat {
					step <- step / 10	# increase resolution
					if(debug) message("Reduced step")
					if(par2.test - step > lower.limits[2]) break
				}
				par2.test <- par2.test - step
			}

			count <- count + 1
			if(last.value > threshold || count > max.trials || step < 1.e-12) break
		}
		if(debug) {
			plot(par.2.prof, pch=19, cex=0.7)
			abline(v=approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y, h=threshold, lwd=2)
		}

		if(count > max.trials || step < 1.e-12) {
			par.2.CIlow <- lower.limits[2]
			warning("Lower CI for 'b' is not accurate, I've given up after ", max.trials, " trials.")
		} else
			par.2.CIlow <- approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y

		# Par 2 CI UPPER
		step <- (par.2.fin - par.2.ini) / len
		par2.test <- par2
		if(debug) par.2.prof <- matrix(ncol=2, nrow=0)
		count <- 0
		last.value <- NA
		last.par <- NA

		repeat {
			limits <- compute.par.limits(lower.limits, upper.limits, c(NA, par2.test), data, parameter=1)

			if(limits[2] > limits[1] && is.finite(log.dist.ci(par1, par2.test, data))) {
	#				tryCatch({
				prev.value <- last.value
				last.value <- optim(log.dist.ci, par = par1, b = par2.test, r = data, lower=limits[1], upper=limits[2], method = "Brent")$value
				if(debug) par.2.prof <- rbind(par.2.prof, c(par2.test, last.value))
	#				}, warning=function(w) browser())
			}
			prev.par <- last.par
			last.par <- par2.test
			par2.test <- par2.test + step
			count <- count + 1
			if(last.value > threshold || count > max.trials) break;
		}

		if(debug) {
			plot(par.2.prof, pch=19, cex=0.7)
			abline(v=approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y, h=threshold, lwd=2)
#			readline("ENTER to go on.")
		}

		if(count > max.trials) {
			par.2.CIupp <- Inf
			warning("Upper CI for 'b' is not accurate, I've given up after ", max.trials, " trials.")
		} else
			par.2.CIupp <- approx(c(prev.value, last.value), c(prev.par, last.par), xout = threshold)$y
	} else {
	  par.2.CIlow <- NA
	  par.2.CIupp <- NA
	}

	return(c(
		par1.CIlow=par.1.CIlow,
		par1.CIupp=par.1.CIupp,
		par2.CIlow=par.2.CIlow,
		par2.CIupp=par.2.CIupp
	))
}
