#' dispfit: A package to Fit Distributions to Dispersal Data
#'
#' The dispfit package fits several pre-defined distributions to dispersal or movement data,
#' computing several estimators: AIC, AICc, BIC, Chi-squared, and Kolgomorov-Smirnov.
#' It also estimates the parameter(s) value(s) and SE of each distribution,
#' as well as its Mean, Variance, Skewness, and Kurtosis.
#'
#'@details
#' The dispfit package fits 11 well-known distributions for estimating dispersal kernels
#' (Clark et al., 1999; Nathan et al., 2012) (Table 1).
#' The simplest functions considered are the single‐parameter Rayleigh and Exponential,
#' which are particularly popular in mathematical developments of theory concerning spatial dynamics
#' (O'Dwyer & Green 2010; Gilbert et al. 2014; Harsch et al. 2014).
#' The remaining 9 functions are two-parameter distributions which are often referred to better
#' represent real dispersal kernels than Rayleigh and Exponential functions
#' (Bullock and Clarke, 2000; Clark et al., 1999).
#' Apart from these 11 distributions, dispfit also includes three univariate functions
#' derived from the so-called generalized extreme value distribution (GEV),
#' and commonly used in extreme value statistics (REFs).
#'
#' @section Functions:
#' dispfit includes two main functions:
#' \tabular{ll}{
#'   \code{\link{dispersal.kernel}} \tab Fits several pre-defined distributions to dispersal or movement data. \cr
#'   \code{\link{dispersal.kernel.plot}} \tab Plots the distributions previously fitted by dispersal.kernel. \cr
#' }
#'
#' @section Distributions:
#' If a single point is the origin site of n dispersers drawn from a large population in continuous space,
#' then the dispersal distance of each dispersal is the Euclidian distance between the origin and end points.
#' The r dispersal distances of all dispersals reflect a continuous parametric distribution k,
#' or probability density function (pdf), that characterizes said population.
#' A dispersal kernel k is then defined as the pdf of the distribution of the final
#' location of a dispersal in relation to the source. Dispersal kernels implemented
#' in dispfit consider a two-dimensional space.
#' Under the assumption of rotational symmetry, direction can be suppressed so that k_(R,θ) (r,θ)
#' depends solely on one meaningful variable, r.
#' Such radially symmetric kernels have been referred to as ‘two-dimensional dispersal kernels’
#' (Cousens and Rawlinson, 2001) or ‘dispersal location kernels’, k_L (r) (Nathan et al., 2012),
#' and correspond to the dispersal kernel formulation adopted in dispfit.
#' The alternative ‘one-dimensional dispersal kernel’ (Cousens and Rawlinson, 2001)
#' or ‘dispersal distance kernel’,  k_D (r) {\cite{(Nathan et al., 2012)}},
#' provides the probability that a dispersion ends at a distance r away from the origin,
#' irrespective of direction, and is related to the ‘dispersal location kernel’
#' such that \ifelse{html}{\out{<i>k<sub>D</sub>(r) = 2&pi; r k<sub>L</sub>(r)</i>.}}
#'
#' \describe{
#' \item{Rayleigh}{\deqn{f(r) = (1/(\pi a^2)) exp(-(r/a)^2)}}
#' \item{Exponential}{\deqn{f(r) = (1/(2\pi a r)) exp(-r/a)}}
#' \item{Generalized Normal}{\deqn{f(r) = (b/(2\pi (a^2) \Gamma(2/b))) exp(-(r/a)^b)}}
#' \item{Bivariate Student’s t (2\emph{Dt)}}{\deqn{f(r) =  ((b-1) / (\pi (a^2))) ((1 + (r^2)/(a^2))^(-b))}}
#' \item{Geometric}{\deqn{f(r) = (((b - 2)(b - 1)) / (2\pi (a^2))) ((1 + (r / a)) ^ -b)}}
#' item{Logistic}{\deqn{f(r) = (b / (2\pi (a^2) \Gamma(2/b) \Gamma(1-(2/b)) )) ((1 + ((r^b) / (a^b)))^(-1))}}
#' \item{Lognormal}{\deqn{f(r) = (1 / (((2\pi) ^ (3/2)) (b (r ^ 2)))) exp(-(log(r / a)^2) / (2 * (b ^ 2)))}}
#' \item{Wald}{\deqn{f(r) = (\sqrt(b)/\sqrt(8 (\pi^3) (r^5))) exp(-(b ((r - a)^2))/(2 (a^2) * r))}}
#' \item{Weibull}{\deqn{f(r) = (b/(2\pi a^b)) (r^(b-2)) exp(-(r^b/a^b))}}
#' \item{Gamma}{\deqn{f(r) = (1 / (2\pi (a^2) \Gamma(b))) ((r/a)^(b-2)) * exp(-r/a)}}
#' \item{Log-sech}{\deqn{f(r) = (1 / ((\pi^2) b (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))}}
#' }
#'
#'
#' @author António Proença-Ferreira, \email{antoniomiguelpferreira@@gmail.com}
#' @docType package
#' @name dispfit
