#' dispfit: A package to Fit Distributions to Dispersal Data
#'
#' The dispfit package fits several pre-defined distributions to dispersal or movement data,
#' computing several estimators: AIC, AICc, BIC, Chi-squared, and Kolgomorov-Smirnov.
#' It also estimates the parameter(s) value(s) and CI of each distribution,
#' as well as its Mean, Variance, Skewness, and Kurtosis.
#'
#'@details
#' The dispfit package fits 9 well-known distributions for estimating dispersal kernels
#' (Clark et al., 1999; Nathan et al., 2012) (Table 1).
#' The simplest functions considered are the single‐parameter Rayleigh and Exponential,
#' which are particularly popular in mathematical developments of theory concerning spatial dynamics
#' (O'Dwyer & Green 2010; Gilbert et al. 2014; Harsch et al. 2014).
#' The remaining 7 functions are two-parameter distributions which are often referred to better
#' represent real dispersal kernels than Rayleigh and Exponential functions
#' (Bullock and Clarke, 2000; Clark et al., 1999).
#'
#' @section Package functions:
#' dispfit includes two main functions:
#' \tabular{ll}{
#'   \code{\link{dispersal.kernel}} \tab Fits several pre-defined distributions to dispersal or movement data. \cr
#'   \code{\link{plot.dispfit}} \tab Plots the distributions previously fitted by dispersal.kernel. \cr
#' }
#'
#' @section Probability density function:
#' Assuming that a single point is the origin site of all dispersers,
#' then the dispersal distance of each disperser is the Euclidian distance between the origin and its end point.
#' The dispersal distances of all dispersers reflect a continuous parametric distribution,
#' or probability density function (pdf), that characterizes the studied population.
#' A dispersal kernel is then defined as the pdf of the distribution of the values
#' of the Euclidean distances between the source and the final location of a dispersal event.
#' There are several characterizations of a dispersal kernel,
#' for instance Nathan et al. (2012) distinguish between “dispersal distance kernel, KD”,
#' and “dispersal location kernel, KL”.
#'
#' @section Distributions
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
#' @authors António Proença-Ferreira, \email{antoniomiguelpferreira@@gmail.com}
#' Luís Borda-de-Água
#' Miguel Porto
#' António Mira
#' Francisco Moreira
#' Ricardo Pita
#' @docType package
#'
"_PACKAGE"
#> [1] "_PACKAGE"
