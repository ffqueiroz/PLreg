#' @keywords internal
pslash <- function(q, nu){
  aux <-  function(x){
    s_aux <- x^2/2
    beta_aux <- nu + (1/2)
    gama_aux <- zipfR::Igamma(beta_aux, s_aux)
    r <- (nu/sqrt(2*pi))*(1/(s_aux^beta_aux))*zipfR::Igamma(beta_aux, s_aux)
    return(r)
  }

  acumu_aux <- function(q) stats::integrate(aux, lower=-Inf, upper=q)$value
  v.acumu_aux <- Vectorize(acumu_aux)
  cdf <- v.acumu_aux(q)
  cdf
}

#' @keywords internal
qslash <- function(p, nu){

  qtf <- function(input){
    p <- input

    if (!is.na(p)){
      obj <- function(q){
        pslash(q, nu) - p
      }

      nleqslv::nleqslv(stats::qnorm(p) / (0.5^(1/nu)), obj)$x
    }else{
      numeric(0)
    }
  }

  q0 <- rep(0, length(p[p == 0.5]) )
  q <- vector()
  if (length(p[p != 0.5]) < 1){
    q <- numeric(0)
  } else{
    q[p != 0.5 & p != 0 & p != 1] <-
      as.numeric(apply(matrix(p[p != 0.5 & p != 0 & p != 1], ncol = 1),
                       1, qtf))
  }

  q[p == 0] <- -Inf
  q[p == 1] <- Inf
  qtf <- c(q0, q)
  index <- c(which(p == 0.5), which(p != 0.5))

  qtf[sort(index, index.return = T)$ix]
}


#' Power Logit Distributions
#'
#' Density, distribution function, quantile function and random generation
#' for power logit distributions.

#' @name PL
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu vector of medians.
#' @param sigma vector of dispersion parameters.
#' @param lambda vector of skewness parameters.
#' @param zeta vector of extra parameters.
#' @param family string that specifies the family used to define the power logit distribution. The family is
#'     \code{NO}, \code{TF}, \code{LO}, \code{PE}, \code{SHN}, \code{Hyp} and \code{SLASH} for normal,
#'     Student-t, type II logistic, power exponential, sinh-normal, hyperbolic and slash
#'     distribution, respectively.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p). Default is FALSE.
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#'
#' @details If \code{zeta} is not specified, it assumes the default value 2. \cr
#'     The power logit distribution has density \deqn{f(y; \mu, \sigma, \lambda) = \lambda r(z^2)/(\sigma y (1-y^\lambda)),}
#'     for \eqn{y \in (0,1)}, in which \eqn{z = [logit(y^\lambda) - logit(\mu^\lambda)]/\sigma}, \eqn{r(\cdot)} is the density generator
#'     and \eqn{\mu \in (0,1)}, \eqn{\sigma>0} and \eqn{\lambda>0} are the median, dispersion and skewness of the distribution. \cr
#'     It is possible to consider \eqn{\lambda=0} to obtain the limiting case, the log-log distribution. This distribution has
#'     density \deqn{f(y; \mu, \sigma, \lambda) = r(z^2)/(\sigma y (-log(y))),} for \eqn{y \in (0,1)},
#'     in which \eqn{z = [-log(-log(y)) - (-log(-log(y)))]/\sigma}. \cr
#'     The \code{family} argument defines the density generator \eqn{r(\cdot)}, which may depend on an extra parameter (\code{zeta}).
#' @return \code{dPL} gives the density, \code{pPL} gives the distribution function, \code{qPL} gives the quantile function,
#'     and \code{rPL} generates random variables.
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022) \emph{The Power Logit Distribution}. ...
#' @examples
#' dPL(0.2, mu = 0.3, sigma = 1, lambda=1, zeta = 2, family = "PE")
#' mu = 0.3; sigma = 1; lambda = 2
#' set.seed(1)
#' PLsample = rPL(1000, mu, sigma, lambda, family = "SN",  zeta = 2.5)
#' hist(PLsample, prob = TRUE, breaks = 15, main = "", las = 1)
#' curve(dPL(x, mu, sigma, lambda, family = "SN",  zeta = 2.5),
#'       from = 0.01, to = 0.8, add = TRUE, col = "red")
#' rug(PLsample)
#'
#' x = seq(0.01, 0.9,0.01)
#' y = dPL(x, mu, sigma, lambda, family = "Hyp",  zeta = 2)
#' plot(x, y, type = "l", lwd = 2, las = 1)
#'
#' x1 = seq(0.01, 0.4, 0.01)
#' y1 = dPL(x1, mu, sigma, lambda, family = "Hyp",  zeta = 2)
#' polygon(c(x1, 0.4, 0), c(y1, 0, 0), col = "lightblue")
#' text(mu-0.025, 1, paste("P(Y<0.4) = ",
#'                         round(pPL(0.4, mu, sigma, lambda,
#'                         family = "Hyp",
#'                         zeta = 2),2)),
#'      font = 1, cex = 0.8)
#'
#' plot(x, pPL(x, mu, sigma, lambda, family = "PE",  zeta = 1.3),
#'      type = "l", las = 1, lwd = 2,
#'      ylab = expression(P(Y<y)),
#'      xlab = "y")
#' p = pPL(0.5, mu, sigma, lambda, family = "PE",  zeta = 1.3)
#' q = qPL(p, mu, sigma, lambda, family = "PE",  zeta = 1.3)
#' points(q, p, pch = 16, col = 2, cex = 1.5)
#' text(0.55, 0.83, paste("(", 0.5, ",", round(p, 2), ")"), font = 2,
#'      cex = 0.8, col = "red")
#' @importFrom stats dnorm qnorm rnorm pnorm dlogis plogis qlogis runif dt pt qt rbeta qbeta
#' @export
dPL <- function(x, mu, sigma, lambda, zeta = 2, family, log = FALSE){

  if (any(mu <= 0) | any(mu >= 1))
    stop(paste("mu must be beetwen 0 and 1 ", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(lambda < 0))
    stop(paste("lambda must be greater than or equal to 0", "\n", ""))
  if (any(zeta <= 0))
    stop(paste("zeta must be positive ", "\n", ""))
  if (any(x < 0) | any(x > 1))
    stop(paste("x must be beetwen (0, 1)", "\n", ""))

  if(lambda == 0){
    z <- (-1/sigma)*(log(log(x)/log(mu)))
  }else{
    z <- (1/sigma)*(VGAM::logitlink(x^lambda) - VGAM::logitlink(mu^lambda))
  }

  if(family == "NO"){
    r <- dnorm(z)
  }
  if(family == "TF"){
    r <- dt(z, zeta)
  }
  if(family == "LO"){
    r <- dlogis(z)
  }
  if(family == "SN"){
    r <- (2*cosh(sqrt(z^2))*exp(-2*sinh(sqrt(z^2))*sinh(sqrt(z^2))/zeta^2)/(sqrt(2*pi)*zeta))
  }
  if(family == "PE"){
    r <- gamlss.dist::dPE(z, mu = 0, sigma = 1, nu = zeta)
  }
  if(family == "Hyp"){
    r <- GeneralizedHyperbolic::dhyperb(z, mu = 0, delta = 1, alpha = zeta, beta = 0)
  }
  if(family == "SLASH"){
    s_aux <- z^2/2
    beta_aux <- zeta + (1/2)
    if(any(s_aux == 0)){
      s_aux[s_aux == 0] <- 0.0001
    }
    r <- (zeta/sqrt(2*pi))*(1/(s_aux^beta_aux))*zipfR::Igamma(beta_aux, s_aux)
  }

  if(lambda == 0){
    log.lik <- ifelse(x <= 0 | x >= 1, -Inf, - log(abs(sigma*x*log(x))) + log(r))
  }else{
    log.lik <- ifelse(x <= 0 | x >= 1, -Inf, log(lambda) - (log(sigma) + log(1 - x^lambda) + log(x) - log(r)))
  }

  if (log == FALSE)
    fy <- exp(log.lik)
  else fy <- log.lik
  fy
}

#' @rdname PL
#' @export
pPL <- function(q, mu, sigma, lambda, zeta = 2, family, lower.tail = TRUE, log.p = FALSE){


  if (any(mu <= 0) | any(mu >= 1))
    stop(paste("mu must be beetwen 0 and 1 ", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(lambda < 0))
    stop(paste("lambda must be greater than or equal to 0", "\n", ""))
  if (any(zeta <= 0))
    stop(paste("zeta must be positive ", "\n", ""))

  if(lambda == 0){
    z <- (-1/sigma)*(log(log(q)/log(mu)))
  }else{
    z <- (1/sigma)*(VGAM::logitlink(q^lambda) - VGAM::logitlink(mu^lambda))
  }

  if(family == "NO"){
    cdf <- pnorm(z, lower.tail = TRUE, log.p = FALSE)
  }
  if(family == "TF"){
    cdf <- pt(z, zeta, lower.tail = TRUE, log.p = FALSE)
  }
  if(family == "LO"){
    cdf <- plogis(z, lower.tail = TRUE, log.p = FALSE)
  }
  if(family == "SN"){
    cdf <-  pnorm((2/zeta)*sinh(z), lower.tail = TRUE, log.p = FALSE)
  }
  if(family == "PE"){
    cdf <- gamlss.dist::pPE(z, mu = 0, sigma = 1, nu = zeta, lower.tail = TRUE, log.p = FALSE)
  }
  if(family == "Hyp"){
    cdf <- GeneralizedHyperbolic::phyperb(z, mu = 0, delta = 1, alpha = zeta, beta = 0)
  }
  if(family == "SLASH"){
    cdf <- pslash(z, nu = zeta)
  }

  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}

#' @rdname PL
#' @export
qPL <- function(p, mu, sigma, lambda, zeta = 2, family, lower.tail = TRUE, log.p = FALSE){

  if (any(mu <= 0) | any(mu >= 1))
    stop(paste("mu must be beetwen 0 and 1 ", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(lambda < 0))
    stop(paste("lambda must be greater than or equal to 0", "\n", ""))
  if (any(zeta <= 0))
    stop(paste("zeta must be positive ", "\n", ""))

  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  if(family == "NO"){
    z_p <- qnorm(p, lower.tail = TRUE)
  }
  if(family == "TF"){
    z_p <- qt(p, zeta, lower.tail = TRUE)
  }
  if(family == "LO"){
    z_p <- qlogis(p, lower.tail = TRUE)
  }
  if(family == "SN"){
    z_p <- asinh((zeta/2)*qnorm(p, lower.tail = TRUE))
  }
  if(family == "PE"){
    z_p <- gamlss.dist::qPE(p, mu = 0, sigma = 1, nu = zeta, lower.tail = TRUE)
  }
  if(family == "Hyp"){
    z_p <- GeneralizedHyperbolic::qhyperb(p, mu = 0, delta = 1, alpha = zeta, beta = 0, lower.tail = TRUE)
  }
  if(family == "SLASH"){
    z_p <- qslash(p, nu = zeta)
  }

  if(lambda == 0){
    suppressWarnings(q <- exp(log(mu)*exp(-sigma*z_p)))
  }else{
    suppressWarnings(q <- ((mu^lambda)*exp(sigma*z_p)/(1 - (mu^lambda)*(1 - exp(sigma*z_p))))^(1/lambda))
  }
  q
}

#' @rdname PL
#' @export
rPL <- function(n, mu, sigma, lambda, zeta = 2, family){

  if (any(mu <= 0) | any(mu >= 1))
    stop(paste("mu must be beetwen 0 and 1 ", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(lambda < 0))
    stop(paste("lambda must be greater than or equal to 0", "\n", ""))
  if (any(zeta <= 0))
    stop(paste("zeta must be positive ", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)

  p <- runif(n)
  if(family == "SLASH"){
    z_p <- (qbeta(p, zeta, 1, lower.tail = TRUE)^(-1/2))*qnorm(p, lower.tail = TRUE)
    if(lambda == 0){
      suppressWarnings(r <- exp(log(mu)*exp(-sigma*z_p)))
    }else{
      suppressWarnings(r <- ((mu^lambda)*exp(sigma*z_p)/(1 - (mu^lambda)*(1 - exp(sigma*z_p))))^(1/lambda))
    }
  }else{
    r <- qPL(p, mu = mu, sigma = sigma, lambda = lambda, zeta = zeta, family = family)
  }
  r
}


