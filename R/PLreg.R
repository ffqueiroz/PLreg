#' Auxiliary for Controlling PL Fitting
#'
#'Parameters that control fitting of power logit regression models using \code{\link{PLreg}}.
#' @name PLregcontrol
#' @param lambda numeric indicating the value of the skewness parameter lambda (if \code{NULL},
#'      lambda will be estimated).
#' @param method character specifying the \code{method} argument passed to \code{\link{optim}}.
#' @param maxit,trace,... arguments passed to \code{\link{optim}}
#' @param start an optional vector with starting values for median and dispersion submodels (starting value for lambda
#'      must not be included).
#' @details The \code{PLreg.control} controls the fitting process of power logit models. Almost all the arguments
#'     are passed on directly to \code{\link{optim}}, which is used to estimate the parameters.
#'     Starting values for median and dispersion submodels may be supplied via \code{start}. If the
#'     estimation process is to be performed with a fixed skewness parameter, a value must be specified
#'      in \code{lambda}. If \code{lambda = 0}, a log-log regression model
#'      will be estimated.
#'
#' @return A list with components named as the arguments.
#' @export
#' @seealso \code{\link{PLreg}}
#' @examples
#' data("PeruVotes")
#'
#'fitPL <- PLreg(votes ~ HDI | HDI, data = PeruVotes,
#'               family = "TF", zeta = 5, control = PLreg.control(lambda = 1))
#'summary(fitPL)
#'
PLreg.control <- function(lambda = NULL, method = "BFGS", maxit = 2000, trace = FALSE, start = NULL, ...){
  val <- list(lambda = lambda, method = method, maxit = maxit
              , trace = trace, start = start)
  val <- c(val, list(...))
  if (!is.null(val$fnscale))
    warning("fnscale must not be modified")
  val$fnscale <- -1
  if (is.null(val$reltol))
    val$reltol <- .Machine$double.eps^(1/2)
  val
}

#' @keywords internal
make.dmu.deta <- function(linkstr){
   switch(linkstr, "logit" = {
     logit_link <- make.link("logit")
     function(eta) logit_link$mu.eta(eta) * (1 - 2 * logit_link$linkinv(eta))
   },
   "probit" = function(eta) -eta * pmax(dnorm(eta), .Machine$double.eps),
   "cauchit" = function(eta) -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps),
   "cloglog" = function(eta) pmax((1 - exp(eta)) * exp(eta - exp(eta)), .Machine$double.eps),
   "loglog" = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
   "log" = function(eta) pmax(exp(eta), .Machine$double.eps),
   "sqrt" = function(eta) rep.int(2, length(eta)),
   "1/mu^2" = function(eta) 3/(4 * eta^2.5),
   "inverse" = function(eta) 2/(eta^3)
   )
}


#' Power Logit Regression Models for Bounded Variables
#'
#' \code{PLreg} is used to fit power logit regression model for continuous and bounded variables via maximum likelihood approach.
#' Both median and dispersion of the response variable are modeled through
#' parametric functions.
#' @name PLreg
#' @param formula a symbolic description of the model. See details for further information.
#' @param data,subset,na.action arguments controlling formula processing via \code{\link{model.frame}}.
#' @param family a description of the symmetric distribution to be used for generating the power logit model.
#'     Supported families include "\code{NO}", "\code{LO}", "\code{TF}", "\code{PE}", "\code{Hyp}", \code{SHN}"
#'      and "\code{SLASH}", which correspond to the power logit normal, type II logistic,
#'      Student-t, power exponential, hyperbolic, sinh-normal, and slash distributions, respectively.
#' @param zeta a numeric value or numeric vector that represents the extra parameter of the distribution. For the
#'     PL-NO and PL-LO models, no extra parameter is needed.
#' @param link an optional character that specifies the link function of the median submodel (mu).
#'     The "\code{logit}", "\code{probit}", "\code{cloglog}", "\code{cauchit}",
#'     "\code{loglog}" functions are supported. The \code{logit} function is the default.
#' @param link.sigma an optional character that specifies the link function of the dispersion submodel (sigma).
#'     The "\code{log}", "\code{sqrt}" functions are supported. The default is \code{log}.
#' @param type character specifying the type of estimator for the skewness parameter.
#'     Currently, penalized maximum likelihood ("\code{pML}") and maximum likelihood ("\code{ML}") are supported.
#'     If the skewness parameter is fixed, \code{ML} type is used.
#' @param control a list of control arguments specified via \code{\link{PLreg.control}}.
#' @param model,y,x logicals. If \code{TRUE} the corresponding components of the fit
#'     (model frame, response, model matrix) are returned.  For \code{\link{PLreg.fit}}, \code{y} must
#'     be the numeric response vector (with values in (0,1)).
#' @param X numeric regressor matrix for the median submodel.
#' @param S numeric regressor matrix for the dispersion submodel.
#' @param ... arguments passed to \code{\link{PLreg.control}}.
#'
#' @details The power logit regression models, proposed by Queiroz and Ferrari (2021), is useful in
#'     situations when the response variable is continuous and bounded on the unit interval (0, 1).
#'     The median and the dispersion parameters are modeled through parametric link
#'     functions. The models depend on a skewness parameter (called \eqn{\lambda}). When the skewness parameter is fixed
#'     and equal to 1, the power logit models coincide with the GJS regression models
#'     (Lemonte and Bazan, 2016). Queiroz and Ferrari (2021)  suggest using a penalized maximum
#'     likelihood method to estimate the parameters. This method is implemented in
#'     \code{PLreg} by default when \eqn{\lambda} is not fixed. If convergence is not reached,
#'     maximum likelihood estimation is performed. The estimation
#'     process uses \code{\link{optim}}. If no starting values are specified,
#'     the \code{PLreg} function uses those suggested by Queiroz and Ferrari (2021).
#'     This function also fits the log-log regression models by setting \eqn{\lambda}
#'     at zero (\eqn{\lambda = 0} represents \eqn{\lambda \rightarrow 0^+}).\cr \cr
#'     The formulation of the model has the same structure as in the usual functions
#'     \code{\link{lm}} and \code{\link{glm}}. The argument
#'     \code{formula} could comprise of three parts (separated by the symbols "\eqn{~}" and "\eqn{|}"),
#'     namely: observed response variable in the unit interval, predictor of the median submodel,
#'     with link function \code{link} and predictor of the dispersion submodel, with \code{link.sigma}
#'     link function. If the model has constant dispersion, the third part may be omitted and the link function for sigma
#'     is "\code{log}" by default. The skewness parameter \code{lambda} may be
#'     treated as fixed or not (default). If \code{lambda} is fixed, its value
#'     must be specified in the \code{control} argument. \cr \cr
#'     Some methods are available for objects of class "\code{PLreg}",
#'     see \code{\link{plot.PLreg}}, \code{\link{summary.PLreg}},
#'     \code{\link{coef.PLreg}}, \code{\link{vcov.PLreg}}, and
#'     \code{\link{residuals.PLreg}}, for details and other methods.
#'
#' @return \code{PLreg} returns an object of class "\code{PLreg}" with the following
#'     components (the \code{PLreg.fit} returns elements up to \code{v}).
#'   \item{coefficients}{a list with the "\code{median}", "\code{dispersion}" and
#'       "\code{skewness}" (if \code{lambda = NULL}) coefficients.}
#'   \item{residuals}{a vector of the raw residuals (the difference between the
#'       observed and the fitted response).}
#'   \item{fitted.values}{a vector with the fitted values of the median submodel.}
#'   \item{optim}{a list with the output from \code{optim}. When lambda is not fixed,
#'       if \code{type = "pML"}, the output refers to the iterative process of
#'       the median and dispersion parameters only and, if \code{type = "ML"},
#'        on the maximization of the likelihood for all the parameters.}
#'   \item{family}{a character specifying the \code{family} used.}
#'   \item{method}{the method argument passed to the optim call.}
#'   \item{control}{the control arguments passed to the optim call.}
#'   \item{start}{a vector with the starting values used in the iterative process.}
#'   \item{nobs}{number of observations.}
#'   \item{df.null}{residual degrees of freedom in the null model
#'       (constant median and dispersion), i.e., \eqn{n-3}.}
#'   \item{df.residual}{residual degrees of freedom in the fitted model.}
#'   \item{lambda}{value of the skewness parameter lambda
#'       (\code{NULL} when lambda is not fixed).}
#'   \item{loglik}{log-likelihood of the fitted model.}
#'   \item{vcov}{covariance matrix of all the parameters.}
#'   \item{pseudo.r.squared}{pseudo R-squared value.}
#'   \item{Upsilon.zeta}{an overall goodness-of-fit measure.}
#'   \item{link}{a list with elements "\code{median}" and "\code{dispersion}" containing the
#'       link objects for the respective models.}
#'   \item{converged}{logical indicating successful convergence of the
#'       iterative process.}
#'   \item{zeta}{a numeric specifying the value of zeta used in the estimation
#'       process.}
#'   \item{type}{a character specifying the estimation method used.}
#'   \item{v}{a vector with the v(z) values for all the observations (see Queiroz and
#'       Ferrari(2021)).}
#'   \item{call}{the original function call.}
#'   \item{formula}{the formula used.}
#'   \item{terms}{a list with elements "\code{median}", "\code{dispersion}" and "\code{full}" containing
#'       the term objects for the respective models.}
#'   \item{levels}{a list with elements "\code{median}", "\code{dispersion}" and "\code{full}" containing
#'       the levels of the categorical regressors.}
#'   \item{contrasts}{a list with elements "\code{median}" and "\code{dispersion}"
#'       containing the contrasts corresponding to levels from the respective models.}
#'   \item{model}{the full model frame (if \code{y = TRUE}).}
#'   \item{y}{the response variable (if \code{y = TRUE}).}
#'   \item{x}{a list with elements "\code{median}" and "\code{dispersion}" with the matrices from
#'       the median and dispersion submodels (if \code{x = TRUE}).}
#' @export
#' @author Francisco Felipe de Queiroz (\email{ffelipeq@@outlook.com}) and Silvia L. P. Ferrari.
#' @seealso \code{\link{summary.PLreg}}, \code{\link{PLreg.control}}, \code{\link{residuals.PLreg}}
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022). Power logit regression 
#'       for modeling bounded data. \emph{arXiv}:2202.01697. \cr \cr
#'       Lemonte, A. J. and Bazan, J. L. (2015). New class of Johnson SB distributions
#'       and its associated regression model for rates and proportions. \emph{Biometrical Journal}. 58:727-746.
#' @examples
#'#### Body fat data
#'data("bodyfat_Aeolus")
#'
#'#Initial model with zeta = 2
#'fit1 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
#'              family = "PE", zeta = 2)
#'summary(fit1)
#'# Choosing the best value for zeta
#'# extra.parameter(fit1, lower = 1, upper = 4, grid = 15)
#'
#'# Using zeta = 1.7
#'fit2 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
#'              family = "PE", zeta = 1.7)
#'summary(fit2)
#'
#'# Fixing lambda = 1
#'fit3 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
#'              family = "PE", zeta = 1.7,
#'              control = PLreg.control(lambda = 1))
#'summary(fit3)
#'
#'# Comparing the AIC and Upsilon values between fit2 and fit3
#'AIC(fit2) < AIC(fit3) # TRUE
#'fit2$Upsilon.zeta < fit3$Upsilon.zeta #TRUE
#'
#'#### Firm cost data
#'data("Firm")
#'
#'fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost,
#'               data = Firm,
#'               family = "SLASH",
#'               zeta = 2.13)
#'summary(fitPL)
#'#extra.parameter(fitPL, lower = 1.2, upper = 4, grid = 10)
#'#plot(fitPL, type = "standardized")
#'#envelope(fitPL, type = "standardized")
#'\donttest{
#'fitPL_wo72 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost,
#'                    data = Firm[-72,],
#'                    family = "SLASH",
#'                    zeta = 2.13)
#'fitPL_wo15 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost,
#'                    data = Firm[-15,],
#'                    family = "SLASH",
#'                    zeta = 2.13)
#'fitPL_wo16 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost,
#'                    data = Firm[-16,],
#'                    family = "SLASH",
#'                    zeta = 2.13)
#'
#'coef.mu      <- coef(fitPL)[1:3]
#'coef.mu_wo72 <- coef(fitPL_wo72)[1:3]
#'coef.mu_wo15 <- coef(fitPL_wo15)[1:3]
#'coef.mu_wo16 <- coef(fitPL_wo16)[1:3]
#'
#'plot(Firm$indcost, Firm$firmcost,
#'     pch = "+",
#'     xlab = "indcost",
#'     ylab = "firmcost")
#'#identify(Firm$indcost, Firm$firmcost)
#'covariate = matrix(c(rep.int(1, 1000),
#'                     rep(median(Firm$sizelog), 1000),
#'                     seq(0, 1.22, length.out = 1000)),
#'                   ncol = 3)
#'lines(covariate[,3],
#'      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu)),
#'      type = "l")
#'lines(covariate[,3],
#'      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo72)),
#'      type = "l", lty = 2, col = "blue")
#'lines(covariate[,3],
#'      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo15)),
#'      type = "l", lty = 3, col = "red")
#'lines(covariate[,3],
#'      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo16)),
#'      type = "l", lty = 4, col = "green")
#'parameters = c("pML",
#'               "pML w/o 72",
#'               "pML w/o 15",
#'               "pML w/o 16")
#'legend(x = 0.5,
#'       y = 0.8,
#'       legend = parameters,
#'       col = c("black", "blue", "red", "green"),
#'       lty = c(1, 2, 3, 4),
#'       cex = 0.6)}
#'
#' @importFrom stats .getXlevels as.formula cor delete.response make.link model.matrix
#'     model.response optim sd terms var dcauchy fitted median model.frame na.omit
#'     printCoefmat qqnorm quantile residuals
PLreg <- function(formula, data, subset, na.action,
                   family = c("NO", "LO", "TF", "PE", "SN", "SLASH", "Hyp"),
                   zeta = NULL, link = c("logit", "probit", "cloglog", "cauchit",
                            "loglog"), link.sigma = NULL,
                   type = c("pML", "ML"),
                   control = PLreg.control(...), model = TRUE, y = TRUE,
                   x = FALSE, ...)
{
  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2] < 2) {
    formula <- Formula::as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  } else {
    if (length(formula)[2] > 2) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts.", call. = FALSE) # RHS right-hand side
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  mtS <- delete.response(terms(formula, data = data, rhs = 2))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  S <- model.matrix(mtS, mf)
  if (length(Y) < 1)
    stop("Empty model", call. = FALSE)
  if (!(min(Y) > 0 & max(Y) < 1))
    stop("Invalid dependent variable, all observations must be in (0, 1).", call. = FALSE)
  n <- length(Y)
  family <- match.arg(family)
  if (family == "SLASH" | family == "TF" | family ==
      "SN" | family == "Hyp" | family == "PE") {
    if (is.null(zeta)) {
      stop("For the family of distributions specified by the user an extra parameter is required.", call. = FALSE)
    } else {
      if (zeta <= 0) stop("Invalid extra parameter; zeta must be positive.", call. = FALSE)
    }
  }else zeta <- 2
  type <- match.arg(type)
  if (is.character(link))
    link <- match.arg(link)
  if (is.null(link.sigma))
    link.sigma <- "log"
  if (is.character(link.sigma))
    link.sigma <- match.arg(link.sigma, c("log", "sqrt"))
  val <- PLreg.fit(X = X, y = Y, S = S, zeta = zeta, family = family,
                   link = link, link.sigma = link.sigma, type = type, control = control)
  val$call <- cl
  val$formula <- oformula
  val$terms <- list(median = mtX, dispersion = mtS, full = mt)
  val$levels <- list(median = .getXlevels(mtX, mf),
                     dispersion = .getXlevels(mtS, mf),
                     full = .getXlevels(mt, mf))
  val$contrasts <- list(median = attr(X, "contrasts"), dispersion = attr(S,"contrasts"))
  if (model)
    val$model <- mf
  if (y)
    val$y <- Y
  if (x)
    val$x <- list(median = X, dispersion = S)
  class(val) <- "PLreg"
  return(val)
}

#' @rdname PLreg
PLreg.fit <- function(X, y, S = NULL, family, type = "pML", zeta = zeta, link = "logit",
                       link.sigma = "log", control = PLreg.control())
{
  n <- NROW(X)
  p <- NCOL(X)
  if (is.null(colnames(X))) {
    if(p == 1) {
      colnames(X) <- "(Intercept)"
    } else {
      colnames(X)[1] <- "(Intercept)"
      for(i in 2:p){
        colnames(X)[1] <- paste("V", i, sep="")
      }
    }
  }
  if (is.null(S)) {
    q <- 1
    S <- matrix(1, ncol = q, nrow = n)
    colnames(S) <- "(Intercept)"
    rownames(S) <- rownames(X)
    sigma_const <- TRUE
  } else {
    q <- NCOL(S)
    if (q < 1L) stop("dispersion regression needs to have at least one parameter", call. = FALSE)
    sigma_const <- (q == 1) && isTRUE(all.equal(as.vector(S[, 1]), rep.int(1, n)))
  }
  if (is.character(link)) {
    linkstr <- link
    if (linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      linkobj$dmu.deta <- make.dmu.deta(linkstr)
    } else {
      linkobj <- structure(list(linkfun = function(mu) -log(-log(mu)),
                                linkinv = function(eta) pmax(pmin(exp(-exp(-eta)),
                                                                  1 - .Machine$double.eps), .Machine$double.eps),
                                mu.eta = function(eta) {
                                  eta <- pmin(eta, 700)
                                  pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                                }, dmu.deta = function(eta) pmax(exp(-exp(-eta) -
                                                                       eta) * expm1(-eta), .Machine$double.eps),
                                valideta = function(eta) TRUE, name = "loglog"),
                           class = "link-glm")
    }
  } else {
    linkobj <- link
    linkstr <- link$name
    if (is.null(linkobj$dmu.deta))
      warning("link needs to provide dmu.deta component", call. = FALSE)
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta

  if (is.character(link.sigma)) {
    sigma_linkstr <- link.sigma
    sigma_linkobj <- make.link(sigma_linkstr)
    sigma_linkobj$dmu.deta <- make.dmu.deta(sigma_linkstr)
  } else {
    sigma_linkobj <- link.sigma
    sigma_linkstr <- link.sigma$name
    if (is.null(sigma_linkobj$dmu.deta))
      warning("link.sigma needs to provide dmu.deta component." , call. = FALSE)
  }

  sigma_linkfun  <- sigma_linkobj$linkfun
  sigma_linkinv  <- sigma_linkobj$linkinv
  sigma_mu.eta   <- sigma_linkobj$mu.eta
  sigma_dmu.deta <- sigma_linkobj$dmu.deta

  logity <- linkfun(y)
  ocontrol <- control

  lambda_fix <- control$lambda

  method <- control$method
  start <- control$start
  control$lambda <- control$method <- control$start <- NULL

  type.used <- type

  if (is.null(start)) {
    beta  <- solve(t(X)%*%X)%*%t(X)%*%logity
    sigma <- rep.int(0, q)
    sigma[1L] <- sd(logity)
    if (!isTRUE(sigma_linkinv(sigma[1]) > 0)) {
      warning("No valid starting value for dispersion parameter found, using 1 instead", call. = FALSE)
      sigma[1L] <- 1
    }
    start <- list(median = beta, dispersion = sigma)
  }

  if (is.list(start)) start <- do.call("c", start)

  v.function   <- function(mu, sigma, lambda) {

    if(lambda == 0){
      z <- (-1/sigma)*(log(log(y)/log(mu)))
    }else{
      z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
    }

    if(family == "NO"){
      vz  <- rep(1, n)
      dvz <- rep(0, n)
    }
    if(family == "TF"){
      vz  <- (zeta + 1)/(zeta + z^2)
      dvz <- -2*(zeta + 1)*z/((zeta + z^2)^2)
    }
    if(family == "LO"){
      vz  <- (1 - exp(-abs(z)))/(abs(z)*(1 + exp(-abs(z))))
      dvz <- (2*abs(z)*exp(-abs(z)) + exp(-2*abs(z)) - 1)/(((z^2)^3/2)*((1 + exp(-abs(z)))^2))
    }
    if(family == "SN"){
      vz  <- 4*sinh(z)*cosh(z)/(zeta^2*z) - tanh(z)/z
      dvz <- ((cosh(z)*z - sinh(z))/z^2)*(4*cosh(z)/zeta^2 - 1/cosh(z)) +
        (sinh(z)/z)*(4*sinh(z)/zeta^2 + sinh(z)/(cosh(z)^2))
    }
    if(family == "PE"){
      pzeta <- sqrt(2^(-2/zeta)*gamma(1/zeta)*(gamma(3/zeta)^(-1)))
      vz    <- (zeta*(z^2)^(zeta/2 - 1))/(2*pzeta^zeta)
      dvz   <- ((zeta^2/2 - zeta)*(pzeta^(-zeta))*(z^2)^(zeta/2))/(z^3)
    }
    if(family == "Hyp"){
      vz  <- zeta/sqrt(1 + z^2)
      dvz <- -(z*zeta)/(1 + z^2)^(3/2)
    }
    if(family == "SLASH"){
      s_aux <- z^2/2
      beta_aux <- zeta + (1/2)
      gama_aux <- zipfR::Igamma(beta_aux, s_aux)
      vz  <- (2/z^2)*zipfR::Igamma(zeta + (3/2), s_aux)/gama_aux
      dvz <-  (-2/z)*vz + (2/z)*(1/gama_aux^2)*exp(-s_aux)*(s_aux^beta_aux)*(gama_aux*(1 - beta_aux/s_aux)
                                                                             + exp(-s_aux)*(s_aux^(beta_aux - 1)) )
    }
    list(z = z, v = vz, dv = dvz)
  }
  logL         <- function(theta){
    beta  <- theta[seq.int(length.out = p)]
    tau   <- theta[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)
    lambda <- exp(theta[seq.int(length.out = 1) + p + q])

    if(lambda == 0){
      z <- (-1/sigma)*(log(log(y)/log(mu)))
    }else{
      z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
    }


    if (any(!is.finite(z)))
      NaN
    else {
      ll <- suppressWarnings(dPL(y, mu = mu, sigma = sigma,
                                 lambda = lambda, zeta = zeta, family = family, log = TRUE))
      if (any(!is.finite(ll)))
        NaN
      else sum(ll)
    }
  }
  logLfixed    <- function(theta, lambda){

    beta  <- theta[seq.int(length.out = p)]
    tau   <- theta[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)

    if(lambda == 0){
      z <- (-1/sigma)*(log(log(y)/log(mu)))
    }else{
      z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
    }

    if (any(!is.finite(z)))
      NaN
    else {
      ll <- suppressWarnings(dPL(y, mu = mu, sigma = sigma,
                                 lambda = lambda, zeta = zeta, family = family, log = TRUE))
      if (any(!is.finite(ll)))
        NaN
      else sum(ll)
    }
  }
  U            <- function(theta){

    beta  <- theta[seq.int(length.out = p)]
    tau   <- theta[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)
    lambda <- exp(theta[seq.int(length.out = 1) + p + q])

    vdv <- v.function(mu, sigma, lambda)
    v   <- vdv$v
    z   <- vdv$z

    d1dot  <- as.vector(1/mu.eta(eta.1))
    d2dot  <- as.vector(1/sigma_mu.eta(eta.2))

    T1 <- diag(1/d1dot)
    T2 <- diag(1/d2dot)

    W <- diag(as.vector(z*v))

    mu_star     <- as.vector(lambda/(sigma*mu*(1 - mu^lambda)))
    sigma_star  <- as.vector((z^2*v - 1)/sigma)
    lambda_star <- as.vector((1/lambda) + ((y^lambda)*(log(y))/(1-y^lambda)) -
                               z*v*(1/sigma)*((log(y)/(1 - y^lambda)) -
                                                (log(mu)/(1 - mu^lambda))))

    U <- c(t(X)%*%W%*%T1%*%mu_star, t(S)%*%T2%*%sigma_star, t(rep.int(1, n))%*%lambda_star*lambda)
    U
  }
  Ufixed       <- function(theta, lambda){

    beta  <- theta[seq.int(length.out = p)]
    tau   <- theta[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)

    vdv <- v.function(mu, sigma, lambda)
    v   <- vdv$v
    z   <- vdv$z

    d1dot  <- as.vector(1/mu.eta(eta.1))
    d2dot  <- as.vector(1/sigma_mu.eta(eta.2))

    T1 <- diag(1/d1dot)
    T2 <- diag(1/d2dot)

    W <- diag(as.vector(z*v))
    if(lambda == 0){
      mu_star <- as.vector(-1/(sigma*mu*log(mu)))
    }else{
      mu_star <- as.vector(lambda/(sigma*mu*(1 - mu^lambda)))
    }
    sigma_star  <- as.vector((z^2*v - 1)/sigma)

    Ufixed <- c(t(X)%*%W%*%T1%*%mu_star, t(S)%*%T2%*%sigma_star)
    Ufixed
  }
  logLp        <- function(gama){
    aux <- optim(start, logLfixed, gr = Ufixed, lambda = exp(gama), method = method, control = control)
    theta.hat <- aux$par

    beta  <- theta.hat[seq.int(length.out = p)]
    tau   <- theta.hat[seq.int(length.out = q) + p]

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)
    lambda <- exp(gama)

    vdv <- v.function(mu, sigma, lambda)
    z   <- vdv$z
    v   <- vdv$v
    dv  <- vdv$dv
    z_lambda <- (1/sigma)*(log(y)/(1 - y^lambda) - log(mu)/(1 - mu^lambda))
    z_lambda.lambda <- (1/sigma)*(((y^lambda)*(log(y)^2)/(1 - y^lambda)^2) -
                                    ((mu^lambda)*(log(mu)^2)/(1 - mu^lambda)^2))
    c1 <- dv*z + v
    J_lambda.lambda <- - sum((-1/lambda^2) + ((y^lambda)*((log(y))^2)/((1-y^lambda)^2)) - v*z*z_lambda.lambda - (z_lambda^2)*c1)

    if (any(!is.finite(z)))
      NaN
    else {
      ll <- suppressWarnings(dPL(y, mu = mu, sigma = sigma,
                                 lambda = lambda, zeta = zeta, family = family, log = TRUE))
      if (any(!is.finite(ll)))
        NaN
      else sum(ll) + (1/2)*log(J_lambda.lambda)
    }
  }
  Jfunction    <- function(beta, tau, lambda){

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)

    d1dot  <- as.vector( 1/mu.eta(eta.1) )
    dd1dot <- as.vector( - dmu.deta(eta.1)*d1dot^3 )

    d2dot  <- as.vector(1/sigma_mu.eta(eta.2))
    dd2dot <- as.vector(-sigma_dmu.deta(eta.2)*d2dot^3)

    z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))

    vdv <- v.function(mu, sigma, lambda)
    v   <- vdv$v
    dv  <- vdv$dv

    mu_star     <- as.vector(lambda/(sigma*mu*(1 - mu^lambda)))
    sigma_star  <- as.vector((z^2*v - 1)/sigma)
    lambda_star <- as.vector((1/lambda) + ((y^lambda)*(log(y))/(1-y^lambda)) -
                               z*v*(1/sigma)*((log(y)/(1 - y^lambda)) -
                                                (log(mu)/(1 - mu^lambda))))
    w1 <- (mu_star^2/lambda)*(sigma*(1 - (mu^lambda)*(1 + lambda))*z*v +
                                + lambda*(v + z*dv))*(1/d1dot) +
      + mu_star*z*v*dd1dot/d1dot^2
    w2 <- -((1/sigma^2) - (z^2/sigma^2)*(3*v + z*dv))*(1/d2dot) +
      + sigma_star*(dd2dot/d2dot^2)
    w3 <- 1/lambda^2 - ((y^lambda)*((log(y))^2)/((1 - y^lambda)^2)) +
      + (1/sigma^2)*(v + z*dv)*((log(y)/(1 - y^lambda)) -
                                  (log(mu)/(1 - mu^lambda)))^2 +
      (1/sigma)*z*v*((((y^lambda)*(log(y)^2))/(1 - y^lambda)^2) -
                       (((mu^lambda)*(log(mu)^2))/(1 - mu^lambda)^2))
    w4 <- (mu_star/sigma)*z*(2*v + z*dv)
    w5 <- -mu_star*(((1 - (mu^lambda)*(1 - lambda*log(mu)))/(lambda*(1 - mu^lambda)))*z*v +
                      + (1/sigma)*((log(y)/(1 - y^lambda)) -
                                     (log(mu)/(1 - mu^lambda)))*(v + z*dv))
    w6 <- ( - 1/sigma^2)*((log(y)/(1 - y^lambda)) - (log(mu)/(1 - mu^lambda)))*z*(2*v + z*dv)

    T1 <- diag(1/d1dot)
    T2 <- diag(1/d2dot)

    W1 <- diag(as.vector(w1))
    W2 <- diag(as.vector(w2))
    W3 <- diag(as.vector(w3))
    W4 <- diag(as.vector(w4))
    W5 <- diag(as.vector(w5))
    W6 <- diag(as.vector(w6))

    Jbeta.beta     <- t(X)%*%W1%*%T1%*%X
    Jtau.tau       <- t(S)%*%W2%*%T2%*%S
    Jlambda.lambda <- t(rep.int(1, n))%*%W3%*%rep.int(1, n)
    Jbeta.tau      <- t(X)%*%W4%*%T1%*%T2%*%S
    Jbeta.lambda   <- t(X)%*%W5%*%T1%*%rep.int(1, n)
    Jtau.lambda    <- t(S)%*%W6%*%T1%*%rep.int(1, n)

    J <- rbind(cbind(Jbeta.beta, Jbeta.tau, Jbeta.lambda),
               cbind(t(Jbeta.tau), Jtau.tau, Jtau.lambda),
               cbind(t(Jbeta.lambda), t(Jtau.lambda), Jlambda.lambda))
    J
  }
  J0function    <- function(beta, tau){

    eta.1 <- as.vector(X %*% beta)
    eta.2 <- as.vector(S %*% tau)

    mu     <- linkinv(eta.1)
    sigma  <- sigma_linkinv(eta.2)

    d1dot  <- as.vector( 1/mu.eta(eta.1) )
    dd1dot <- as.vector( - dmu.deta(eta.1)*d1dot^3 )

    d2dot  <- as.vector(1/sigma_mu.eta(eta.2))
    dd2dot <- as.vector(-sigma_dmu.deta(eta.2)*d2dot^3)

    z <- (-1/sigma)*(log(log(y)/log(mu)))

    vdv <- v.function(mu, sigma, lambda)
    v   <- vdv$v
    dv  <- vdv$dv

    mu_star     <- as.vector(-1/(sigma*mu*log(mu)))
    sigma_star  <- as.vector((z^2*v - 1)/sigma)

    w1 <- - (mu_star^2)*(sigma*(log(mu) + 1)*z*v - v - z*dv)*(1/d1dot) + mu_star*z*v*dd1dot/d1dot^2
    w2 <- - ((1/sigma^2) - (z^2/sigma^2)*(3*v + z*dv))*(1/d2dot) + sigma_star*(dd2dot/d2dot^2)

    w3 <-  (mu_star/sigma)*z*(2*v + z*dv)

    T1 <- diag(1/d1dot)
    T2 <- diag(1/d2dot)

    W1 <- diag(as.vector(w1))
    W2 <- diag(as.vector(w2))
    W3 <- diag(as.vector(w3))

    Jbeta.beta     <- t(X)%*%W1%*%T1%*%X
    Jtau.tau       <- t(S)%*%W2%*%T2%*%S
    Jbeta.tau      <- t(X)%*%W3%*%T1%*%T2%*%S

    J <- rbind(cbind(Jbeta.beta, Jbeta.tau),
               cbind(t(Jbeta.tau), Jtau.tau))
    J
  }

  if(is.null(lambda_fix)){
    if(type == "ML"){
      theta.opt <- optim(par = c(start,0), fn = logL, gr = U, method = method, control = control, hessian = TRUE)
      beta   <- theta.opt$par[seq.int(length.out = p)]
      tau    <- theta.opt$par[seq.int(length.out = q) + p]
      lambda <- exp(theta.opt$par[seq.int(length.out = 1) + p + q])
      if (theta.opt$convergence == 0) {
        converged <- TRUE
      }else{
        converged <- FALSE
        warning("Optimization failed to converge.", call. = FALSE)
      }
      vcov <- solve(Jfunction(beta, tau, lambda))
    }else{
      gama.opt  <- tryCatch(optim( c(0), logLp, method = method, control = control), error=function(e) {e})
      fixed.opt <- tryCatch(optim(start[seq.int(length.out = p + q)], logLfixed, gr = Ufixed,
                                  method = method, lambda = exp(gama.opt$par), control = control), error=function(e) {e})
      if(BBmisc::is.error(gama.opt) | BBmisc::is.error(fixed.opt)){
        warning("Optimization failed to converge with type = pML, using type = ML instead.", call. = FALSE)
        theta.opt <- optim(par = c(start,0), fn = logL, gr = U, method = method, control = control)
        beta   <- theta.opt$par[seq.int(length.out = p)]
        tau    <- theta.opt$par[seq.int(length.out = q) + p]
        lambda <- exp(theta.opt$par[seq.int(length.out = 1) + p + q])
        if (theta.opt$convergence == 0) {
          converged <- TRUE
        }else{
          converged <- FALSE
          warning("Optimization failed to converge.", call. = FALSE)
        }
        type.used <- "ML"
        vcov <- solve(Jfunction(beta, tau, lambda))
      }else{
        lambda <- exp(gama.opt$par)
        beta   <- fixed.opt$par[seq.int(length.out = p)]
        tau    <- fixed.opt$par[seq.int(length.out = q) + p]
        if (gama.opt$convergence == 0 & fixed.opt$convergence == 0) {
          converged <- TRUE
        }else{
          converged <- FALSE
          warning("optimization failed to converge.", call. = FALSE)
        }
        vcov <- solve(Jfunction(beta, tau, lambda))
        if(any(diag(vcov) <= 0)){
          warning("Instability in the observed information matrix, using type = ML.", call. = FALSE)
          theta.opt <- optim(par = c(start,0), fn = logL, method = method, control = control)
          beta   <- theta.opt$par[seq.int(length.out = p)]
          tau    <- theta.opt$par[seq.int(length.out = q) + p]
          lambda <- exp(theta.opt$par[seq.int(length.out = 1) + p + q])
          if (theta.opt$convergence == 0) {
            converged <- TRUE
          }else{
            converged <- FALSE
            warning("Optimization failed to converge.", call. = FALSE)
          }
          type.used <- "ML"
          vcov <- solve(Jfunction(beta, tau, lambda))
        }
      }
    }
  }else{
    lambda    <- lambda_fix
    type.used <- "ML"
    if(lambda < 0 ) stop("invalid skewness parameter; lambda must be greater than or equal to 0.", call. = FALSE)
    theta.opt <- optim(start[seq.int(length.out = p + q)], logLfixed, gr = Ufixed,
                       method = method, lambda = lambda, control = control, hessian = TRUE)
    beta   <- theta.opt$par[seq.int(length.out = p)]
    tau    <- theta.opt$par[seq.int(length.out = q) + p]
    if (theta.opt$convergence == 0) {
      converged <- TRUE
    }else{
      converged <- FALSE
      warning("optimization failed to converge.", call. = FALSE)
    }
    if(lambda == 0){
      vcov <- solve(J0function(beta, tau))
    }else{
      vcov <- solve(Jfunction(beta, tau, lambda)[seq.int(length.out = p + q), seq.int(length.out = p + q)])
    }
  }
  eta.1 <- as.vector(X %*% beta)
  eta.2 <- as.vector(S %*% tau)
  mu     <- linkinv(eta.1)
  sigma  <- sigma_linkinv(eta.2)
  Upsilon <- function(zeta){
    fda <- sort(pPL(y, mu = mu, sigma = sigma, lambda = lambda, zeta = zeta, family = family))
    Upsilon_zeta <- mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
    Upsilon_zeta
  }

  optim.fit <- if(type.used == "pML") fixed.opt else theta.opt
  ll        <- logL(c(beta, tau, log(lambda)))

  Ups.zeta <- Upsilon(zeta)

  pseudor2 <- ifelse(var(eta.1) * var(linkfun(y)) <= 0, NA, cor(eta.1, linkfun(y))^2)
  v <- v.function(mu, sigma, lambda)$v
  names(beta) <- colnames(X)
  names(tau) <- if(sigma_const) "(sigma)" else colnames(S)
  names(lambda) <- "(lambda)"

  if(is.null(lambda_fix)){
    rownames(vcov) <- colnames(vcov) <- c(colnames(X), if(sigma_const) "(sigma)" else paste("(sigma)",
                                                                               colnames(S), sep = "_"), "(lambda)")
  }else{
    rownames(vcov) <- colnames(vcov) <- c(colnames(X), if(sigma_const) "(sigma)" else paste("(sigma)",
                                                                               colnames(S), sep = "_"))
  }

  val <- list(coefficients = list(median = beta, dispersion = tau, skewness = lambda),
              residuals = y - mu, fitted.values = structure(mu, .Names = names(y)),
              optim = optim.fit , family = family, method = method, control = ocontrol,
              start = start, nobs = n, df.null = n - ifelse(is.null(lambda_fix), 3, 2),
              df.residual = n - p - q - ifelse(is.null(lambda_fix), 1, 0),
              lambda = lambda_fix, loglik = ll, vcov = vcov,
              pseudo.r.squared = pseudor2, Upsilon.zeta = Ups.zeta,
              link = list(median = linkobj, dispersion = sigma_linkobj), converged = converged,
              zeta = zeta, type = type.used, v = v)
  return(val)
}


#' Methods for PLreg Objects
#'
#' Some S3 Methods for PLreg regression models.
#'
#' @name methodsPLreg
#' @param object,x fitted model object of class "\code{PLreg}".
#' @param type character specifying the type of residuals to be included in the
#'      summary output, see \code{\link{residuals.PLreg}}.
#' @param model character specifying for which component of the model the
#'      coefficients/covariance are extracted.
#' @param ... currently not used.
#' 
#' @details A set of methods for objects of class "\code{PLreg}", including methods
#'      for the functions \code{\link{summary}} and \code{\link{vcov}},
#'      which print the estimated coefficients along with some other information and
#'      presents the covariance matrix, respectively. The \code{summary} also presents
#'      the partial Wald tests for the model parameters. Finally, \code{summary} returns
#'      an object of class "\code{summary.PLreg}" containing information to be printed
#'      using the \code{print} method.
#'
#' @export
#' @seealso \code{\link{PLreg}}
#' @examples
#' data("bodyfat_Aeolus")
#'
#'fitPL <- PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
#'               family = "SN", zeta = 1.6)
#'fitPL
#'summary(fitPL)
#'coef(fitPL, model = "median")
#'vcov(fitPL)
#'logLik(fitPL)
#'AIC(fitPL)
summary.PLreg <- function(object, type = "standardized", ...)
{

  ## residuals
  type <- match.arg(type, c("quantile", "deviance", "standardized"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type

  ## extend coefficient table
  p <- length(object$coefficients$median)
  q <- length(object$coefficients$dispersion)
  cf <- as.vector(do.call("c", object$coefficients))
  cf <- if(is.null(object$lambda)) cf else cf[-(p + q + 1)]
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  cf <- if(is.null(object$lambda)) list(median = cf[seq.int(length.out = p), , drop = FALSE],
                                        dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE],
                                        skewness = cf[seq.int(length.out = 1) + p + q, 1:2, drop = FALSE]) else
                                          list(median = cf[seq.int(length.out = p), , drop = FALSE],
                                               dispersion = cf[seq.int(length.out = q) + p, , drop = FALSE])

  rownames(cf$median) <- names(object$coefficients$median)
  rownames(cf$dispersion) <- names(object$coefficients$dispersion)
  if(is.null(object$lambda)) rownames(cf$skewness) <- names(object$coefficients$skewness)
  object$coefficients <- cf

  ## number of iterations
  mytail <- function(x) x[length(x)]
  object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$count))))

  ## AIC
  object$AIC <- -2*object$loglik + 2*(p + q + ifelse(is.null(object$lambda), 1, 0))

  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.PLreg"
  object
}

#' @rdname methodsPLreg
#' @export
print.PLreg <- function(x, ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  digits <- 4
  if(!x$converged) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$median)) {
      cat(paste("Coefficients (median model with ", x$link$median$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$median, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in median model)\n\n")
    if(length(x$coefficients$dispersion)) {
      cat(paste("Sigma coefficients (dispersion model with ", x$link$dispersion$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$dispersion, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in dispersion model)\n\n")
    if(is.null(x$lambda)){
      cat(paste("Lambda coefficients:\n", sep = ""))
      print.default(format(x$coefficients$skewness, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }
  }
  invisible(x)
}

#' @rdname methodsPLreg
#' @export
print.summary.PLreg <- function(x, ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  digits <- 4
  if(!x$converged) {
    cat("Model did not converge\n")
  } else {
    types <- c("quantile", "deviance", "standardized")
    Types <- c("Quantile residuals", "Deviance residuals", "Standardized residuals")
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$median)) {
      cat(paste("\nCoefficients (median model with ", x$link$median$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$median, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in median model)\n")

    if(NROW(x$coefficients$dispersion)) {
      cat(paste("\nSigma coefficients (dispersion model with ", x$link$dispersion$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$dispersion, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in dispersion model)\n")

    if(is.null(x$lambda)) {
      cat(paste("\nLambda coefficient:\n", sep = ""))
      printCoefmat(x$coefficients$skewness, digits = digits, signif.legend = FALSE)
    } else {
      if(x$lambda == 0) {
        cat(paste("\nFixed skewness parameter (limiting case lambda -> 0).\n"))
      } else {
        cat(paste("\nFixed skewness parameter (lambda = ", x$lambda, ").\n", sep = ""))
      }
    }

    aux <- x$coefficients[c("median", "dispersion")]

    if(getOption("show.signif.stars") & any(do.call("rbind", aux)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    if(!is.null(x$lambda) && x$lambda == 0){
      if(x$family == "NO" | x$family == "LO"){
        cat("\nFamily: log-log -", x$family, switch(x$family,
                                                    "NO" = "(log-log normal)",
                                                    "LO" = "(log-log type II logistic)"))
      }else{
        cat("\nFamily: log-log -", x$family, "(", x$zeta, ")", switch(x$family,
                                                                      "TF"    = "(log-log Student-t)",
                                                                      "PE"    = "(log-log power exponential)",
                                                                      "Hyp"   = "(log-log hyperbolic)",
                                                                      "SLASH" = "(log-log slash)",
                                                                      "SN"    = "(log-log sinh-normal)"))
      }
    }else{
      if(x$family == "NO" | x$family == "LO"){
        cat("\nFamily: PL -", x$family, switch(x$family,
                                               "NO" = "(Power logit normal)",
                                               "LO" = "(Power logit type II logistic)"))
      }else{
        cat("\nFamily: PL -", x$family, "(", x$zeta, ")", switch(x$family,
                                                                 "TF"    = "(Power logit Student-t)",
                                                                 "PE"    = "(Power logit power exponential)",
                                                                 "Hyp"   = "(Power logit hyperbolic)",
                                                                 "SLASH" = "(Power logit slash)",
                                                                 "SN"    = "(Power logit sinh-normal)"))
      }
    }


    cat("\nEstimation method:", x$type, switch(x$type,
                                               "ML" = "(maximum likelihood)",
                                               "pML" = "(penalized maximum likelihood)"))
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    if(!is.na(x$Upsilon.zeta)) cat("\nUpsilon statistic:", formatC(x$Upsilon.zeta, digits = digits))
    if(!is.na(x$AIC)) cat("\nAIC:", formatC(x$AIC, digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
  }

  invisible(x)
}

#' @rdname methodsPLreg
#' @export
coef.PLreg <- function(object, ...) {
  cf <- object$coefficients
  name1 <- names(cf$median)
  name2 <- names(cf$dispersion)
  name3 <- names(cf$skewness)
  cf <- c(cf$median, cf$dispersion, cf$skewness)
  names(cf) <- c(name1, paste("(sigma)", name2, sep = "_"), name3)
  if(is.null(object$lambda)) cf else cf[-length(cf)]
}

#' @rdname methodsPLreg
#' @export
vcov.PLreg <- function(object, ...) {
  object$vcov
}

#' @rdname methodsPLreg
#' @export
logLik.PLreg <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}


#' @rdname methodsPLreg
#' @export
model.matrix.PLreg <- function(object, model = c("median", "dispersion"), ...) {
  model <- match.arg(model)
  val <- if(!is.null(object$x[[model]])) object$x[[model]]
  else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(val)
}

#' Residuals Method for PLreg Objects
#'
#' The function provides three types of residuals for power logit regression models:
#' quantile, deviance and standardized.
#'
#' @param object fitted model object of class "\code{PLreg}".
#'
#' @param type character specifying the type of residuals to be used.
#' @param ... currently not used.
#'
#' @details The \emph{quantile residuals} is based on Dunn and Smyth (1996) idea. The
#'     residuals are well-defined for all the distributions in the power logit class and
#'     have, approximately, a standard normal distribution in large samples if the model is
#'     correctly specified. \cr \cr
#'     The \emph{deviance residuals} are based on the log-likelihood contributions
#'     of each observation in the sample. The distribution of this residual is unknown
#'     (both exact and asymptotic),  except for the power logit normal model, which is,
#'     approximately, standard normal.\cr \cr
#'     The \emph{standardized residuals} are a standardized form of the ordinary
#'     residual. These residuals take into account the diagonal elements of the \eqn{H}
#'     matrix, being useful for detecting leverage observations. The distribution
#'     of the standardized residuals is unknown.
#'     
#' @return \code{residuals} method for object of class "\code{PLreg}" returns a vector
#'     with the residuals of the type specified in the \code{type} argument.
#' 
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022). Power logit regression 
#'       for modeling bounded data. \emph{arXiv}:2202.01697. \cr \cr
#'       Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
#'      \emph{Journal of Computational and Graphical Statistics}, 5:236-244.
#' @seealso \code{\link{PLreg}}, \code{\link{plot.PLreg}}, \code{\link{envelope}}, \code{\link{influence}}
#' @examples
#'data("PeruVotes")
#'fitPL <- PLreg(votes ~ HDI | HDI, data = PeruVotes,
#'               family = "TF", zeta = 5)
#'
#'res_quantile = residuals(fitPL, type = "quantile")
#'res_standardized = residuals(fitPL, type = "standardized")
#'
#'plot(res_standardized, pch = "+", ylim = c(-6, 6))
#'abline(h = -3, lty = 2)
#'abline(h = 3, lty = 2)
#'
#'qqnorm(res_quantile)
#'qqline(res_quantile)
#' @export
residuals.PLreg <- function(object,
                            type = c("quantile", "deviance", "standardized"), ...)
{
  type <- match.arg(type)
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  X <- if(is.null(object$x$median)) model.matrix(object, model = "median") else object$x$median
  S <- if(is.null(object$x$dispersion)) model.matrix(object, model = "dispersion") else object$x$dispersion
  mu <- object$fitted.values
  lambda <- object$coefficients$skewness
  sigma <- object$link$dispersion$linkinv(as.vector(S %*% object$coefficients$dispersion))
  family <- object$family
  zeta <- object$zeta

  res <- switch(type,
                "quantile" = {
                  qnorm(pPL(y, mu, sigma, lambda, zeta = zeta, family = family))
                },

                "deviance" = {
                  Li <- function(sat){
                    mu.dev <- if(sat == 1) y else mu
                    dPL(y, mu = mu.dev, sigma, lambda, zeta = zeta, family = family, log = TRUE)
                  }
                  sign(y - mu)*sqrt(2*abs(Li(sat = 1) - Li(sat = 0)))

                },

                "standardized" = {
                  if(family == "TF" & zeta <= 2) stop("The standardized residual is not well-defined for PL-t family with extra parameter less or equal to 2.", call. = FALSE)
                  if(family == "PE" & zeta < 1/2) stop("The standardized residual is not well-defined for PL-PE family with extra parameter < 0.5.", call. = FALSE)
                  if(family == "SLASH" & zeta < 1) stop("The standardized residual is not well-defined for PL-slash family with extra parameter < 1", call. = FALSE)
                  xidr.function <- function(mu, sigma, lambda) {

                    if(lambda == 0){
                      z <- (-1/sigma)*(log(log(y)/log(mu)))
                    }else{
                      z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
                    }

                    if(family == "NO"){
                      dr  <- 1
                      xih <- 1
                    }
                    if(family == "TF"){
                      dr  <- (zeta + 1)/(zeta + 3)
                      xih <- ifelse(zeta > 2, zeta/(zeta - 2), NA)
                    }
                    if(family == "LO"){
                      dr  <- 1/3
                      xih <- pi^2/3
                    }
                    if(family == "SN"){
                      dr   <- 2 + 4/(zeta^2) - (sqrt(2*pi)/zeta)*(1-2*(pnorm(sqrt(2)/zeta, mean = 0, sd = sqrt(2)/2)-0.5))*exp(2/(zeta^2))
                      dshn <- function(z) 2*cosh(z)*exp(-(2/zeta^2)*(sinh(z))^2)/(zeta*sqrt(2*pi))
                      fgf  <- function(z) dshn(z)*z^2
                      xih  <- 2*stats::integrate(fgf, 0, 20)$value
                    }
                    if(family == "PE"){
                      dr  <- ifelse(zeta > 1/2, (zeta^2)*(gamma(1/zeta)^(-2))*gamma(3/zeta)*gamma((2*zeta - 1)/zeta), NA)
                      xih <- 1
                    }
                    if(family == "Hyp"){
                      aux <-  function(x){
                        (sqrt(x^2 - 1)/x)*exp(-zeta*x)
                      }
                      h1  <- stats::integrate(f = aux, 1, Inf)$value
                      dr  <- (zeta^2)*h1/besselK(x = zeta, nu = 1, expon.scaled = FALSE)
                      xih <- besselK(x = zeta, nu = 2, expon.scaled = FALSE)/(zeta*besselK(x = zeta,
                                                                                           nu = 1, expon.scaled = FALSE))
                    }
                    if(family == "SLASH"){
                      dr  <- 4*(zeta*(zeta + 1/2)*((zeta + 3/2)*(zeta + 5/2) + zeta + 1))/((zeta+1)*((zeta
                                                                                                      + 3/2)^2)*(zeta + 5/2))
                      xih <- ifelse(zeta > 1, zeta/(zeta-1), NA)
                    }
                    list(dr = dr, xih = xih)
                  }

                  xidr <- xidr.function(mu, sigma, lambda)
                  xih <- xidr$xih
                  dr <- xidr$dr

                  T1 <- diag(object$link$median$mu.eta(as.vector(X %*% object$coefficients$median)))

                  if(lambda == 0){
                    z <- (-1/sigma)*(log(log(y)/log(mu)))
                    mu_star <- as.vector(-1/(mu*log(mu)))
                    D_star <- diag(mu_star)
                    D <- D_star%*%T1%*%X
                  }else{
                    z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
                    mu_star <- as.vector(1/(mu*(1-mu^lambda)))
                    D_star <- diag(mu_star)
                    D <- lambda*D_star%*%T1%*%X
                  }

                  Sigma <- diag(sigma^2)
                  Hd <- solve(Sigma^(1/2))%*%D%*%solve(t(D)%*%solve(Sigma)%*%D)%*%t(D)%*%solve(Sigma^(1/2))

                  z/(sqrt(xih)*sqrt(1 - ((dr*xih)^(-1))*diag(Hd)))
                })

  return(res)
}


#' Influence Diagnostics for PLreg Objects
#'
#' The \code{influence} function provides two influence measures and the generalized
#' leverage for power logit regression models.
#'
#' @param model fitted model object of class "\code{PLreg}".
#'
#' @param graph logical. If \code{graph = TRUE} the plots are shown, if
#'     \code{graph = FALSE} the plots are not shown. Default is \code{graph = TRUE}.
#' @param ... currently not used.
#' @return \code{influence} returns a list with three objects:
#'     \item{case.weights}{The values of \eqn{h_{max}} eigenvector based on case
#'     weights perturbation scheme (see Queiroz and Ferrari (2021)).}
#'     \item{totalLI}{The total local influence (see Lesaffre and Verbeke (1998))}
#'     \item{GL}{The diagonal elements of the generalized leverage matrix.}
#' @seealso \code{\link{PLreg}}, \code{\link{residuals.PLreg}}, \code{\link{envelope}},
#'     \code{\link{plot.PLreg}}
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022). Power logit regression 
#'       for modeling bounded data. \emph{arXiv}:2202.01697.
#' @examples
#'data("Firm")
#'
#'fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost,
#'               data = Firm, family = "SLASH", zeta = 2.13)
#'\donttest{
#'influence_measures = influence(fitPL, graph = FALSE)
#'plot(influence_measures$case.weights, type = "h", ylim = c(0,1))
#'plot(influence_measures$totalLI, type = "h", ylim = c(0,6))
#'plot(Firm$sizelog, influence_measures$GL, pch = "+")}
#'
#' @export
influence <- function(model, graph = TRUE, ...){

  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  X <- if(is.null(model$x$median)) model.matrix(model, model = "median") else model$x$median
  S <- if(is.null(model$x$dispersion)) model.matrix(model, model = "dispersion") else model$x$dispersion
  mu     <- model$fitted.values
  lambda <- model$coefficients$skewness
  sigma  <- model$link$dispersion$linkinv(as.vector(S %*% model$coefficients$dispersion))
  family <- model$family
  zeta   <- model$zeta
  J      <- solve(model$vcov)
  beta   <- model$coefficients$median
  tau    <- model$coefficients$dispersion

  p <- NCOL(X)
  q <- NCOL(S)
  n <- NROW(X)

  eta.1 <- as.vector(X %*% beta)
  eta.2 <- as.vector(S %*% tau)

  d1dot  <- as.vector( 1/model$link$median$mu.eta(eta.1) )
  dd1dot <- as.vector( - model$link$median$dmu.deta(eta.1)*d1dot^3 )

  d2dot  <- as.vector(1/model$link$dispersion$mu.eta(eta.2))
  dd2dot <- as.vector(-model$link$dispersion$dmu.deta(eta.2)*d2dot^3)

  if(lambda == 0){
    z       <- (-1/sigma)*(log(log(y)/log(mu)))
    mu_star <- as.vector(-1/(sigma*mu*log(mu)))
    Dy      <- diag(as.vector(-1/(sigma*y*log(y))))
  }else{
    z       <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
    mu_star <- as.vector(lambda/(sigma*mu*(1 - mu^lambda)))
    Dy      <- diag(as.vector(lambda/(sigma*y*(1 - y^lambda))))
  }
  v <- model$v

  sigma_star <- as.vector((z^2*v - 1)/sigma)

  W <- diag( as.vector(z*v) )

  D_beta <- diag(mu_star)
  D_tau  <- diag(sigma_star)

  T1 <- diag(1/d1dot)
  T2 <- diag(1/d2dot)

  Deltacw      <- rbind(t(X)%*%W%*%T1%*%D_beta, t(S)%*%T2%*%D_tau)
  case.weights <- abs(eigen(-t(Deltacw)%*%solve(-J[seq.int(length.out = p + q),seq.int(length.out = p + q)])%*%Deltacw)$vec[,1])

  totalLI <- NULL
  for(i in 1:n){
    totalLI[i] <- 2*abs(t(Deltacw[,i])%*%solve(-J[seq.int(length.out = p + q),seq.int(length.out = p + q)])%*%Deltacw[,i])
  }
  totalLI <- abs((totalLI - mean(totalLI))/sd(totalLI))

  # generalized leverage

  v.function <- function(mu, sigma, lambda) {

    if(lambda == 0){
      z <- (-1/sigma)*(log(log(y)/log(mu)))
    }else{
      z <- (1/sigma)*(VGAM::logitlink(y^lambda) - VGAM::logitlink(mu^lambda))
    }

    if(family == "NO"){
      vz  <- rep(1, n)
      dvz <- rep(0, n)
    }
    if(family == "TF"){
      vz  <- (zeta + 1)/(zeta + z^2)
      dvz <- -2*(zeta + 1)*z/((zeta + z^2)^2)
    }
    if(family == "LO"){
      vz  <- (1 - exp(-abs(z)))/(abs(z)*(1 + exp(-abs(z))))
      dvz <- (2*abs(z)*exp(-abs(z)) + exp(-2*abs(z)) - 1)/(((z^2)^3/2)*((1 + exp(-abs(z)))^2))
    }
    if(family == "SN"){
      vz  <- 4*sinh(z)*cosh(z)/(zeta^2*z) - tanh(z)/z
      dvz <- ((cosh(z)*z - sinh(z))/z^2)*(4*cosh(z)/zeta^2 - 1/cosh(z)) +
              (sinh(z)/z)*(4*sinh(z)/zeta^2 + sinh(z)/(cosh(z)^2))
    }
    if(family == "PE"){
      pzeta <- sqrt(2^(-2/zeta)*gamma(1/zeta)*(gamma(3/zeta)^(-1)))
      vz    <- (zeta*(z^2)^(zeta/2 - 1))/(2*pzeta^zeta)
      dvz   <- ((zeta^2/2 - zeta)*(pzeta^(-zeta))*(z^2)^(zeta/2))/(z^3)
    }
    if(family == "Hyp"){
      vz  <- zeta/sqrt(1 + z^2)
      dvz <- -(z*zeta)/(1 + z^2)^(3/2)
    }
    if(family == "SLASH"){
      s_aux <- z^2/2
      beta_aux <- zeta + (1/2)
      gama_aux <- zipfR::Igamma(beta_aux, s_aux)
      vz  <- (2/z^2)*zipfR::Igamma(zeta + (3/2), s_aux)/gama_aux
      dvz <-  (-2/z)*vz + (2/z)*(1/gama_aux^2)*exp(-s_aux)*(s_aux^beta_aux)*(gama_aux*(1 - beta_aux/s_aux)
                                                                             + exp(-s_aux)*(s_aux^(beta_aux - 1)) )
    }
    list(z = z, v = vz, dv = dvz)
  }
  v.dv <- v.function(mu, sigma, lambda)
  v <- v.dv$v
  dv <- v.dv$dv

  Wbeta   <- diag( c( v + z*dv ) )
  Wtau    <- diag( c( (1/sigma)*( 2*z*v + (z^2)*dv )  ) )
  if(dim(J)[1] == p + q){
    L_beta <- cbind(T1%*%X, matrix(rep.int(0, n*q), ncol = q))
    Lbetay <- rbind(t(X)%*%T1%*%Dy%*%Wbeta%*%D_beta,
                   t(S)%*%T2%*%Dy%*%Wtau)
    GL <- L_beta%*%solve(J)%*%Lbetay
  }else{
    Wlambda <- diag( (sigma*(y^lambda)*(1 + lambda*log(y) - y^lambda)/((1 - y^lambda)*lambda)) -
      ((1/sigma)*((log(y)/(1 - y^lambda)) - (log(mu))/(1 - mu^lambda))*(v + z*dv)) -
      (z*v*(((y^lambda)*(lambda*log(y) - 1) + 1)/((1 - y^lambda)*lambda))) )
    L_beta <- cbind(T1%*%X, matrix(rep.int(0, n*q), ncol = q), rep.int(0, n))
    Lbetay <- rbind(t(X)%*%T1%*%Dy%*%Wbeta%*%D_beta,
                   t(S)%*%T2%*%Dy%*%Wtau,
                   t(matrix(rep.int(1, n), ncol = 1))%*%Dy%*%Wlambda)
    GL <- L_beta%*%solve(J)%*%Lbetay
  }

  if(graph == TRUE){
    plot(case.weights, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Local influence")
    cat(paste("\nClick on points you would like to identify and press Esc."))
    identify(seq_len(n), case.weights, cex = 0.8)
    plot(totalLI, type = "h", main = "Case-weight perturbation", las = 1, xlab = "Index", ylab = "Total local influence")
    identify(seq_len(n), totalLI, cex = 0.8)
    plot(diag(GL), pch = "+", las = 1, cex = 0.8, main = "Generalized leverage", xlab = "Index", ylab = expression(GL[ii]))
    identify(seq_len(n), diag(GL), cex = 0.8)
  }

  list(case.weights = case.weights, totalLI = totalLI, GL = diag(GL))
}

#' Diagnostic Plots for PLreg Objects
#'
#' This function provides plots for diagnostic analysis of
#' power logit regression models.
#'
#' @param x fitted model object of class "\code{PLreg}".
#'
#' @param which numeric specifying a subset of plots (numbers between 1 and 7).
#'     Default is \code{which = 1:4}.
#' @param type character specifying the type of residuals,
#'     see \code{\link{residuals.PLreg}}. Default is \code{type = "standardized"}.
#' @param pch,las,cex,... graphical parameters (see \code{\link[graphics]{par}})
#'
#' @details The \code{plot} method for \code{\link{PLreg}} objects provides 7 types
#'     of diagnostic plots in the following order.
#'     \describe{
#'         \item{Residuals vs indexes of obs.}{An index plot of the residuals
#'             versus indexes of observations.}
#'         \item{Case-weight perturbation}{An index plot of local influence based on the
#'             case-weight perturbation scheme.}
#'         \item{Generalized leverage}{A dispersion diagram of the generalized leverage
#'             versus the predicted values.}
#'         \item{Residuals vs linear predictor}{A dispersion diagram of the residuals versus
#'             the linear predictors.}
#'         \item{Normal probability plot}{A normal probability plot of the residuals.}
#'         \item{Predicted vs observed values}{A dispersion diagram of the predicted values
#'             versus the observed values.}
#'         \item{Residuals vs v(z) function}{A dispersion diagram of the \eqn{v(z)} function
#'             versus the residuals. For some power logit models, the \eqn{v(z)} function
#'             may be interpreted as weights in the estimation process. If \code{family = "NO"},
#'             the \eqn{v(z)} function is constant.}
#'      }
#'      The \code{which} argument can be used to select a subset of the implemented plots.
#'      Default is \code{which = 1:4}.
#' @return \code{plot} method for \code{\link{PLreg}} objects returns 7 types
#'     of diagnostic plots. 
#' @importFrom graphics abline identify mtext par points text title
#' @seealso \code{\link{PLreg}}, \code{\link{residuals.PLreg}}, \code{\link{envelope}}, \code{\link{influence}}
#' @examples
#'data("Firm")
#'
#'fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
#'               family = "SLASH", zeta = 2.13)
#'par(mfrow = c(3,3))
#'plot(fitPL, type = "standardized")
#'par(mfrow = c(1, 1))
#' @export
plot.PLreg <- function(x, which = 1:4, type = "standardized", pch = "+",
                       las = 1, cex = 0.8, ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 7))
    stop("`which' must be in 1:7")

  main = ""
  types <- c("quantile", "deviance", "standardized")
  Types <- c("Quantile residuals", "Deviance residuals", "Standardized residuals")
  type <- match.arg(type, types)
  Type <- Types[type == types]

  sub.caption <- paste(deparse(x$call), collapse = "\n")

  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  res           <- residuals(x, type = type)
  influence.GL  <- influence(x, graph = FALSE)
  cw <- influence.GL$case.weights
  GL <- influence.GL$GL
  v  <- x$v

  n <- length(res)
  p <- length(x$coefficients$median)

  op <- par(ask = TRUE)
  on.exit(par(op))

  op2 <- par(mar=c(7,4,3,2)+.1)
  on.exit(par(op2))
  show <- rep(FALSE, 7)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1

  if(show[1]) {
     plot(1:n, res, xlab = "Index", ylab = Type, main = main, pch = pch, las = las, cex = cex)
     mtext("Residuals vs indexes of obs.", 3, 0.25)
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
     abline(h = 0, lty = 3, col = "gray")
  }

  if(show[2]) {
     plot(1:n, cw, xlab = "Obs. number", ylab = "Local influence", type = "h", main = main, las = las)
     mtext("Case-weight perturbation", 3, 0.25)
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
  }

  if(show[3]) {
     plot(fitted(x), GL, pch = pch, las = las, cex = cex,
          xlab = "Predicted values", ylab = "Generalized leverage")
     mtext("Generalized leverage", 3, 0.25)
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
  }
  if(show[4]) {
     plot(x$link$median$linkfun(x$fitted.values), res, xlab = "Linear predictor", ylab = Type,
          main = main, pch = pch, las = las, cex = cex, ...)
     mtext("Residuals vs linear predictor", 3, 0.25)
     abline(h = 0, lty = 3, col = "gray")
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
  }
  if(show[5]) {
     qqnorm(res, xlab = "Normal quantiles", ylab = Type,  pch = pch, las = las, cex = cex, main = main, ...)
     mtext("Normal probability plot of residuals", 3, 0.25)
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
  }
  if(show[6]) {
     plot(y, fitted(x),  xlab = "Observed values", ylab = "Predicted values", main = main,
          pch = pch, las = las, cex = cex, ...)
     mtext("Predicted vs observed values", 3, 0.25)
     abline(0, 1, lty = 2, col = "gray")
     if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
  }
  if(show[7]) {
    if(x$family == "NO"){
      warning("The v(z) function is constant for this family.")
      plot(res, v, xlab = Type, ylab = expression(v(z)), main = main,
           pch = pch, las = las, cex = cex, ...)
      mtext("v(z) function vs residuals", 3, 0.25)
      if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
    }else{
      plot(res, v, xlab = Type, ylab = expression(v(z)), main = main,
          pch = pch, las = las, cex = cex, ...)
      mtext("v(z) function vs residuals", 3, 0.25)
      if(one.fig) title(sub = sub.caption, cex.sub = 0.9, line = 5, ...)
    }
  }
  invisible()
}

#' Normal Probability Plots with Simulated Envelope of Residuals for PLreg Objects
#'
#' \code{envelope} is used to display normal probability plots with simulated
#' envelope of residuals for the power logit models. Currently, three types of
#' residuals are supported: quantile, deviance and standardized residuals.
#'
#' @param object fitted model object of class "\code{PLreg}".
#' @param type character specifying the type of residuals to be used,
#'     see \code{\link{residuals.PLreg}}. Default is \code{type = "standardized"}.
#' @param rep a positive integer representing the number of iterations to calculate
#'     the simulated envelopes. Default is \code{rep=40}.
#' @param conf a numeric value in the interval (0,1) that represents the confidence
#'     level of the simulated envelopes. Default is \code{conf=0.95}.
#' @param xlab character specifying the label for \eqn{x} axis (optional). Default is
#'     "\code{Quantile N(0,1)}".
#' @param ylab character specifying the label for \eqn{y} axis (optional). Default is
#'     the name of the used residual.
#' @param main character specifying the overall title for the plot.
#' @param ylim,xlim numeric values, specifying the left/lower limit and the right/upper
#'     limit of the scale.
#' @param envcol character specifying the color of the envelope.
#'
#' @details The \code{envelope} uses the idea of Atkinson (1985) to create normal
#'     probability plots with simulated envelope. Under the correct model,
#'     approximately 100*\code{conf} of the residuals are expected to be inside
#'     the envelope.
#' 
#' @return \code{envelope} returns normal probability plot with simulated envelopes
#'     for the residuals. 
#'
#' @export
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022). Power logit regression 
#'       for modeling bounded data. \emph{arXiv}:2202.01697. \cr \cr
#'      Atkinson, A. C. (1985) Plots, transformations and regression: an introduction
#'      to graphical methods of diagnostic regression analysis.
#'      \emph{Oxford Science Publications}, Oxford.
#' @importFrom methods missingArg
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{PLreg}}, \code{\link{residuals.PLreg}}
#' @examples
#' data("Firm")
#'
#'fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
#'               family = "SLASH", zeta = 2.13)
#'summary(fitPL)
#'\donttest{
#'envelope(fitPL, type = "standardized")
#'envelope(fitPL, type = "quantile")
#'envelope(fitPL, type = "deviance")}
#'
envelope <- function(object,type = c("quantile", "deviance", "standardized"),
                     rep = 40, conf = 0.95, xlab, ylab, main, envcol, ylim, xlim){
  type <- match.arg(type)
  rep <- max(30,floor(rep))

  if(rep != floor(rep) | rep <= 0) stop("The rep argument must be a positive integer.",call.=FALSE)
  if(conf <= 0  | conf >= 1) stop("The conf argument must be within the interval (0, 1).", call.=FALSE)

  family  <- object$family
  zetahat <- object$zeta
  n       <- object$nobs
  X       <- if(is.null(object$x$median)) model.matrix(object, model = "median") else object$x$median
  S       <- if(is.null(object$x$dispersion)) model.matrix(object, model = "dispersion") else object$x$dispersion
  muhat     <- object$fitted.values
  lambdahat <- object$coefficients$skewness
  sigmahat  <- object$link$dispersion$linkinv(as.vector(S %*% object$coefficients$dispersion))

  resRid <- residuals(object, type = type)

  resid_env <- matrix(0, n, rep)
  i <- 1
  bar <- txtProgressBar(min = 0, max = rep, initial = 0, width = 50, char = "+", style = 3)
  while(i <= rep){
    tryCatch({
      y_env <- rPL(n, muhat, sigmahat, lambdahat, zetahat, family = family)
      val <- suppressWarnings(PLreg.fit(X = X, y = y_env, S = S, zeta = zetahat, family = family,
                                        link = object$link$median, link.sigma = object$link$dispersion,
                                        type = object$type, control = object$control))

      mu     <- val$fitted.values
      sigma  <- val$link$dispersion$linkinv(as.vector(S %*% val$coefficients$dispersion))
      lambda <- val$coefficients$skewness
      zeta   <- val$zeta

      res_env <- switch(type,
                        "quantile" = {
                          qnorm(pPL(y_env, mu, sigma, lambda, zeta = zeta, family = family))
                        },

                        "deviance" = {
                          Li <- function(sat){
                            mu.dev <- if(sat == 1) y_env else mu
                            dPL(y_env, mu = mu.dev, sigma, lambda, zeta = zeta, family = family, log = TRUE)
                          }
                          sign(y_env - mu)*sqrt(2*abs(Li(sat = 1) - Li(sat = 0)))

                        },

                        "standardized" = {
                          if(family == "TF" & zeta <= 2) stop("The standardized residual is not well-defined for PL-t family with extra parameter less or equal to 2.", call. = FALSE)
                          if(family == "PE" & zeta < 1/2) stop("The standardized residual is not well-defined for PL-PE family with extra parameter < 0.5.", call. = FALSE)
                          if(family == "SLASH" & zeta < 1) stop("The standardized residual is not well-defined for PL-slash family with extra parameter < 1", call. = FALSE)
                          xidr.function <- function(mu, sigma, lambda) {

                            if(lambda == 0){
                              z <- (-1/sigma)*(log(log(y_env)/log(mu)))
                            }else{
                              z <- (1/sigma)*(VGAM::logitlink(y_env^lambda) - VGAM::logitlink(mu^lambda))
                            }

                            if(family == "NO"){
                              dr  <- 1
                              xih <- 1
                            }
                            if(family == "TF"){
                              dr  <- (zeta + 1)/(zeta + 3)
                              xih <- ifelse(zeta > 2, zeta/(zeta - 2), NA)
                            }
                            if(family == "LO"){
                              dr  <- 1/3
                              xih <- pi^2/3
                            }
                            if(family == "SN"){
                              dr   <- 2 + 4/(zeta^2) - (sqrt(2*pi)/zeta)*(1-2*(pnorm(sqrt(2)/zeta, mean = 0, sd = sqrt(2)/2)-0.5))*exp(2/(zeta^2))
                              dshn <- function(z) 2*cosh(z)*exp(-(2/zeta^2)*(sinh(z))^2)/(zeta*sqrt(2*pi))
                              fgf  <- function(z) dshn(z)*z^2
                              xih  <- 2*stats::integrate(fgf, 0, 20)$value
                            }
                            if(family == "PE"){
                              dr  <- ifelse(zeta > 1/2, (zeta^2)*(gamma(1/zeta)^(-2))*gamma(3/zeta)*gamma((2*zeta - 1)/zeta), NA)
                              xih <- 1
                            }
                            if(family == "Hyp"){
                              aux <-  function(x){
                                (sqrt(x^2 - 1)/x)*exp(-zeta*x)
                              }
                              h1  <- stats::integrate(f = aux, 1, Inf)$value
                              dr  <- (zeta^2)*h1/besselK(x = zeta, nu = 1, expon.scaled = FALSE)
                              xih <- besselK(x = zeta, nu = 2, expon.scaled = FALSE)/(zeta*besselK(x = zeta,
                                                                                                   nu = 1, expon.scaled = FALSE))
                            }
                            if(family == "SLASH"){
                              dr  <- 4*(zeta*(zeta + 1/2)*((zeta + 3/2)*(zeta + 5/2) + zeta + 1))/((zeta+1)*((zeta
                                                                                                              + 3/2)^2)*(zeta + 5/2))
                              xih <- ifelse(zeta > 1, zeta/(zeta-1), NA)
                            }
                            list(dr = dr, xih = xih)
                          }

                          xidr <- xidr.function(mu, sigma, lambda)
                          xih <- xidr$xih
                          dr <- xidr$dr

                          T1 <- diag(val$link$median$mu.eta(as.vector(X %*% val$coefficients$median)))

                          if(lambda == 0){
                            z <- (-1/sigma)*(log(log(y_env)/log(mu)))
                            mu_star <- as.vector(-1/(mu*log(mu)))
                            D_star <- diag(mu_star)
                            D <- D_star%*%T1%*%X
                          }else{
                            z <- (1/sigma)*(VGAM::logitlink(y_env^lambda) - VGAM::logitlink(mu^lambda))
                            mu_star <- as.vector(1/(mu*(1-mu^lambda)))
                            D_star <- diag(mu_star)
                            D <- lambda*D_star%*%T1%*%X
                          }

                          Sigma <- diag(sigma^2)
                          Hd <- solve(Sigma^(1/2))%*%D%*%solve(t(D)%*%solve(Sigma)%*%D)%*%t(D)%*%solve(Sigma^(1/2))

                          z/(sqrt(xih)*sqrt(1 - ((dr*xih)^(-1))*diag(Hd)))
                        })

      resid_env[,i] <- sort(res_env)
      setTxtProgressBar(bar,i)
      i = i + 1
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  liml <- apply(resid_env, 1, quantile, prob=(1-conf)/2)
  limu <- apply(resid_env, 1, quantile, prob=(1-(1-conf)/2))
  mm   <- apply(resid_env, 1, median)

  close(bar)
  cat("\n")
  types = c("quantile", "deviance", "standardized")
  Types <- c("Quantile residuals", "Deviance residuals", "Standardized residuals")
  faixaRid <- range(resRid , liml , limu)

  if(missingArg(xlab)) xlab <- "Quantile N(0,1)"
  if(missingArg(ylab)) ylab <- Types[types == match.arg(type, types)]
  if(missingArg(main)) main <- ""
  if(missingArg(envcol) || !is.character(envcol)) envcol <- "black"
  if(missingArg(ylim)) ylim <- faixaRid
  if(missingArg(xlim)){
    qqnormInt <- function(y,IDENTIFY = TRUE){
      qqnorm(y, pch = "+", las = 1, ylim = ylim, xlab = xlab, ylab = ylab, main = main, cex = 0.8,
             lwd = 3, las = 1) -> X
      cat(paste("\nClick on points you would like to identify and press Esc."))
      if(IDENTIFY) return(identify(X, cex = 0.8))
      invisible(X)
    }
    qqnormInt(resRid)


    par(new = TRUE)
    qqnorm(liml, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
           col = envcol)
    par(new = TRUE)
    qqnorm(limu, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 1, main = "",
           col = envcol)
    par(new = TRUE)
    qqnorm(mm, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, lty = 2, main = main)
  }else{
    qqnormInt <- function(y,IDENTIFY = TRUE){
      qqnorm(y, pch = "+", las = 1, ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, main = main, cex = 0.8,
             lwd = 3, las = 1) -> X
      cat(paste("\nClick on points you would like to identify and press Esc."))
      if(IDENTIFY) return(identify(X, cex = 0.8))
      invisible(X)
    }
    qqnormInt(resRid)

    par(new = TRUE)
    qqnorm(liml, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, xlim = xlim, lty = 1, main = "",
           col = envcol)
    par(new = TRUE)
    qqnorm(limu, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, xlim = xlim, lty = 1, main = "",
           col = envcol)
    par(new = TRUE)
    qqnorm(mm, axes = FALSE, xlab = "", ylab = "", type = "l", ylim = ylim, xlim = xlim, lty = 2, main = main)
  }

}


#' Procedure to Select the Extra Parameter for PLreg Objects
#'
#' The \code{extra.parameter} function is used to select the extra parameter
#' of some power logit models. It provides plots of -2\code{logLik} and the
#' Upsilon measure (see Queiroz and Ferrari (2021)) versus \eqn{\zeta},
#' the extra parameter.
#'
#' @param object fitted model object of class "\code{PLreg}".
#' @param lower a numeric value representing the lower limit of the interval for
#'     the extra parameter.
#' @param upper a numeric value representing the upper limit of the interval for
#'     the extra parameter.
#' @param grid a positive integer representing the number of points in the plots.
#'     Default is \code{grid=10}. If \code{grid} is less than 10, then \code{grid=10}.
#' @param graph logical. If \code{graph = TRUE} the plots are shown, if
#'     \code{graph = FALSE} the plots are not shown. Default is \code{graph = TRUE}.
#' @return \code{extra.parameter} returns a list with five objects:
#'     \item{zeta.Ups}{The selected zeta based on the Upsilon measure.}
#'     \item{zeta.loglik}{The selected zeta based on -2\code{logLik}.}
#'     \item{zeta.values}{The values of zeta used in the graphs.}
#'     \item{Upsilon.values}{-2\code{logLik} evaluated at each value of zeta.}
#'     \item{loglik.values}{Upsilon measure evaluated at each value of zeta.}
#' @seealso \code{\link{PLreg}}
#' @references Queiroz, F. F. and Ferrari, S. L. P. (2022). Power logit regression 
#'       for modeling bounded data. \emph{arXiv}:2202.01697.
#' @examples
#'data("bodyfat_Aeolus")
#'
#'#Initial model with zeta = 2
#'fit1 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
#'              family = "PE", zeta = 2)
#'summary(fit1)
#'# Choosing the best value for zeta
#'\donttest{
#'extra.parameter(fit1, lower = 1, upper = 4, grid = 15)}
#' @export
extra.parameter <- function(object, lower, upper, grid = 10, graph = TRUE){
  if(object$family == "NO" | object$family == "LO" ){
    stop("This model does not depend on extra parameters.")
  }
  new <- as.list(object$call)
  grid <- max(floor(abs(grid)),10)

  zetas <- seq(lower, upper, length = grid)
  conv  <- matrix(0,grid,1)
  Ups   <- matrix(0,grid,1)
  ll    <- matrix(0,grid,1)
  i <- 1
  bar <- txtProgressBar(min = 1, max = grid, initial = 0, width = 50, char = "+", style = 3)
  while(i <= grid){
    new$zeta <- zetas[i]
    val <- suppressWarnings(try(eval(as.call(new), envir = parent.frame()), silent = TRUE))
    if(is.list(val)){
      Ups[i] <- val$Upsilon.zeta
      ll[i]  <- -2*val$loglik
      conv[i] <- 1
    }
    i <- i + 1
    setTxtProgressBar(bar,i)
  }
  close(bar)
  cat("\n")
  zetas <- zetas[conv == 1]
  Ups   <- as.matrix(Ups[conv==1])
  loglik.zeta <- as.matrix(ll[conv == 1])

  if(graph == TRUE){
    op <- par(ask = TRUE)
    on.exit(par(op))

    plot(zetas, Ups, type = "l", xlim = range(zetas), ylim = range(Ups),
         xlab = expression(zeta), ylab = expression(Upsilon(zeta)), las = 1)
    mtext("Behaviour of Upsilon", 3, 0.25)
    points(zetas, Ups, pch = 16, cex = 0.6)
    abline(v = zetas[Ups == min(Ups)], lty = 2, col = "gray")
    txt <- bquote(zeta == .(paste(round(zetas[Ups == min(Ups)], 2))))
    text(zetas[Ups == min(Ups)] + (upper - lower)/grid, range(Ups)[2], txt   )

    plot(zetas, loglik.zeta, type = "l", xlim = range(zetas), ylim = range(loglik.zeta),
         xlab = expression(zeta), ylab = "-2*log-likelihood", las = 1)
    points(zetas, loglik.zeta, pch = 16, cex = 0.6)
    mtext("Behaviour of -2*log-Likelihood", 3, 0.25)
    abline(v = zetas[loglik.zeta == min(loglik.zeta)], lty = 2, col = "gray")
    txt <- bquote(zeta == .(paste(round(zetas[loglik.zeta == min(loglik.zeta)], 2))))
    text(zetas[loglik.zeta == min(loglik.zeta)] + (upper - lower)/grid, range(loglik.zeta)[2], txt)
  }

  cat(paste("\nEstimates for zeta are: \nzeta.Ups = ", round(zetas[Ups == min(Ups)], 2),
            "\nzeta.loglik = ", round(zetas[loglik.zeta == min(loglik.zeta)], 2), sep = ""))

  invisible(list(zeta.Ups = round(zetas[Ups == min(Ups)], 2),
              zeta.loglik = round(zetas[loglik.zeta == min(loglik.zeta)], 2),
              zeta.values = zetas, Upsilon.values = Ups,
              loglik.values = loglik.zeta))
}



