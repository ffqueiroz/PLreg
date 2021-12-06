
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PLreg

<!-- badges: start -->
<!-- badges: end -->

The **PLreg** package allows fitting power logit regression models.
Diagnostic tools associated with the fitted model, such as the
residuals, local influence measures, leverage measures, and goodness-of-
fit statistics, are implemented.

## Installation

You can install the development version of PLreg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ffqueiroz/PLreg")
```

## Main functions

#### <tt>dPL</tt>, <tt>pPL</tt>, <tt>qPL</tt>, and <tt>rPL</tt>

Currently, the **PLreg** package include 7 members of the power logit
class of distributions: the power logit normal, power logit Student-t,
power logit type II logistic, power logit power exponential, power logit
sinh-normal, power logit hyperbolic and power logit slash distributions.
The package provides the <tt>dPL</tt>, <tt>pPL</tt> and <tt>qPL</tt>
functions to compute the probability density function, cumulative
distribution function and quantile function of the power logit
distribution. Also, the <tt>rPL</tt> function may be used to generate
random sample of variables with power logit distribution. The basic
usages of these functions are:

``` r
dPL(x, mu, sigma, lambda, zeta = 2, family, log = FALSE)

pPL(q, mu, sigma, lambda, zeta = 2, family, lower.tail = TRUE, log.p = FALSE)

qPL(p, mu, sigma, lambda, zeta = 2, family, lower.tail = TRUE, log.p = FALSE)

rPL(n, mu, sigma, lambda, zeta = 2, family)
```

#### <tt>PLreg</tt>

The main function of the **PLreg** package is represented by
<tt>PLreg()</tt>, which allows to fit proportional data with power logit
regression model; this explains the name. The arguments of this function
are:

``` r
PLreg(formula, data, subset, na.action, family = c("NO", "LO", "TF", "PE", "SN", "SLASH", "Hyp"), 
      zeta = NULL, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"), 
      link.sigma = NULL, type = c("pML", "ML"), control = PLreg.control(...), 
      model = TRUE, y = TRUE, x = FALSE, ...)
```

The <tt>PLreg()</tt> function returns an object of class
“<tt>PLreg</tt>”, similar to “<tt>betareg</tt>” and “<tt>glm</tt>”
objects, which some methods available. The <tt>summary()</tt> method
presents a standard output, with coefficient estimates, standard errors,
partial Wald tests and p values for the regression coefficients, the
overall goodness-of-fit measure, the pseudo *R*<sup>2</sup>, etc.. The
<tt>type</tt> argument in <tt>summary()</tt> specifies the type of
residuals included in the output, currently three residuals are
supported: <tt>“standardized”, “quantile”</tt> and <tt>“deviance”</tt>.
The <tt>plot()</tt> method draws graphs for diagnostic and influence
analyses.

#### <tt>extra.parameter</tt>

An important function in the **PLreg** package is
<tt>extra.parameter()</tt>. This function can be used to estimate the
extra parameter of some power logit models. The basic usage is:

``` r
extra.parameter(object, lower, upper, grid = 10)
```

## Example

``` r
library(PLreg)
#> 
#> Attaching package: 'PLreg'
#> The following object is masked from 'package:stats':
#> 
#>     influence
## basic example code
```

In the following, a example is presented to illustrate the capacities of
**PLreg** package. We use the <tt>bodyfat\_Aeolus</tt> data set,
available in the package.

    help(bodyfat_Aeolus, package = "PLreg")

We start out with a model where <tt>percentfat</tt> depends on the sex
of the sample bat (<tt>sex</tt>), the hibernation time (<tt>days</tt>)
and the year that the bat was sample (<tt>year</tt>). After fitting some
models in power logit class, we choose the power logit power exponential
regression model with constant dispersion as the best model. In order to
fit this model, we first choose the best value for *ζ* by using the
<tt>extra.parameter()</tt> function. In this case, fit a power logit
model with any fixed zeta, e.g., *ζ* = 2, and then use this object in
<tt>extra.parameter()</tt> function as in the following.

``` r
fitPL_PE_start <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
              family = "PE", zeta = 2)
extra.parameter(fitPL_PE_start, lower = 1, upper = 2.5)
```

Then, fit the model with the choosed value.

``` r
fitPL_PE <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
              family = "PE", zeta = 1.7)
summary(fitPL_PE)
#> 
#> Call:
#> PLreg(formula = percentfat ~ days + sex + year, data = bodyfat_Aeolus, 
#>     family = "PE", zeta = 1.7)
#> 
#> Standardized residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.0006 -0.6000  0.0224  0.6353  2.5352 
#> 
#> Coefficients (median model with logit link):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.1615018  0.0647074 -17.950   <2e-16 ***
#> days        -0.0092638  0.0005363 -17.273   <2e-16 ***
#> sexM        -0.0499217  0.0570871  -0.874    0.382    
#> year2016     0.5176694  0.0607028   8.528   <2e-16 ***
#> 
#> Sigma coefficients (dispersion model with log link):
#>         Estimate Std. Error z value Pr(>|z|)    
#> (sigma) -1.94446    0.07869  -24.71   <2e-16 ***
#> 
#> Lambda coefficient:
#>           Estimate Std. Error
#> (lambda) 0.0004124      0.067
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Family: PL - PE ( 1.7 ) (Power logit power exponential)
#> Estimation method: pML (penalized maximum likelihood)
#> Log-likelihood:   320 on 6 Df
#> Pseudo R-squared: 0.6756
#> Upsilon statistic: 0.04852
#> AIC:  -628
#> Number of iterations in BFGS optimization: 20
```

The goodness of fit is assessed using diagnostic graphs through the plot
method.

``` r
plot(fitPL_PE, which = 1:4)
```

Further details and examples on the R package **PLreg** can be found
using the help on R, typing:

    help("PLreg")

## Reference

Queiroz, F.F. and Ferrari, S.L.P. (2022). Power logit regression for
modeling bounded data.
