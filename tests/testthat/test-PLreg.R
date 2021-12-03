#### IID case

mu <- 0.3; sigma <- 1; lambda <- 2
set.seed(1)
PLsample <- rPL(1000, mu, sigma, lambda, family = "SN",  zeta = 2.5)
hist(PLsample, prob = TRUE, breaks = 15, main = "", las = 1)
curve(dPL(x, mu, sigma, lambda, family = "SN",  zeta = 2.5),
      from = 0.01, to = 0.8, add = TRUE, col = "red")
rug(PLsample)

x <- seq(0.01, 0.9,0.01)
y <- dPL(x, mu, sigma, lambda, family = "Hyp",  zeta = 2)
plot(x, y, type = "l", lwd = 2, las = 1)

x1 <- seq(0.01, 0.4, 0.01)
y1 <- dPL(x1, mu, sigma, lambda, family = "Hyp",  zeta = 2)
polygon(c(x1, 0.4, 0), c(y1, 0, 0), col = "lightblue")

text(mu-0.025, 1, paste("P(Y<0.4) = ",
                       round(pPL(0.4, mu, sigma, lambda, family = "Hyp",  zeta = 2),2)),
     font = 1, cex = 0.8)

plot(x, pPL(x, mu, sigma, lambda, family = "PE",  zeta = 1.3),
     type = "l", las = 1, lwd = 2, ylab = expression(P(Y<y)), xlab = "y")
p <- pPL(0.5, mu, sigma, lambda, family = "PE",  zeta = 1.3)
q <- qPL(p, mu, sigma, lambda, family = "PE",  zeta = 1.3)
points(q, p, pch = 16, col = 2, cex = 1.5)
text(0.55, 0.83, paste("(", 0.5, ",", round(p, 2), ")"), font = 2,
     cex = 0.8, col = "red")


#### IID case
data("bodyfat_Aeolus")

fitPL <- PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
               family = "SN", zeta = 1.6)
summary(fitPL)

extra.parameter(fitPL, lower = 1, upper = 4)


hist(bodyfat_Aeolus$percentfat, prob = TRUE)
curve(dPL(x, mu = fitPL$fitted.values[1],
          sigma = exp(fitPL$coefficients$dispersion),
          lambda = fitPL$coefficients$skewness,
          zeta = fitPL$zeta, family = fitPL$family),
      from = 0.001, to = 0.4, add = TRUE, col = "red")



plot(ecdf(bodyfat_Aeolus$percentfat),verticals = TRUE,  lwd = 2, col = "gray",
     main = "", ylab="Cumulative density function", xlab = "y")
stripchart(bodyfat_Aeolus$percentfat, add = TRUE, at = 0, col = "black", pch=3)
curve(pPL(x, mu = fitPL$fitted.values[1],
          sigma = exp(fitPL$coefficients$dispersion),
          lambda = fitPL$coefficients$skewness,
          zeta = fitPL$zeta, family = fitPL$family),
      from = 0.001, to = 0.4, add = TRUE, col = "red")




#### Errors


expect_error(PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
                   family = "normal",
                   zeta = 4))
expect_error(PLreg(percentfat ~ 1, data = bodyfat_Aeolus, family = "TF"))
expect_error(PLreg(percentfat ~ 1, family = "TF"))
expect_error(PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
               family = "TF", zeta = 0))
expect_error(summary(PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
               family = "TF", zeta = 2)))
fitPL <- PLreg(percentfat ~ 1, data = bodyfat_Aeolus,
               family = "NO")
expect_error(extra.parameter(fitPL, lower = 1, upper = 4))

#### Different link functions
data("PeruVotes")
fitPL <- PLreg(votes ~ HDI | HDI, data = PeruVotes,
               family = "TF", zeta = 5,
               link = "probit", link.sigma = "sqrt")
summary(fitPL)

#### Different residuals
data("PeruVotes")
fitPL <- PLreg(votes ~ HDI | HDI, data = PeruVotes,
               family = "TF", zeta = 5)

res_quantile <- residuals(fitPL, type = "quantile")
res_standardized <- residuals(fitPL, type = "standardized")

plot(res_standardized, pch = "+", ylim = c(-6, 6))
abline(h = -3, lty = 2)
abline(h = 3, lty = 2)

qqnorm(res_quantile)
qqline(res_quantile)

#### Body fat data
data("bodyfat_Aeolus")


#Initial model with zeta = 2
fit1 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
                     family = "PE", zeta = 2)
summary(fit1)
# Choosing the best value for zeta
extra.parameter(fit1, lower = 1, upper = 4, grid = 15)

# Using zeta = 1.7
fit2 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
              family = "PE", zeta = 1.7)
summary(fit2)

# Fixing lambda = 1
fit3 <- PLreg(percentfat ~ days + sex + year, data = bodyfat_Aeolus,
              family = "PE", zeta = 1.7, control = PLreg.control(lambda = 1))
summary(fit3)

# Comparing the AIC and Upsilon values between fit2 and fit3
AIC(fit2) < AIC(fit3)
fit2$Upsilon.zeta < fit3$Upsilon.zeta


#### Peru votes
data("PeruVotes")

fitPL <- PLreg(votes ~ HDI | HDI, data = PeruVotes,
               family = "TF", zeta = 5, control = PLreg.control(lambda = 1))
summary(fitPL)



plot(fitPL, type = "standardized")

#### Firm
data("Firm")

fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
               family = "SLASH", zeta = 2.13)
summary(fitPL)
extra.parameter(fitPL, lower = 1.2, upper = 4, grid = 10)
plot(fitPL, type = "standardized")
envelope(fitPL, type = "standardized")

fitPL_wo72 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-72,],
               family = "SLASH", zeta = 2.13)
fitPL_wo15 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-15,],
                    family = "SLASH", zeta = 2.13)
fitPL_wo16 <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-16,],
                    family = "SLASH", zeta = 2.13)

coef.mu      <- coef(fitPL)[1:3]
coef.mu_wo72 <- coef(fitPL_wo72)[1:3]
coef.mu_wo15 <- coef(fitPL_wo15)[1:3]
coef.mu_wo16 <- coef(fitPL_wo16)[1:3]

plot(Firm$indcost, Firm$firmcost, pch = "+", xlab = "indcost", ylab = "firmcost")
#identify(Firm$indcost, Firm$firmcost)
covariate <- matrix(c(rep.int(1, 1000), rep(median(Firm$sizelog), 1000),
                     seq(0, 1.22, length.out = 1000)),
                   ncol = 3)
lines(covariate[,3], as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu)),
      type = "l")
lines(covariate[,3], as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo72)),
      type = "l", lty = 2, col = "blue")
lines(covariate[,3], as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo15)),
      type = "l", lty = 3, col = "red")
lines(covariate[,3], as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo16)),
      type = "l", lty = 4, col = "green")
parameters <- c("pML",
             "pML w/o 72",
             "pML w/o 15",
             "pML w/o 16")
legend(x = 0.5,
       y = 0.8,
       legend = parameters,
       col = c("black", "blue", "red", "green"),
       lty = c(1, 2, 3, 4),
       cex = 0.6)



#### influence diagnostics
data("Firm")

fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
               family = "SLASH", zeta = 2.13)

influence_measures <- influence(fitPL, graph = FALSE)
plot(influence_measures$case.weights, type = "h", ylim = c(0,1))
plot(influence_measures$totalLI, type = "h", ylim = c(0,6))
plot(Firm$sizelog, influence_measures$GL, pch = "+")


### Plot function

data("Firm")

fitPL <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
               family = "NO", zeta = 2.13)
par(mfrow = c(3,3))
plot(fitPL, type = "standardized")
par(mfrow = c(1, 1))
