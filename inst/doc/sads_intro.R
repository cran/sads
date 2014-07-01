### R code from vignette source 'sads_intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: R setup
###################################################
options(width=60, continue=" ")


###################################################
### code chunk number 2: installation-github (eval = FALSE)
###################################################
## library(devtools)
## install_github('sads', 'piklprado')


###################################################
### code chunk number 3: load-sads-package
###################################################
library(sads)


###################################################
### code chunk number 4: Loading datasets
###################################################
data(moths)# William's moth data
data(ARN82.eB.apr77)# Arntz et al. benthos data


###################################################
### code chunk number 5: Tabulating species in octaves
###################################################
(moths.oc <- octav(moths))
(arn.oc <- octav(ARN82.eB.apr77))


###################################################
### code chunk number 6: Ploting-octaves
###################################################
plot(moths.oc)


###################################################
### code chunk number 7: Biomass-octave-plot
###################################################
plot(arn.oc)


###################################################
### code chunk number 8: Rank-abundance tables
###################################################
head(moths.rad <- rad(moths))
head(arn.rad <- rad(ARN82.eB.apr77))


###################################################
### code chunk number 9: radplot1
###################################################
plot(moths.rad, ylab="Number of individuals")


###################################################
### code chunk number 10: radplots
###################################################
plot(arn.rad, ylab="Biomass")


###################################################
### code chunk number 11: Fitting a logseries model
###################################################
(moths.ls <- fitsad(moths,'ls'))


###################################################
### code chunk number 12: Operations on fitsad object
###################################################
summary(moths.ls)
coef(moths.ls)
logLik(moths.ls)
AIC(moths.ls)


###################################################
### code chunk number 13: Profiling and intervals
###################################################
moths.ls.prf <- profile(moths.ls)
confint(moths.ls.prf) # conf intervals


###################################################
### code chunk number 14: Ploting-profiles
###################################################
par(mfrow=c(1,2))
plotprofmle(moths.ls.prf)# log-likelihood profile
plot(moths.ls.prf)# z-transformed profile
par(mfrow=c(1,1))


###################################################
### code chunk number 15: Plot-of-predicted-values
###################################################
par(mfrow=c(2,2))
plot(moths.ls)
par(mfrow=c(1,1))


###################################################
### code chunk number 16: Fitting two other models
###################################################
(moths.pl <- fitsad(x=moths,sad="poilog"))#default is zero-truncated
(moths.ln <- fitsad(x=moths,sad="lnorm", trunc=0.5)) # lognormal truncated at 0.5


###################################################
### code chunk number 17: Model selection table
###################################################
AICtab(moths.ls, moths.pl, moths.ln, base=TRUE)


###################################################
### code chunk number 18: Predicted values for octaves
###################################################
head(moths.ls.oc <- octavpred(moths.ls))
head(moths.pl.oc <- octavpred(moths.pl))
head(moths.ln.oc <- octavpred(moths.ln))


###################################################
### code chunk number 19: Octaves-plot
###################################################
plot(moths.oc)
lines(moths.ls.oc, col="blue")
lines(moths.pl.oc, col="red")
lines(moths.ln.oc, col="green")
legend("topright", 
       c("Logseries", "Poisson-lognormal", "Truncated lognormal"), 
       lty=1, col=c("blue","red", "green"))


###################################################
### code chunk number 20: Predicted values - radplots
###################################################
head(moths.ls.rad <- radpred(moths.ls)) 
head(moths.pl.rad <- radpred(moths.pl))
head(moths.ln.rad <- radpred(moths.ln))


###################################################
### code chunk number 21: Rad-plots
###################################################
plot(moths.rad)
lines(moths.ls.rad, col="blue")
lines(moths.pl.rad, col="red")
lines(moths.ln.rad, col="green")
legend("topright", 
       c("Logseries", "Poisson-lognormal", "Truncated lognormal"), 
       lty=1, col=c("blue","red", "green"))


###################################################
### code chunk number 22: rsad-example1
###################################################
set.seed(42)# fix random seed to make example reproducible
(samp1 <- rsad(S = 10, frac = 0.1, sad = "lnorm", zeroes=TRUE,
               ssize = 2, meanlog = 3, sdlog = 1.5))


###################################################
### code chunk number 23: rsad-example2
###################################################
(samp2 <- rsad(S = 100, frac=0.1, sad="lnorm", 
              meanlog=5, sdlog=2))


###################################################
### code chunk number 24: rsad-poilog-fit
###################################################
(samp2.pl <- fitsad(samp2, 'poilog'))
## checking correspondence of parameter mu
coef(samp2.pl)[1] - log(0.1)


###################################################
### code chunk number 25: rsad-repeated-samples
###################################################
results <- matrix(nrow=250,ncol=2)
for(i in 1:250){
    x <- rsad(S = 100, frac=0.1, sad="lnorm", 
              meanlog=5, sdlog=2)
    y <- fitsad(x, "poilog")
    results[i,] <- coef(y)
}
results[,1] <- results[,1]-log(0.1)


###################################################
### code chunk number 26: rsads_bias
###################################################
##Mean of estimates
apply(results,2,mean)
## relative bias
(c(5,2)-apply(results,2,mean))/c(5,2)


###################################################
### code chunk number 27: rsads_bias
###################################################
##Mean of estimates
apply(results,2,sd)
## relative precision
apply(results,2,sd)/apply(results,2,mean)


###################################################
### code chunk number 28: rsads-bias-plots
###################################################
par(mfrow=c(1,2))
plot(density(results[,1]-log(0.1)))
abline(v=c(mean(results[,1]-log(0.1)),5), col=2:3)
plot(density(results[,2]))
abline(v=c(mean(results[,2]), 2), col=2:3)
par(mfrow=c(1,1))


