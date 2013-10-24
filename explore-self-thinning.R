## Determine if a forest has reached the self-thinning line:
##  self-thinning line described by:
##
##   ln(Bt) = a0 + a1*ln(Nt)
##
## where Bt = bole volume, Nt = Stem density, t = census time
##
######
## This equation can describe the forest before and after reaching the self-thinning liine:
##
##     ln(Bt) = a0 + a1*ln(Nt) - [a0 + a1*ln(N0)-ln(B0)][N0/Nt]^a2
##
## where B0 and N0 are total bole volume and stem density at initial measurement
##
## data, in long format, includes dead trees
## NOTE: Use only Doug Firs
source("functions.R")


#### Lets see what the data looks like:
dflong <- read.csv("~/Allometry/data/long-dfonly-derived.csv")
dfplot <- read.csv("sdp-plot-level.csv")

## plot basal area, colored by plot number
dflong$pplot <- as.factor(dflong$pplot)
ggplot(dflong, aes(plotba, color = pplot)) + geom_histogram()

## plot basal area, colored by time
dflong$time = as.factor(dflong$time)
ggplot(dflong, aes(plotba, color = time)) + geom_histogram()
ggplot(dflong, aes(plotba)) + geom_histogram() + facet_wrap(~time)

## single plot example
single <- subset(dflong, pplot == 1.1)
ggplot(single, aes(plotba, color = time)) + geom_density()

## stem density
ggplot(dflong, aes(plotden)) + geom_histogram()
ggplot(dflong, aes(plotden, color=time)) + geom_histogram()
ggplot(dflong, aes(plotden, color=time)) + geom_density()
## single plot
ggplot(single, aes(plotden, color=time)) + geom_density()
ggplot(single, aes(plotden, color=time)) + geom_histogram()

## Looking at both variables together
ggplot(dflong, aes(plotden, plotba)) + geom_point()
ggplot(dflong, aes(plotden, plotba, color=time)) + geom_point() + geom_smooth(se = FALSE)

## individual plot trajectories
ggplot(dflong, aes(plotden, plotba, color=pplot, group=pplot)) + geom_point() +
    geom_path(arrow = arrow())

## log transform
ggplot(dflong, aes(log(plotden), log(plotba), color=pplot, group=pplot)) + geom_point() +
    geom_path(arrow = arrow())

ggplot(dflong, aes(log(plotden), log(plotbv), color=pplot, group=pplot)) + geom_point() +
    geom_path(arrow = arrow())

## Fitting to the data
## linear self-thinning line, fisrt with one plot
single <- subset(dfplot, install==2 & plot == 6 & !is.na(pplotden) & !is.na(pplotbv))
plot(single$plotden, single$plotbv)

thinning <- function(a0, a1, plotden) {
    a0 + a1*log(plotden)
}

thin <- lm(log(plotbv) ~ log(plotden), data = single)
plot(log(single$plotden), log(single$plotbv))
abline(thin)

## fit with MLE, nonlinear
## non-linear model to incorporate before and after reaching thinning line
nlthin <- function(params, Nt, N0, B0)
{
    a0 = params["a0"]
    a1 = params["a1"]
    a2 = params["a2"]
    a0 + a1*log(Nt) - (a0 + a1*log(N0)-log(B0))*(N0/Nt)^a2
}

normNLL <- function(params, x) {
    sd = params[["sd"]]
    mu = nlthin(params, Nt = single$plotden, N0 = single$pplotden,
    B0 = single$pplotbv)
    -sum(dnorm(x, mean = mu, sd = sd, log = TRUE))
}

## optim fit
stime <- Sys.time()
fit <- optim(fn = normNLL,
             par = c(sd = 1, a0 = 15, a1 = -1.5, a2 = -100),
             x = log(single[,"plotbv"]), method = "Nelder-Mead",
             control = list(maxit = 100000))
optim.time <- Sys.time() - stime

## graph results
params = fit$par
pred <- nlthin(params, Nt = single$plotden, N0 = single$pplotden, B0 = single$pplotbv)
plot(log(single$plotden), pred, col="red", main="Log Plot BV vs Log Plot Density, \nPredicted and Observed")
points(log(single$plotden), log(single$plotbv), col = "blue")
lines(log(single$plotden), pred, col ="red")
lines(log(single$plotden), log(single$plotbv), col = "blue")
abline(a = fit$par["a0"], b = fit$par["a1"], col= "green")
legend("bottomleft", legend = c("predicted","observed","thin line"),
       col=c("red","blue","green"), lty = 1)


## NLS, using MLE params
## fitting to a single plot
levels(dfplot$pplot)
ins = 5
p = 8
st = list(a0 = 11, a1 = -.75, a2 = -85)
single <- subset(dfplot, install==ins & plot==p & !is.na(pplotden) &
                 time > 76 & !is.na(pplotbv))
plot(log(single$plotden), log(single$plotbv))
abline(lm(log(single$plotbv)~log(single$plotden)))

nlmod <- nls(log(plotbv) ~ a0 + a1*log(plotden) - (a0 + a1*log(pplotden) - log(pplotbv))*
             (pplotden/plotden)^a2, data = single,
             start = st,
             control = nls.control(warnOnly=TRUE, maxiter = 1000))

## graph results
params = coef(nlmod)
pred <- nlthin(params, Nt = single$plotden, N0 = single$pplotden, B0 = single$pplotbv)
plot(log(single$plotden), log(single$plotbv), col = "black")
lines(log(single$plotden), pred, col="red")
##curve(thinning(coef(fitsingle)[1], coef(fitsingle)[2], x), add=TRUE)
abline(a = params["a0"], b = params["a1"], col = "green")

## Using Port algorithm
single <- subset(dfplot, install==2 & plot==11 & !is.na(pplotden) & !is.na(pplotbv))
lcols <- c("plotba","pplotba","plotbv","pplotbv")
tst <- single
tst[,lcols] <- apply(tst[,lcols],2,log)
err <- rnorm(n=nrow(tst),mean=-.051,sd = 0.01)
tst$lplotden <- log(tst$plotden)
tst$lpplotden <- log(tst$pplotden)
tst$lplotden <- tst$lplotden+err
tst$lpplotden <- tst$lpplotden + err
plot(tst$lplotden, tst$plotbv)

n1 <- nls(plotbv ~ a0 + a1*lplotden - (a0 + a1*lpplotden - pplotbv)*
          (pplotden/plotden)^a2, data = tst,
          start=c(a0=30,a1=-1.5,a2=-110),
          control = nls.control(warnOnly = TRUE))
abline(a = coef(n1)["a0"], b = coef(n1)["a1"], col = "green")
