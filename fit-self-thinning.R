## Fit the self-thinning line to all plots, producing an
##  interactive graph during the process
##  Goal: use NLS fitting with Gauss-Newton (a more general version of
##   the algorithm Scott used, 'Marquardt') algorithm (closed-form)
##   wherever possible, if gauss-newton fails, try Nelder-Mead,
##   then SANN if all else fails...
source("functions.R")
dfplot <- read.csv("sdp-plot-level.csv")
## dfplot$pplot <- factor(dfplot$pplot)
## dfplot$install <- factor(dfplot$install)
## dfplot$plot <- factor(dfplot$plot)

## functions for fitting, nonlinear thin line and likelihood function
##  with normal distribution
nlthin <- function(params, Nt, N0, B0)
{
    a0 = params["a0"]
    a1 = params["a1"]
    a2 = params["a2"]
    a0 + a1*log(Nt) - (a0 + a1*log(N0)-log(B0))*(N0/Nt)^a2
}

normNLL <- function(params, x) {
    sd = params[["sd"]]
    mu = nlthin(params, Nt = get(x$plotden), N0 = get(x$pplotden), B0 = get(x$pplotbv))
    -sum(dnorm(x, mean = mu, sd = sd, log = TRUE))
}

## Make list of possible starting parameters
a0 <- seq(9,30,by=3)
a1 <- seq(-.25,-1.5, by = -.25)
a2 <- seq(-5000,2000,100)
params <- expand.grid(a0 = a0,a1 = a1,a2 = a2)
params <- params[params$a2 != 0,]

## Make a pass through all the plots using the same starting parameters,
##  store plots where the fit failed, where the estimates were bad (defined
##  by a1 > -.1) in two error lists: errors and badest.  Store data.frame with
##  successfully fitted plots and the three parameters.
dftest <- dfplot
nls.errors <- c() # stores plots that failed completely
badest <- c() # stores plots where a1 > -0.1 (i.e. bad fit)
thinline <- ddply(dftest, .(install, plot), .fun = function(x) {
    x <- droplevels(x)
    pars <- list(a0 = 11, a1 = -.75, a2 = -85)
    fit <- NULL
    try(fit <- nls(log(plotbv) ~ a0 + a1*log(plotden) -
                   (a0 + a1*log(pplotden) - log(pplotbv))*
                   (pplotden/plotden)^a2, data = x,
                   start = pars,
                   control = nls.control(warnOnly=TRUE)),silent = TRUE)
    if(is.null(fit) || coef(fit)[["a1"]] > -.3) {
        for(i in 1:nrow(params)) {
            if(i%%10==0)print(i) # track progress
            if(!is.null(fit) && coef(fit)[["a1"]] <= -0.3) break;
            if(is.null(fit) || coef(fit)[["a1"]] > -0.3) {
                pars <- list(a0=params[i,][["a0"]],a1=params[i,][["a1"]],
                             a2=params[i,][["a2"]])
                try(fit <- nls(log(plotbv) ~ a0 + a1*log(plotden) - ## round 2 with NLS
                               (a0 + a1*log(pplotden) - log(pplotbv))*
                               (pplotden/plotden)^a2, data = x,
                               start = pars,
                               control = nls.control(warnOnly=TRUE)),silent = TRUE)
            }
        }
    }
    if(is.null(fit)) { # if all NLS fits fail, resort to Nelder-Mead and log failure
        fit3 <- NULL
        nls.errors <<- c(nls.errors, levels(x$pplot))
        try( # try MLE fit -- Nelder-Mead algorithm
            fit3 <- optim(fn = normNLL,
                          par = c(sd = 1, a0 = 15, a1 = -1.5, a2 = -1000),
                          x = log(single[,"plotbv"]), method = "Nelder-Mead",
                          control = list(maxit = 100000)),silent = TRUE)
        if(!is.null(fit3)) {
            data.frame = c(
            install = mean(x$install),
            plot = mean(x$plot),
            a0 = fit3$par["a0"],
            a1 = fit3$par["a1"],
            a2 = fit3$par["a2"])
        }
    }
    if(!is.null(fit) & coef(fit)[["a1"]] > -0.1) { badest <<- c(badest, levels(x$pplot)) }
    if(!is.null(fit) & coef(fit)[["a1"]] <= -0.1) {
        data.frame = c(
        install = mean(x$install),
        plot = mean(x$plot),
        a0 = coef(fit)["a0"],
        a1 = coef(fit)["a1"],
        a2 = coef(fit)["a2"])
    }
    ##print(paste0("Install: ",mean(x$install)," ,Plot: ",mean(x$plot)))
})

## save fits data.frame
names(thinline) <- c("install","plot","a0","a1","a2")
write.csv(thinline, "thinline-fits.csv", row.names=FALSE)

## Lets look at the fits and the data
tl <- read.csv("thinline-fits.csv")
par(ask = TRUE)
showstuff <- ddply(dfplot, .(install, plot), .fun = function(x) {
    plot(log(x$plotden), log(x$plotbv), col="black",
         main = paste("Install:", mean(x$install),", Plot:", mean(x$plot),
         " BV vs Log Plot Density, Predicted and Observed"))
    params <- list(a1 = tl[tl$install==mean(x$install) & tl$plot==mean(x$plot),]$a1,
                   a0 = tl[tl$install==mean(x$install) & tl$plot==mean(x$plot),]$a0)
    abline(a = params[["a0"]], b = params[["a1"]])
})

## manually fit install 2, plot 11
##  This plot has a vertical line for the first 5 periods,
##  Going to just use a straight line fit to the last two periods
p11 <- subset(dfplot, install == 2 & plot == 11 & time > 72)
ggplot(p11, aes(log(plotden),log(plotbv))) + geom_point()
fit <- NULL
pars <- list(a0 = 20, a1 = -2, a2 = -10000)
fit <- nls(log(plotbv) ~ a0 + a1*log(plotden) -
           (a0 + a1*log(pplotden) - log(pplotbv))*
           (pplotden/plotden)^a2, data = p11,
           start = pars,
           control = nls.control(warnOnly=TRUE))


p11.fit <- lm(log(p11$plotbv) ~ log(p11$plotden))
tl[tl$install==2 & tl$plot==11,c("a0","a1","a2")] <- c(coef(p11.fit)[[1]],
                   coef(p11.fit)[[2]], NA)

## save fits data.frame again
write.csv(tl, "thinline-fits.csv", row.names=FALSE)




