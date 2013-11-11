## Calculate SDP values based on fitted thinning lines
##  SDP values are ratio of (est. log(PlotBV))/(obs. log(PlotBV))
##  - calculate an sdp value for each sample period
tl <- read.csv("~/work/data/data/sdp/thinline-fits.csv")
sdps <- ddply(dfplot, .(install, plot, time), .fun = function(x) {
    x <- droplevels(x)
    a0 <- tl[tl$install==x$install & tl$plot==x$plot,]$a0
    a1 <- tl[tl$install==x$install & tl$plot==x$plot,]$a1
    est <- a0 + a1*log(x$plotden)
    data.frame = c(install = x$install, plot = x$plot, time=x$time,
    sdp = log(x$plotbv)/est)
})

## create classes for fits
## 1 = >0.995
## 2 = >0.95 <=0.995
## 3 = <=0.95
sdps <- merge(sdps, tl, all =TRUE)
quantile(sdps$sdp)
sdps$sdpclass <- cut(sdps$sdp, breaks = c(0,.95,.995,1.5))

## add sdps and sdpclass to dflong
dflong <- read.csv("~/work/data/data/long-dfonly-derived.csv")
times <- ddply(dflong, .(install,plot,time), .fun = function(x){
    x <- droplevels(x)
    data.frame = c(
    install = mean(x$install),
    plot = mean(x$plot),
    time = mean(x$time),
    num = nrow(x))
})
dflong$sdp <- rep(sdps$sdp, times = times$num)
dflong$sdpclass <- rep(sdpclass, times = times$num)

write.csv(dflong, "~/work/data/data/long-dfonly-derived.csv", row.names = FALSE)

