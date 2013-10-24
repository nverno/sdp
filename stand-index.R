## look at SI for doug fir data
source("functions.R")

dflong <- read.csv("long-dfonly-derived.csv")
dflong <- dflong[order(dflong$install, dflong$plot, dflong$time),]
dflong$si <- factor(dflong$si)
levels(dflong$sdpclass) <- c(3:1)

## SI vs growth
ggplot(dflong, aes(si, bvgrowth)) + geom_point()
ggplot(dflong, aes(si, bvgrowth, group = pplot)) + geom_boxplot()
ggplot(dflong, aes(si)) + geom_histogram()

## si and sdpclass
table(dflong$si)
table(dflong$sdpclass)

ggplot(dflong, aes(si)) + geom_freqpoly()

## bvgrowth vs interaction of si and sdpclass
dflong$sisdpclass <- interaction(dflong$si, dflong$sdpclass)
ggplot(dflong, aes(y = bvgrowth, x = sisdpclass)) + geom_boxplot()

dflong$ratio <- dflong$bvgrowth/dflong$priorbv
dflong <- dflong[dflong$ratio < 100,]
dfonly <- dflong[dflong$spec == "FD",]
ggplot(dfonly, aes(y = ratio, x = si, fill = sdpclass)) + geom_boxplot()

## sample sizes when fitting with sdpclass and si
bclong <- read.csv("long-bc-derived.csv")
samps <- ddply(bclong, .(sdpclass, si), function(x) nrow(x))
range(samps[,3])

