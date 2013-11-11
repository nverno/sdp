## get a path graph showing all the sdp trajectories for the different plots
library(ggplot2)
source("~/work/functions/functions.R")
dfplot <- read.csv("~/work/data/data/long-bc-derived.csv")

## overlay graph of each plot over time
plot(log(dfplot$plotden), log(dfplot$plotbv), col="black",
         main = "Plot level self-thinning lines", pch = "")

showstuff <- ddply(dfplot, .(install, plot), .fun = function(x) {
    lines(log(x$plotden), log(x$plotbv), col="black",
         main = paste("Install:", mean(x$install),", Plot:", mean(x$plot),
         " BV vs Log Plot Density, Predicted and Observed"))
    params <- list(a1 = tl[tl$install==mean(x$install) & tl$plot==mean(x$plot),]$a1,
                   a0 = tl[tl$install==mean(x$install) & tl$plot==mean(x$plot),]$a0)
    abline(a = params[["a0"]], b = params[["a1"]])
})


## create dummy variable for group
library(grid)
dfplot <- dfplot[order(dfplot$install, dfplot$plot, dfplot$time),]
ggplot(dfplot, aes(log(plotden), log(plotbv), group=pplot, color = pplot)) +
    geom_path(arrow = arrow(), lwd = 2) +
    ggtitle("Plot level self-thinning lines for BC data") +
    xlab("Log plot density") + ylab("Log plot bole volume")
ggsave("~/work/sdp/visuals/plot-level-sdp-lines.pdf",width=16, height=8)
