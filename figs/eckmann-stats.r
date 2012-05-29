
library(brem)
library(reshape2)
library(ggplot2)

load(paste("data/eckmann-small.rdata",sep=""))
load(paste("figs/eckmann-sm-fit.rdata",sep=""))

# Plot observed stats for total and a few pshifts for two blocks
df <- dyad.ps(train,N)
dimnames(df)[[3]] <- c("ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab","total")
df <- melt(df)
z <- fit$params$z  # use latest clustering
df$i <- z[df$Var1]
df$j <- z[df$Var2]
df <- subset(df,Var3 %in% c("total","ab-ba","ab-ay","ab-by"))
df$Var3 <- factor(as.character(df$Var3),c("total","ab-ba","ab-ay","ab-by"))
df <- subset(df,(i==1 & j==3) | (i==3 & j==1))
q3 <- qplot(Var3,value,data=df,geom="boxplot",outlier.size=0.5) + facet_wrap( ~ i+j) + theme_bw() + labs(x="",y="count for a given dyad") + scale_y_continuous(limits=c(0,35))

pdf("figs/eckmann-small/example-obs-stats.pdf",width=5,height=4)
q3
dev.off()

# Boxplot of posterior samples for each of the above blocks
effs <- c("intercept","ab-ba","ab-by","ab-xa","ab-xb","ab-ay","ab-ab")
b <- lapply(fit$samples,function(s) s$beta)
b <- melt(b)
b <- subset(b,Var1 %in% fit$priors$px & L1 > 20 &
            ((Var2 == 1 & Var3 == 3) | (Var2 ==3 & Var3 == 1)))
b$stat <- factor(effs[b$Var1],effs)

pdf("figs/eckmann-small/example-estimates.pdf",width=5,height=4)
qplot(stat,value,data=b,geom="boxplot",outlier.size=0.5)+facet_wrap(~Var2 + Var3) + theme_bw() + xlab("parameter")
dev.off()
