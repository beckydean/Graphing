#####################################################
# Order the x axis of a graph according to ascending median values
#####################################################
# e.g. Figure 1, Dean et al 2018 Evolution
#####################################################
# Dummy data
values <- rnorm(1000, 0, 1)
categories <- rep(c('a', 'b', 'c', 'd', 'e'), 20)
data <- data.frame(values, categories)

# Get median values
mediandata <- aggregate(data$values, list(data$categories), median, na.rm=T)
names(mediandata) <- c("Category", "MedianValue")

# Order data with the median values
median.ord <- mediandata[order(mediandata$MedianValue, mediandata$Category), ]
median.ord$PlotOrder <- as.factor(seq(1:5))

# Merge with original data
dataOrder1 <- merge(data, median.ord, by.x = 'categories', by.y = 'Category')

# Plot
plot(values ~ PlotOrder, data = dataOrder1, ylab = 'Average value', xlab = 'Category in ascending order', xaxt = 'n', cex.lab = 1.5, col = 'grey90')
axis(1, at = seq(1:5), labels = median.ord$Category, las = 2)
mtext(side=3, line=1, adj=0, 'A')

#####################################################
# Binning an average value for each x axis incrememnt
#####################################################
# e.g Figure 3, 5 & 6 Dean et al 2018 Evolution
# And Figure 1 & 2 Dean & Mank 2016 American Naturalist
#####################################################
# Make dummy data
xval <- rnorm(1000, 0, 1)
yval <- rnorm(1000, 5, 10)
d <- data.frame(xval, yval)

minValue <- -2 #set the range of data along the x-axis
maxValue <- 2
step <- 0.2 #set the size of the bins - big number = fewer, smaller bins

#################
# the function binit gets mean value for each x-axis incremental increase
binit <- function(values, factors, minValue, maxValue, step) {
  breaks <- seq(minValue, maxValue - step, step)
  histcat <- as.factor(sapply(factors, function(x) { sum(x > breaks) }))
  
  newdf <- data.frame(binID = seq(1, length(breaks)), 
                      lower = breaks, 
                      midpoints = breaks + (step / 2))
  
  Responsemean <- tapply(values, histcat, mean) #can summarise each bin as mean or median
  Responselength <- tapply(values, histcat, length)
  Response.df <- data.frame(binID = rownames(Responsemean), 
                            ResponseMean = Responsemean, 
                            N = Responselength, 
                            freq = Responselength / sum(Responselength))
  
  finaldf <- merge(newdf, Response.df, by = "binID", all.x = TRUE)
  finaldf
}
#################

# run binit to get a dataframe of the binned values
df <- binit(d$yval, d$xval, minValue, maxValue, step)

# plot it
par(xpd=FALSE) # this means things can NOT be drawn outside the plotting region
plot(ResponseMean ~ midpoints, data = df, pch=21, 
     bg=rgb(0.1, 0.1, 0.1, 0.3), #red, green, blue, transparency/alpha
     col='black', xlab="", ylab="", 
     cex=log10(df$N), # size of datapoint based upon amount of data in each bin
     mgp=c(1,0.5,0), # (dist between axes titles and axes, dist between axis labels and axes,location of axes) default is c(3,1,0) although 
     tck = -0.02) 
abline(h = 2, lty=3) # horizontal line
abline(v = 0, lty=3) # vertical line
mtext("Sex bias", side = 1, line = 1.5, cex=1.2)
mtext("Spearman's rho", side = 2, line = 2, cex=1.2)
mtext("Male insemination capacity", line=0.2) # figure title
par(xpd=NA) # this means things CAN be drawn outside the plotting region
#someone once wanted me to label in shaded boxes the phenotype I was measuring - as if a figure legend was not good enough for him
rect(-2.05, 6.83, 2.05, 7.5, border='black', lwd=0.5, col=rgb(0.1, 0.1, 0.1, 0.2)) #order ing  rect is left size x-axis min, y-axis min, y-axis max, right side x-axis max 

# model the relationship between x and y axes
model <- lm(ResponseMean ~ midpoints , data=df, weights = freq, na.action = na.exclude)

# add model lines to graph
lines(seq(minValue, maxValue, step), predict(model, newdata = data.frame(midpoints=seq(minValue, maxValue, step))), lty=2, col="tomato3", lwd=3)

#####################################################
# Plotting interaction model fit in ggplots (this approach was suggested by Ben Bolker on stackoverflow)
#####################################################
# e.g. Figure 4 Dean et al 2018 Evolution
#####################################################
library(effects)
library(lme4)
library(ggplot2)
library(gridExtra)

# Is the expression level of male-biased genes more positively correlated with male fitness than that of female-biased genes? value = gene expression (continuous), SexBias = categorical

# run model for each fitness measure
mod1 <- lmer(Scale_LatCop ~ value * SexBias + (1|gene), data = All_Data1, control=lmerControl(optCtrl=list(maxfun=50000)))
mod2 <- lmer(Scale_Prom ~ value * SexBias + (1|gene), data = All_Data1, control=lmerControl(optCtrl=list(maxfun=50000)))

# get model fitted values
ee1 <- Effect(c("value","SexBias"),mod1) 
ee2 <- Effect(c("value","SexBias"),mod2)

# make plot pretty!
theme_set(theme_bw())

#The key is using as.data.frame() to turn the effects object into something useful
EE1 <- as.data.frame(ee1)
EE1$Measure <- "Inverse latency to copulate"
EE2 <- as.data.frame(ee2)
EE2$Measure <- "Male insemination capacity"

# use facet wrap to make one plot for the two fitness measures - make into one dataset
EE3 <- rbind(EE1, EE2)

# make a ggplot
p <- ggplot(EE3,
            aes(value,fit,colour=SexBias,fill=SexBias))+
  ylim(-0.08,0.08)+
  geom_line()+
  scale_color_manual(values=c("red", "blue", "green"))+ # colour of 3 levels of SexBias
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Gene Expression")+
  ylab("Fitted values")+
  theme(plot.title = element_text(hjust = 0)) +
  ## colour=NA suppresses edges of the ribbon
  geom_ribbon(colour=NA,alpha=0.1,
              aes(ymin=lower,ymax=upper)) +
  facet_wrap(~Measure, ncol=2)

# view the plot in new window
dev.new()
p
#####################################################
# Multi-panel figures with labelling
#####################################################
# e.g. Figure 3 & 4 Dean & Mank 2016 American Naturalist
#####################################################
# Dummy data
xval <- rnorm(500, 6, 12)
yval <- rnorm(500, 5, 10)
d <- data.frame(xval, yval)
###############################
# since we are plotting 6x6 times, here is a function for plotting
###############################
plot_it <- function(data){
  axislimit <- c(-50,50)
  
  # make an empty plot so the ablines can sit beneath the points
  plot(0,0, ylim = axislimit, xlim = axislimit, main = '', xaxt = 'n', yaxt = 'n')
  abline(0,1, col=alpha('black', 0.2), lwd=10) # thick, transparent 1:1 line
  abline(v=0, lty=2, col='grey')
  abline(h=0, lty=2, col='grey')
  points(data$xval, data$yval, col='black')
  
  # Add model line
  mod <- lm(data$yval ~ data$xval)
  abline(mod, lty=1, col='black', lwd=2)
  
  adjRsq <- summary(mod)[9][[1]]
  pval <- summary(mod)$coefficients[2,4]
  estimate <- summary(mod)$coefficients[2,1]
  text(-40, 40, signif(pval, digits=1), cex=1)
  
  }
###############################
# produce a 6 x 6 plot comparing data in two different species
dev.new(width=10.8, height=10.8)
par(mfrow=c(6,6),
    oma = c(2,1.5,0,0) + 0.7, #this changes the space around the whole plot
    mar = c(2,1,0.2,1) + 0.3) #this changes the size of the graphs
###############################
# set the parameters
spacing <- c(-50,-25,0,25,50)
speciesLabelSize <- 0.7
axisy <- 0.4 # distance between y-axis and numbering (horizontal adjust)
axisx <- -1.2 # this brings the labelling up close to the x axis
###############################
#row1
plot(0, xaxt='n', yaxt ='n', bty='n', pch='', ylab='', xlab='')
mtext("White Leghorn", line=0.1, cex=speciesLabelSize)
plot_it(d)
mtext("Black Minorcan", line=0.1, cex=speciesLabelSize)
plot_it(d)
mtext("Red Junglefowl", line=0.1, cex=speciesLabelSize)
plot_it(d)
mtext("Rhode Island Red", line=0.1, cex=speciesLabelSize)
plot_it(d)
mtext("Oxford Old English", line=0.1, cex=speciesLabelSize)
plot_it(d)
mtext("Yokahama", line=0.1, cex = speciesLabelSize)
mtext("White Leghorn", side=4, line=0.6, cex=speciesLabelSize)

#row2
plot_it(d)
axis(side=2, las=1, tck=-0.03, at=spacing, hadj=axisy, labels=spacing, cex.axis=1)
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
plot_it(d)
plot_it(d)
plot_it(d)
plot_it(d)
mtext("Black Minorcan", side=4, line=0.6, cex=speciesLabelSize)

#row3
plot_it(d)
axis(side=2, las=1, tck=-0.03, at=spacing, hadj=axisy, labels=spacing, cex.axis=1)
mtext("Change in female expression", side=2, line=1, outer=T, cex=1)
plot_it(d)
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
plot_it(d)
plot_it(d)
plot_it(d)
mtext("Red Junglefowl", side=4, line=0.6, cex=speciesLabelSize)

#row4
plot_it(d)
axis(side=2, las=1, tck=-0.03, at=spacing, hadj=axisy, labels=spacing, cex.axis=1)
plot_it(d)
plot_it(d)
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
plot_it(d)
plot_it(d)
mtext("Rhode Island Red", side=4, line=0.6, cex=speciesLabelSize)

#row5
plot_it(d)
axis(side=2, las=1, tck=-0.03, at=spacing, hadj=axisy, labels=spacing, cex.axis=1)
plot_it(d)
plot_it(d)
plot_it(d)
plot(0, xaxt='n', yaxt = 'n', bty='n', pch='', ylab='', xlab='')
plot_it(d)
mtext("Oxford Old English", side=4, line=0.6, cex=speciesLabelSize)

#row6
plot_it(d)
axis(side=2, las=1, tck=-0.03, at=spacing, hadj=axisy, labels=spacing, cex.axis=1)
axis(side=1, las=1, tck=-0.03, at=spacing, padj=axisx, labels=spacing, cex.axis=1)
mtext("Change in male expression", side=1, line=-0.5, outer=T, cex=1)
plot_it(d)
axis(side=1, las=1, tck=-0.03, at=spacing, padj=axisx, labels=spacing, cex.axis=1)
plot_it(d)
axis(side=1, las=1, tck=-0.03, at=spacing, padj=axisx, labels=spacing, cex.axis=1)
plot_it(d)
axis(side=1, las=1, tck=-0.03, at=spacing, padj=axisx, labels=spacing, cex.axis=1)
plot_it(d)
axis(side=1, las=1, tck=-0.03, at=spacing, padj=axisx, labels=spacing, cex.axis=1)
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
mtext("Yokahama", side=4, line=0.6, cex=speciesLabelSize)

#####################################################
# Heat map
#####################################################
# e.g. Figure 2 Dean et al 2017 Molecular Ecology
#####################################################
# dumy data
w1 <- rnorm(100, 0, 1)
w2 <- rnorm(100, 1, 1)
w3 <- rnorm(100, 2, 1)
w4 <- rnorm(100, 3, 1)
w5 <- rnorm(100, 0, 2)
w6 <- rnorm(100, 1, 2)
w7 <- rnorm(100, 2, 2)
w8 <- rnorm(100, 3, 2)
wrassedata <- data.frame(w1, w2, w3, w4, w5, w6, w7, w8)

library(pvclust) # clustering tree
library(pheatmap) # heatmap

# produce a tree
results <- pvclust(wrassedata, method.hclust = 'complete', method.dist = 'euclidean', nboot = 1000) # method.hclust = 'average',
plot(results)
pvrect(results, alpha=0.95) # which branches collapse

# determine the colour for the heatmap
my_palette <- colorRampPalette(c("black", "grey30", "grey40", 'grey60', "white", 'yellow', "gold", 'gold2', "goldenrod1"))(n=100)

pheatmap(wrassedata, cluster_row=T, cluster_cols=T, show_rownames=F, show_colnames=T, clustering_method='average', color=my_palette, fontsize=13, legend=F, border_color=NA)
#####################################################
# Clustered box plot
#####################################################
# e.g. Figure 4-6 Dean et al 2017 Molecular Ecology
#####################################################
boxplot(value ~ Morph + each_quartile, data = data, 
        outline=F, notch=T, # notched with no outliers shown
        col=c('tomato3', nmcol, satcol, sncol),  
        boxwex=0.8, # make boxes thinner
        at=c(1:4, 6:9, 11:14, 16:19), # put a space of 1 between each cluster
        ylim=c(0,10), xaxt='n',
        ylab=expression("Average log"[2]*' expression'), # subscript
        xlab ='Female expression level quartile', cex.lab=1.2)
axis(side = 1, las=2, tck = -0.02, at = c(2.5,7.5,12.5,17.5), labels = FALSE, cex.axis = 0.8)
text(c(2.5,7.5,12.5,17.5), -0.03, labels=c('1st', '2nd', '3rd', '4th'), pos=1, xpd= TRUE, cex=1) 
# denoting dignificant differences between bars
text(2.5, 0.11, "**", cex=0.7)
text(3, 0.19, "***", cex=0.7)
lines(c(2,4), c(0.18,0.18))
#####################################################
# Polygons as confidence intervals
#####################################################
# e.g. Figure 1 Dean et al 2015 Molecular Biology & Evolution
#####################################################
library(scales)
MGAFG <- c(0.0223, 0.8350, 0.7900)
PCRFG <- c(0.0327, 0.7990, 0.7610)
NMEFG <- c(0.0444, 0.7950, 0.7400)
ACYFG <- c(0.0934, 0.7490, 0.6640)
APLFG <- c(0.0989, 0.7230, 0.6180)
FG <- rbind(MGAFG, PCRFG, NMEFG, ACYFG, APLFG)
FG=as.data.frame(FG)
colnames(FG) <- c('time','GA', 'GZ')

plot(FG$GA ~ FG$time, axes = FALSE,  
     col = 'black', pch = 20, cex=2, cex.lab = 1.6, ylim = c(0.5,0.85),
     xlab = 'Divergence',
     ylab = 'Rank order male gene expression') # 

polygon(c(MS$time, rev(MS$time)), c(UpperSpleenA_m$x, rev(LowerSpleenA_m$x)), col = alpha('grey', 0.3), border=NA)
polygon(c(MS$time, rev(MS$time)), c(UpperGonadA_m$x, rev(LowerGonadA_m$x)), col = alpha('grey', 0.3), border=NA)
lines(FG$GA ~ FG$time, type = 'l', lty = 1, col = 'black', lwd = 2)
points(FG$GZ ~ FG$time, col = 'steelblue', pch = 20, cex=1.9)
lines(FG$GZ ~ FG$time, type = 'l', lty = 1, col = 'steelblue', lwd = 2)
mtext(side=3, line=1, adj=0, 'A')

box()
labs <- c('','MGA','PCR', 'NME', 'ACY', 'APL')
spacing <- c(0, 0.0223, 0.0327, 0.0444, 0.0934, 0.0989)
axis(side = 1, las=2, tck = -0.02, at = spacing, labels = FALSE, cex.axis = 2)

text(spacing, 0.48, labels=labs, srt = 45, pos=1, xpd= TRUE, cex=0.6)

axis(side = 2, at = c(0.5, 0.6,0.7, 0.8), tck = -0.02, las=1, cex.axis = 1.5, hadj=0.8)

#####################################################
# Violin plot with a jitter to show all the datapoints
#####################################################
# Figure 2 Dean & Pizzari in prep
#####################################################
library(ggplot2)
theme_set(theme_bw())
p <- ggplot(data, aes(x = categories, y = values)) +
  geom_violin(trim = TRUE) + 
  theme(text = element_text(size=15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Category") +
  ylab("Values") + 
  geom_point(color ='black', position = position_jitter(w=0.02)) #denotes how much jitter

# Plot violon + jitter + mean values
p + stat_summary(fun.data = mean_sdl, # mean plus/minus 1 standard deviation 
                 fun.args = list(mult=1), # denotes ± 1 S.D (default is 2)
                 geom = "pointrange", 
                 size = 0.8, 
                 col=c("tomato3", "steelblue", "tomato3", "steelblue", "tomato3")) 

# plots mean plus the data jitter 
p + stat_summary(fun.y = mean, 
                 geom = 'point', 
                 size = 4, 
                 col= c("tomato3", "steelblue","tomato3", "steelblue","tomato3"))
