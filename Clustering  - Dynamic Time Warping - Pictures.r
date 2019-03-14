libraryrequiresetwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Hypervolume analyses")
nam <- as.vector(NULL)
for(ss in seq_along(dir(pattern = "_volumes.txt")[-c(7,9,11,14)]))
{
  n <- strsplit(dir(pattern = "_volumes.txt")[-c(7,9,11,14)][ss], split="_volumes.txt")[[1]][1]
  nam <- c(nam, n)
}

species_list <- nam
par(mfrow=c(2,2), mar=c(3,3,3,3))

### Cluster by all metrics together###

for (tt in seq_along(species_list))
{
  Vol <- read.delim(paste(species_list[tt], "_volumes.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(Vol))
  vol2 <- t(Vol)
  Hum_ov <-read.delim(paste(species_list[tt],"_humans_overlap.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(Hum_ov))
  over2 <- t(Hum_ov)
  comb1 <- cbind(Vol, Hum_ov)
  comb_tim <- as.vector(NULL)
  for (i in 1:nrow(comb1))
  {
    tim <- as.numeric(strsplit(strsplit(rownames(comb1)[i], split="_")[[1]][3], split="k")[[1]][1])
    comb_tim <- c(comb_tim, tim)
  }
  comb2 <- data.frame(Time=comb_tim, comb1)
  setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Geographic Overlap per time inerval/Geographic oVerlap 3 Variables")
  maxent_overlap <- read.delim(paste("Geographic_overlap_", species_list[tt], "_humans.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(maxent_overlap))
  colnames(maxent_overlap) <- c("Time", "Geogr_Overlap")
  comb3 <- merge(comb2, maxent_overlap, by.x="Time", by.y="Time", all.x=FALSE, all.y=TRUE, sort=TRUE)
  timeX <- comb3$Time
  comb4 <- t(comb3[,-1])
  colnames(comb4) <- timeX
  if(nrow(comb1) != nrow(maxent_overlap))
  {
    print(comb1)
  }
  assign(paste0(species_list[tt], "_timeseries", sep=""), comb4)
  setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Hypervolume analyses")
}

dist_len <- length(ls(pattern = "_timeseries"))
cl <- matrix(NA, nrow=dist_len, ncol=dist_len)
dimnames(cl) <- list(species_list, species_list)
for (i in seq_along(ls(pattern = "_timeseries")))
{
  for(j in seq_along(ls(pattern = "_timeseries")))
  {
    if (i!=j)
    {
      tim1 <- get(ls(pattern = "_timeseries")[i])
      tim2 <- get(ls(pattern = "_timeseries")[j])
      ## normalized distance matrix of all timeseries (as well different lengths) ###
      cl[i,j] <- dtw(t(tim1), t(tim2), distance.only = FALSE)$normalizedDistance 
    }
  }
}

cl2 <- as.dist(cl) ### it takes the lower triangular matrix
clusters <- hclust(cl2, method="complete") ## Hierarchical clustering
mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58") ## choose colors
op = par(bg = "#E8DDCB")
require("ape")
clus_cut <- cutree(clusters, k=c(2:5))
## Plot as phylogenetic tree - type fan ##
# par(mfrow=c(2,2), mar=c(2.5,1.8,2.6,1.6))
# for (cuts in 1:4)
# {
#   plot(as.phylo(clusters), type = "fan", tip.color = mypal[clus_cut[,cuts]])
#   # label.offset = 1)
# }
# 
par(mfrow=c(2,2), mar=c(3,3,3,3))
for (cuts in 1:4)
{
  cols <- as.numeric(colnames(clus_cut))[cuts]
  plot(clusters, labels=rownames(clus_cut), cex=0.6, hang=-1, col=mypal[1], lwd=2, col.axis=mypal[1], main="Species Clustering", cex.main=0.7, )
  rect.hclust(clusters, k=cols, border=mypal[cuts])
}

library(MASS)

## Non- Metric Multidimennsional Scaling ##
#mydata <- as.dist(cl) ### NMDS needs a distance matrix to represent the low dimension distances between species
# fit <- cmdscale(cl2, k=2) ### non metric MDS. k is the number of dimensions
# x <- fit$points[,1]
# y <- fit$points[,2]
# plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric	MDS",	type="n")
# text(x, y, labels = row.names(cl), cex=.7)

## Use pictures instead of points ##
library(jpeg)
library(png)
fit <- isoMDS(cl2, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

plot(x, y, xlab="", ylab="", main="Nonmetric MDS", cex.main=0.8, type="n", xlim=c(min(fit$points[,1])-1, max(fit$points[,1]+1)),
     ylim=c(min(fit$points[,2])-1, max(fit$points[,2]+1)), bty="l", axes=F, xaxs="i")
axis(side=1, col="black", at=seq(-7,11,2), lwd=1, col.ticks="black", cex.axis=0.7)
axis(side=2, col="black", at=seq(-7,7,2), lwd=1, col.ticks="black", cex.axis=0.7)
setwd("C:/Users/Dest/Desktop/Project 2 PhD/Species_images")
for (i in 1:nrow(fit$points))
{
  png <- readPNG(paste0("C:/Users/Dest/Desktop/Project 2 PhD/Species_images/", species_list[i], ".png", sep=""), native=T)
  rasterImage(png, xleft = fit$points[i,1], ybottom = fit$points[i,2], xright = (fit$points[i,1]+1), ytop = (fit$points[i,2]+1))
}
text(x, y-0.2, labels = row.names(cl), cex=.7, col="steelblue3")
mtext(side = 1, line = 2.4, cex=0.7, "NMDS1")
mtext(side = 2, line = 2.4, cex=0.7, "NMDS2")

### Partial clustering with niche volume and niche overlap ###

library("dtw")
setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Hypervolume analyses")
nam <- as.vector(NULL)
for(ss in seq_along(dir(pattern = "_volumes.txt")[-c(7,9,11,14)]))
{
  n <- strsplit(dir(pattern = "_volumes.txt")[-c(7,9,11,14)][ss], split="_volumes.txt")[[1]][1]
  nam <- c(nam, n)
}

species_list <- nam

for (tt in seq_along(species_list))
{
  Vol <- read.delim(paste(species_list[tt], "_volumes.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(Vol))
  vol2 <- t(Vol)
  Hum_ov <-read.delim(paste(species_list[tt],"_humans_overlap.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(Hum_ov))
  over2 <- t(Hum_ov)
  comb1 <- cbind(Vol, Hum_ov)
  comb_tim <- as.vector(NULL)
  for (i in 1:nrow(comb1))
  {
    tim <- as.numeric(strsplit(strsplit(rownames(comb1)[i], split="_")[[1]][3], split="k")[[1]][1])
    comb_tim <- c(comb_tim, tim)
  }
  comb2 <- data.frame(Time=comb_tim, comb1)
  timeX <- comb2$Time
  comb3 <- t(comb2[,-1])
  colnames(comb3) <- timeX
  assign(paste0(species_list[tt], "_timeseries", sep=""), comb3)
}
dist_len <- length(ls(pattern = "_timeseries"))
cl <- matrix(NA, nrow=dist_len, ncol=dist_len)
dimnames(cl) <- list(species_list, species_list)
for (i in seq_along(ls(pattern = "_timeseries")))
{
  for(j in seq_along(ls(pattern = "_timeseries")))
  {
    if (i!=j)
    {
      tim1 <- get(ls(pattern = "_timeseries")[i])
      tim2 <- get(ls(pattern = "_timeseries")[j])
      cl[i,j] <- dtw(t(tim1), t(tim2), distance.only = FALSE)$normalizedDistance  ## normalized distance matrix of all timeseries (even different lengths)
    }
  }
}

cl2 <- as.dist(cl) 
library(MASS)
library(jpeg)
library(png)
fit <- isoMDS(cl2, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
#pdf(file=paste("Final_clustering_NMDS.pdf"), width=12, height = 9)
#par(oma=c(3,3,3,3))
plot(x, y, xlab="", ylab="", main="Nonmetric MDS", cex.main=0.8, type="n", xlim=c(min(fit$points[,1])-1, max(fit$points[,1]+1)),
     ylim=c(min(fit$points[,2])-1, max(fit$points[,2]+1)), bty="l", axes=F, xaxs="i")
axis(side=1, col="black", at=seq(-7,11,2), lwd=1, col.ticks="black", cex.axis=0.7)
axis(side=2, col="black", at=seq(-7,7,2), lwd=1, col.ticks="black", cex.axis=0.7)
setwd("C:/Users/Dest/Desktop/Project 2 PhD/Species_images")
for (i in 1:nrow(fit$points))
{
  png <- readPNG(paste0("C:/Users/Dest/Desktop/Project 2 PhD/Species_images/", species_list[i], ".png", sep=""), native=T)
  rasterImage(png, xleft = fit$points[i,1], ybottom = fit$points[i,2], xright = (fit$points[i,1]+1), ytop = (fit$points[i,2]+1))
}
text(x, y-0.2, labels = row.names(cl), cex=.7, col="purple")
mtext(side = 1, line = 2.4, cex=0.7, "NMDS1")
mtext(side = 2, line = 2.4, cex=0.7, "NMDS2")
#dev.off()


### Partial clustering with geographic and niche overlap ###

library("dtw")
setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Hypervolume analyses")
nam <- as.vector(NULL)
for(ss in seq_along(dir(pattern = "_volumes.txt")[-c(7,9,11,14)]))
{
  n <- strsplit(dir(pattern = "_volumes.txt")[-c(7,9,11,14)][ss], split="_volumes.txt")[[1]][1]
  nam <- c(nam, n)
}

species_list <- nam

for (tt in seq_along(species_list))
{
 
  Hum_ov <-read.delim(paste(species_list[tt],"_humans_overlap.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(Hum_ov))
  comb_tim <- as.vector(NULL)
  for (i in 1:nrow(Hum_ov))
  {
    tim <- as.numeric(strsplit(strsplit(rownames(Hum_ov)[i], split="_")[[1]][3], split="k")[[1]][1])
    comb_tim <- c(comb_tim, tim)
  }
  comb2 <- data.frame(Time=comb_tim, Hum_ov)
  setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Geographic Overlap per time inerval/Geographic oVerlap 3 Variables")
  maxent_overlap <- read.delim(paste("Geographic_overlap_", species_list[tt], "_humans.txt", sep=""), h=T, sep="\t", stringsAsFactors = FALSE)
  print(nrow(maxent_overlap))
  colnames(maxent_overlap) <- c("Time", "Geogr_Overlap")
  comb3 <- merge(comb2, maxent_overlap, by.x="Time", by.y="Time", all.x=FALSE, all.y=TRUE, sort=TRUE)
  timeX <- comb3$Time
  comb4 <- t(comb3[,-1])
  colnames(comb4) <- timeX
  assign(paste0(species_list[tt], "_timeseries", sep=""), comb4)
  setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Hypervolume analyses")
}

dist_len <- length(ls(pattern = "_timeseries"))
cl <- matrix(NA, nrow=dist_len, ncol=dist_len)
dimnames(cl) <- list(species_list, species_list)
for (i in seq_along(ls(pattern = "_timeseries")))
{
  for(j in seq_along(ls(pattern = "_timeseries")))
  {
    if (i!=j)
    {
      tim1 <- get(ls(pattern = "_timeseries")[i])
      tim2 <- get(ls(pattern = "_timeseries")[j])
      cl[i,j] <- dtw(t(tim1), t(tim2), distance.only = FALSE)$normalizedDistance  ## normalized distance matrix of all timeseries (even different lengths)
    }
  }
}

cl2 <- as.dist(cl) 
library(MASS)
library(jpeg)
library(png)
fit <- isoMDS(cl2, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="", ylab="", main="Nonmetric MDS", cex.main=0.8, type="n", xlim=c(-0.2,0.2),
     ylim=c(-0.2,0.2), bty="l", axes=F, xaxs="i")
axis(side=1, col="black", at=seq(-0.2,0.2,0.1), lwd=1, col.ticks="black", cex.axis=0.7, labels=TRUE)
axis(side=2, col="black", at=seq(-0.2,0.2,0.1), lwd=1, col.ticks="black", cex.axis=0.7)
setwd("C:/Users/Dest/Desktop/Project 2 PhD/Species_images")
for (i in 1:nrow(fit$points))
{
  png <- readPNG(paste0("C:/Users/Dest/Desktop/Project 2 PhD/Species_images/", species_list[i], ".png", sep=""), native=T)
  rasterImage(png, xleft = fit$points[i,1], ybottom = fit$points[i,2], xright = (fit$points[i,1]+ 0.04), ytop = (fit$points[i,2]+0.04))
}
text(x-0.02, y+0.01, labels = row.names(cl), cex=.6, col="orange")
mtext(side = 1, line = 2.4, cex=0.7, "NMDS1")
mtext(side = 2, line = 2.4, cex=0.7, "NMDS2")

