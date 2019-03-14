setwd("D:/Geographic Overlap New humans data")
require(fBasics)
require(gplots)
require(pheatmap)
require(lattice)
require(latticeExtra)
bins<- c(seq(8,21,1), seq(22,46,2))
all_over <- read.delim("All_overlaps.txt", h=T, sep="\t", stringsAsFactors=FALSE, quote="")
all_thresh <- read.delim("All_overlaps_threshold.txt", h=T, sep="\t", stringsAsFactors=FALSE, quote="")

nam <- as.vector(NULL)
for(ss in seq_along(dir(pattern = "_humans_overlaps")[seq(1,21,2)]))
{
  n <- strsplit(dir(pattern = "_humans_overlaps")[seq(1,21,2)][ss], split="_humans_overlaps")[[1]][1]
  nam <- c(nam, n)
}
rownames(all_over) <- nam
colnames(all_over) <- bins
all_over2 <- data.matrix(all_over[-c(5),c(4:27)])

### Plot with levelplot ####
colos <- seqPalette(17, "Blues")
cols2 <- c(colos[2:11])

all_over3 <- t(all_over2[nrow(all_over2):1,])
rows=c(1:24)
cols=c(1:10)
gri <- expand.grid(x=rows, y=cols)
gri$z <- as.vector(all_over3)
myPanel <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...,)
    panel.text(x, y, ifelse(is.na(gri$z),"",round(gri$z,2)), cex=0.9, col="dark red")
}
par(oma=c(4,4,4,4), mar=c(3,3,3,3))

levelplot(all_over3, margin=TRUE, shrink=c(0.8,0.8), at=seq(0,1,0.1),
colorkey=list(space="bottom", col = cols2, at=seq(0,1,0.1), labels=c(as.character(seq(0,1,0.1)))), par.settings=list(layout.heights=list(xlab.key.padding=1.5)),
main = list(paste("Mammal-Humans Overlaps"), cex=0.85),  scales = list(tck = c(0.8,0)), 
xlab="Time(ka)", ylab="Species", col.regions = cols2, panel=myPanel)


#### Plot overlaps with threshold ###
all_thresh <- read.delim("All_overlaps_threshold.txt", h=T, sep="\t", stringsAsFactors=FALSE, quote="")
nam <- as.vector(NULL)
for(ss in seq_along(dir(pattern = "_humans_overlaps")[seq(2,22,2)]))
{
  n <- strsplit(dir(pattern = "_humans_overlaps")[seq(2,22,2)][ss], split="_humans_overlaps")[[1]][1]
  nam <- c(nam, n)
}
rownames(all_thresh) <- nam
colnames(all_thresh) <- bins

all_thresh2 <- data.matrix(all_thresh[-c(5),c(4:27)])

### Plot with levelplot ####
colos <- seqPalette(17, "Greens")
cols2 <- c(colos[2:11])

all_thresh3 <- t(all_thresh2[nrow(all_thresh2):1,])
rows=c(1:24)
cols=c(1:10)
gri <- expand.grid(x=rows, y=cols)
gri$z <- as.vector(all_thresh3)
myPanel <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...,)
    panel.text(x, y, ifelse(is.na(gri$z),"",round(gri$z,2)), cex=0.9, col="brown")
}
par(oma=c(4,4,4,4), mar=c(3,3,3,3))

levelplot(all_thresh3, margin=TRUE, shrink=c(0.8,0.8), at=seq(0,1,0.1),
colorkey=list(space="bottom", col = cols2, at=seq(0,1,0.1), labels=c(as.character(seq(0,1,0.1)))), par.settings=list(layout.heights=list(xlab.key.padding=1.5)),
main = list(paste("Mammal-Humans Overlaps"), cex=0.85),  scales = list(tck = c(0.8,0)), 
xlab="Time(ka)", ylab="Species", col.regions = cols2, panel=myPanel)

### PLot with heatmap ###

my_palette <- colorRampPalette(c("grey90", "orange", "dark red"))(n = 299)

col.breaks <- c(seq(0.0001, 0.1, length=30), seq(0.1001, 0.2, length=30), seq(0.2001, 0.3, length=30), seq(0.3001, 0.4, length=30), seq(0.4001, 0.5, length=30),
seq(0.5001, 0.6, length=30), seq(0.6001, 0.7, length=30), seq(0.7001, 0.8, length=30), seq(0.8001, 0.9, length=30), seq(0.9001, 1, length=30))

distance = dist(all_over2, method = "manhattan")
cluster = hclust(distance, method = "average")

hclustfunc <- function(x) hclust(x, method="average")
distfunc <- function(x) dist(x,method="manhattan")
# obtain the clusters
fit <- hclustfunc(distfunc(all_over2))
clusters <- cutree(fit, 3) 

par(oma=c(4,4,4,4), mar=c(3,3,3,3))

all_over_heat <- heatmap.2(all_over2, cellnote=ifelse(is.na(all_over2),"",round(all_over2,2)), main="Human-Mammal Overlap", notecol="black", density.info="none",
trace="none", col=my_palette, breaks= col.breaks, Rowv=as.dendrogram(cluster),
key=TRUE, keysize=0.85, key.par = list(cex=0.75), Colv=F, dendrogram="row", na.rm=TRUE, cexRow=0.7)
