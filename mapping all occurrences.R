library(maptools)
library(mapproj)
library(rgdal)
library(maps)
library(mapdata)
library(fBasics)


setwd("C:/Users/Dest/Desktop/Project 2 PhD/2nd Chpapter ANALYSES/Human localities again/Geographic Overlap per time inerval/Geographic oVerlap 3 Variables/Geographic Overlap 3 variables resampling")

species_DB <- read.delim("Species_environmental_variables_new_way_different_bins.txt", h=T, sep="\t")
Human_DB <- read.delim("Humans_environmental_variables_new_way_different_bins.txt", h=T, sep="\t")
#species_list <- as.vector(unique(species_DB$Species)) ## get species names - in our case can be time bins as well

spec_Intervals <- sort(as.vector(unique(species_DB$Interval)), decreasing=FALSE)
spec_cols <- seqPalette(length(spec_Intervals) + 7, "Oranges")
spec_cols2 <- spec_cols[c(8:(length(spec_Intervals) + 7))]  
temp_spec <- cbind.data.frame(Interval=spec_Intervals, Interval_col=as.character(spec_cols2))
species_map <- cbind.data.frame(Species=species_DB$Species, lon=species_DB$Longitude, lat=species_DB$Latitude, Interval=species_DB$Interval)
for(i in seq_along(species_map$Interval))
{
  for(co in seq_along(temp_spec$Interval))
  {
    if (species_map$Interval[i] == temp_spec$Interval[co])
    {
      species_map$color[i] <- as.character(temp_spec$Interval_col[co])
    }
  }
}
spec_coords <- mapproject(species_map$lon, species_map$lat, proj="azequidistant", orientation=c(90,0,30))

hum_Intervals <- sort(as.vector(unique(Human_DB$Interval)), decreasing=FALSE)
hum_cols <- seqPalette(length(hum_Intervals) + 7, "Greens")
hum_cols2 <- hum_cols[c(8:(length(hum_Intervals) + 7))]  
temp_hum <- cbind.data.frame(Interval=hum_Intervals, Interval_col=as.character(hum_cols2))
human_map <- cbind.data.frame(Species=Human_DB$Species, lon=Human_DB$Longitude, lat=Human_DB$Latitude, Interval=Human_DB$Interval)
for(i in seq_along(human_map$Interval))
{
  for(co in seq_along(temp_hum$Interval))
  {
    if (human_map$Interval[i] == temp_hum$Interval[co])
    {
      human_map$color[i] <- as.character(temp_hum$Interval_col[co])
    }
  }
}
hum_coords <- mapproject(human_map$lon, human_map$lat, proj="azequidistant", orientation=c(90,0,30))

par(mfrow=c(1,2))
sp_map <- map("worldHires", projection="azequidistant", orientation=c(90,0,30),xlim=c(-180,180), ylim=c(30,90), col="grey25", fill=TRUE)
#box()
points(spec_coords, pch=21, cex=1.2, col="steelblue", bg=species_map$color)
title(main="Mammal assemblages occurrences", cex.main=0.7)
#legend("topleft", legend=temp_spec$Interval, col="steelblue", pt.bg=as.character(temp_spec$Interval_col), pch=21)
par(new=T)
map("worldHires", projection="azequidistant", orientation=c(90,0,30),xlim=c(-180,180), ylim=c(30,90), col="grey25", fill=TRUE)
#box()
points(hum_coords, pch=21, cex=1.2, col="steelblue", bg=human_map$color)
title(main="Modern Human occurrences", cex.main=0.7)
