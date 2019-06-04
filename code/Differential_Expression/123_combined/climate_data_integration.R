### climate_data_integration.R
### IGR, 2019.06.04
### Baby steps into pulling the climate data from our five main villages with the help of a couple of tutorials from the internet, namely:
### https://github.com/kapitzas/WorldClimTiles
### https://kapitzas.github.io/post/worldclimtiles/
### https://emilypiche.github.io/BIO381/raster.html

### I have no idea what I'm doing! 

### Last edit: IGR 2019.06.04
### Initial commit, playing with code, testing things. 





##########################################
### 0. Load dependencies and set paths ###
##########################################

library(raster)
library(WorldClimTiles)

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/spatial_testing/"
edaoutput <- paste0(outputdir, "/eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

setwd(outputdir)

### This pulls data from all of Indonesia in one go, which is neat - I can now make rainfall maps, although I'm not sure how useful they would be. I would like to pull specific regions. 

boundary <- getData("GADM", country = "IDN", level =0) #Indonesia borders
wc_tiles <- tile_name(boundary, "worldclim") #for 0.5 arcmin worldclim tiles of Indonesia
srtm_tiles <- tile_name(boundary, "srtm") #for 90m srtm tiles of Indonesia

# Download tiles
tiles <- tile_get(wc_tiles, var = "prec", name = "worldclim", path = outputdir) #download and load worldclim
tiles <- tile_get(srtm_tiles, var = "prec", name= "srtm", path = outputdir) #download and load srtm

# Merge tiles
final <- tile_merge(tiles)
final <- tile_merge(tiles, fac = 10) #Reprojects data to 10 times smaller res, i.e. to get srtm data from 90m to approx 1km resolution.

#If I plot final, do I get a lovely plot?
pdf("tester.pdf")
    plot(final)
dev.off()

### The example in the post from emilypiche.github.io is perfect:
indoVillages <- data.frame(site=c("Madobag", "Taileleu", "Anakalung", "Wunga", "Korowai"), long=c(99.0837, 99.1371, 119.575, 119.958, 139.673002), lat=c(-1.594, -1.7878, -9.588, -9.385, -5.480278))

# downloading and processing the data from worlclim, a function:
getWC <- function(dataset){
    #pull data
    r2EastIndonesia <- getData("worldclim", var=dataset, res=0.5, lat=-1, lon=120) 
    r2WestIndonesia <- getData("worldclim", var=dataset, res=0.5, lat=-1, lon=99) 

    #configure points
    pointsEast <- SpatialPoints(indoVillages[,2:3], proj4string = r2EastIndonesia@crs)
    pointsWest <- SpatialPoints(indoVillages[,2:3], proj4string = r2WestIndonesia@crs)

    # extracting the values for our city coordinates
    wcDataEast <- extract(r2EastIndonesiawcData, pointsEast) #only Korowai is in this tile
    wcDataWest <- extract(r2WestIndonesiawcData, pointsWest) 

    wcDataNA <- which(is.na(wcDataWest))
    wcDataWest[wcDataNA] <- wcDataEast[wcDataNA]

    # creating a dataframe with our results
    wcDataAll <- cbind.data.frame(indoVillages, wcDataWest)
    return(wcDataAll)
}

tMax <- getWC("tmax")
tMin <- getWC("tmin")
Prec <- getWC("prec")

# # Tmax
# r2EastIndonesiaTmax <- getData("worldclim", var="tmax", res=0.5, lat=-1, lon=120) 
# r2WestIndonesiaTmax <- getData("worldclim", var="tmax", res=0.5, lat=-1, lon=99) 

# # Tmin
# r2EastIndonesiaTmin <- getData("worldclim", var="tmin", res=0.5, lat=-1, lon=120) 
# r2WestIndonesiaTmin <- getData("worldclim", var="tmin", res=0.5, lat=-1, lon=99) 

# # Rainfall
# r2EastIndonesiaRain <- getData("worldclim", var="prec", res=0.5, lat=-1, lon=120) 
# r2WestIndonesiaRain <- getData("worldclim", var="prec", res=0.5, lat=-1, lon=99) 

# # putting our coordinates into the format r::raster likes
# pointsEast <- SpatialPoints(indoVillages[,2:3], proj4string = r2EastIndonesia@crs)
# pointsWest <- SpatialPoints(indoVillages[,2:3], proj4string = r2WestIndonesia@crs)

# # extracting the values for our city coordinates
# TmaxEast <- extract(r2EastIndonesiaTmax, pointsEast) #only Korowai is in this tile
# TmaxWest <- extract(r2WestIndonesiaTmax, pointsWest) 

# TminEast <- extract(r2EastIndonesiaTmin, pointsEast) #only Korowai is in this tile
# TminWest <- extract(r2WestIndonesiaTmin, pointsWest) 

# RainEast <- extract(r2EastIndonesiaRain, pointsEast) #only Korowai is in this tile
# RainWest <- extract(r2WestIndonesiaRain, pointsWest) 

# # Combina data the two dataframes, since we needed different tiles above. Could bruteforce it, but let's be elegant here
# RainNA <- which(is.na(RainWest))
# RainWest[RainNA] <- RainEast[RainNA]

# # creating a dataframe with our results
# RainAll <- cbind.data.frame(indoVillages, RainWest)
# print(Rain2)
