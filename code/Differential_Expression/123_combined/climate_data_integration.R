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
library(ggplot2)
library(reshape)
library(wesanderson)
library(gridExtra)
library(ggpubr)

# Set output directory and create it if it does not exist:
outputdir <- "/data/cephfs/punim0586/igallego/indoRNA/spatial_testing/"
edaoutput <- paste0(outputdir, "/eda/")

if (file.exists(outputdir) == FALSE){
    dir.create(outputdir, recursive=T)
    dir.create(edaoutput, recursive=T)
}

#Colour schemes:
korowai <- wes_palette("Zissou1", 20, type = "continuous")[20]
mentawai <- wes_palette("Zissou1", 20, type = "continuous")[1]
sumba <- wes_palette("Zissou1", 20, type = "continuous")[11]

smb_mtw <- wes_palette("Darjeeling1", 9, type = "continuous")[3]
smb_kor <- wes_palette("Darjeeling1", 9, type = "continuous")[7]
mtw_kor <- "darkorchid4"

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

# #If I plot final, do I get a lovely plot?
# pdf("tester.pdf")
#     plot(final)
# dev.off()

### The example in the post from emilypiche.github.io is perfect:
indoVillages <- data.frame(site=c("Madobag", "Taileleu", "Anakalung", "Wunga", "Korowai"), long=c(99.0837, 99.1371, 119.575, 119.958, 139.673002), lat=c(-1.594, -1.7878, -9.588, -9.385, -5.480278))

# downloading and processing the data from worlclim, a function: 
# only needs two tiles at 0.5 minute resolution, but that's ok. 
getWC <- function(dataset, resolution){
    #pull data
    r2EastIndonesia <- getData("worldclim", var=dataset, res=resolution, lat=-1, lon=120) 
    r2WestIndonesia <- getData("worldclim", var=dataset, res=resolution, lat=-1, lon=99) 

    #configure points
    pointsEast <- SpatialPoints(indoVillages[,2:3], proj4string = r2EastIndonesia@crs)
    pointsWest <- SpatialPoints(indoVillages[,2:3], proj4string = r2WestIndonesia@crs)

    # extracting the values for our city coordinates
    wcDataEast <- extract(r2EastIndonesia, pointsEast) #only Korowai is in this tile
    wcDataWest <- extract(r2WestIndonesia, pointsWest) 

    wcDataNA <- which(is.na(wcDataWest))
    wcDataWest[wcDataNA] <- wcDataEast[wcDataNA]

    # creating a dataframe with our results
    wcDataAll <- cbind.data.frame(indoVillages, wcDataWest)
    return(wcDataAll)
}

tMax <- getWC("tmax", 2.5) # in 10xdegree, need to divide
tMin <- getWC("tmin", 2.5) # in 10xdegree, need to divide
precip <- getWC("prec", 2.5) # in mm
alt <- getWC("alt", 0.5) # more specific because I don't want to average out the altitude gains

tMax[,4:15] <- tMax[,4:15]/10
tMin[,4:15] <- tMin[,4:15]/10

tMaxPlot <- melt(tMax[,c(1,4:15)])
tMinPlot <- melt(tMin[,c(1,4:15)])
precipPlot <- melt(precip[,c(1,4:15)])
# altPlot <- melt(alt[,c(1,4:15)])

pdf(file="climate_vars.pdf")
    tMaxgg <- ggplot(tMaxPlot, aes(x=variable, y=value, group=site, fill=site, color=site)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        labs(title="", y="mean monthly max temperature (1970-2015), 2.5 arc min res", x="")

    tMingg <- ggplot(tMinPlot, aes(x=variable, y=value, group=site, fill=site, color=site)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        labs(title="", y="mean monthly min temperature (1970-2015), 2.5 arc min res", x="")

    precipgg <- ggplot(precipPlot, aes(x=variable, y=value, group=site, fill=site, color=site)) +
        geom_line() +
        geom_point() +
        theme_bw() +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        labs(title="", y="mean monthly precipitation (mm; 1970-2015), 2.5 arc min res", x="")

    print(tMaxgg)
    print(tMingg)
    print(precipgg)
dev.off()

tMax05 <- getWC("tmax", 0.5) # in 10xdegree, need to divide
tMin05 <- getWC("tmin", 0.5) # in 10xdegree, need to divide
precip05 <- getWC("prec", 0.5)

tMax05[,4:15] <- tMax05[,4:15]/10
tMin05[,4:15] <- tMin05[,4:15]/10

tMax05Plot <- melt(tMax05[,c(1,4:15)])
tMin05Plot <- melt(tMin05[,c(1,4:15)])
precip05Plot <- melt(precip05[,c(1,4:15)])

villageCols <- c("Korowai" = korowai, "Taileleu" = mentawai, "Madobag" = "steelblue4", "Wunga" = sumba, "Anakalung" = "goldenrod")

pdf(file="climate_vars_0.5.pdf", width=3.5, height=10)
    tMaxgg <- ggplot(tMax05Plot, aes(x=variable, y=value, group=site, color=site)) +
        geom_line() +
        geom_point(size = 2) +
        theme_bw() +
        theme(legend.position="top") +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        scale_color_manual(name = "Sampling Site",values = villageCols) +
        labs(title="", y="Mean monthly max temperature (C; 1970-2000)", x="") +
        guides(color=guide_legend(nrow=2))

    tMingg <- ggplot(tMin05Plot, aes(x=variable, y=value, group=site, color=site)) +
        geom_line() +
        geom_point(size = 2) +
        theme_bw() +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        scale_color_manual(name = "Sampling Site",values = villageCols) +
        labs(title="", y="Mean monthly min temperature (C; 1970-2000)", x="") +
        guides(color=F)

    precipgg <- ggplot(precip05Plot, aes(x=variable, y=value, group=site, color=site)) +
        geom_line() +
        geom_point(size = 2) +
        theme_bw() +
        scale_x_discrete(labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
        scale_color_manual(name = "Sampling Site",values = villageCols) +
        labs(title="", y="Mean monthly precipitation (mm; 1970-2000)", x="") +
        guides(color=F)

    ggarrange(tMaxgg, tMingg, precipgg, common.legend=T, legend="top", ncol=1, nrow=3)

    # print(tMaxgg)
    # print(tMingg)
    # print(precipgg)
dev.off()
    