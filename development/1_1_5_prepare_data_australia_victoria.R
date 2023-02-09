rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(ncdf4)
library(maptools)
library(fields)


# Define name of the region for the outputs
region="VICTORIA"

# Define time limits
start_year <- 1950
end_year <- 2021
years = start_year:end_year

# Define coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Load world limits
data(wrld_simpl)

# Define data directories
dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Dropbox/model/fire_database/source/"
dir_fire = '/Users/marco/Documents/dati/fire_climate_data/fire/australia/victoria/'

## load data
mydata<- readOGR(dir_fire, "FIRE_HISTORY", verbose = FALSE)
mydata <-
  spTransform(mydata,
              CRS(
                "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
              ))

ba = mydata$AREA_HA
centroids <- getSpPPolygonsLabptSlots(mydata)
lon_fire=centroids[,1]
lat_fire=centroids[,2]
anni = as.numeric(substr(mydata$STRTDATIT,1,4))
mesi = as.numeric(substr(mydata$STRTDATIT,5,6))

aux = mydata@bbox

summary(as.vector(anni))
summary(as.vector(mesi))

# Load grid
fname <- file.path(dir_obs, 'land_sea_mask_1degree.nc4')
obs.nc <- nc_open(fname)
lon = obs.nc$dim$lon$vals 
lat = obs.nc$dim$lat$vals
lon01 = lon
lat01 = lat
ilon = which(lon01 >= aux[1, 1] & lon01 <= aux[1, 2])
ilat = which(lat01 >= aux[2, 1] & lat01 <= aux[2, 2])
lon01 = lon01[ilon]
lat01 = lat01[ilat]
# Create grid
points <- expand.grid(lon01, lat01)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new

## fire  to 1x1 grid
BA = array(data = 0, dim = c(length(lon01), length(lat01), length(years) *
                                12))

for (i in 1:length(lon01)) {
  
  print(paste0(
    'lon ',
    i,
    '/',
    length(lon01)
  ))
  
  
  for (j in 1:length(lat01)) {
    
    
    
    
    idx = which(
      lon_fire >= lon01[i] - 0.5 &
        lon_fire <= lon01[i] + 0.5 &
        lat_fire >= lat01[j] - 0.5 & lat_fire <= lat01[j] + 0.5
    )
    if (length(idx) >= 1) {
      k = 0
      for (iyear in years) {
        for (imonth in 1:12) {
          k = k + 1
          idx2 = which(anni[idx] == iyear & mesi[idx] == imonth)
          if (length(idx2) >= 1) {
            
            BA[i, j, k] = sum(ba[idx[idx2]],na.rm = TRUE)
          }
          
        }
      }
      
    }
  }
}
dim(BA)

image.plot(lon01, lat01, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)
image.plot(lon01, lat01, apply((BA), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

lon=lon01
lat=lat01

# export as RData
save(BA, file = paste0(dir_out, "BA_",region,"_v1.RData"))
save(lon, file = paste0(dir_out, "lon_BA_",region,"_v1.RData"))
save(lat, file = paste0(dir_out, "lat_BA_",region,"_v1.RData"))