rm(list = ls())
graphics.off()
gc()

require(rgdal)
library(ncdf4)
library(maptools)
library(fields)


# Define name of the region for the outputs
region=toupper("Northern_Territory")

# Define time limits
start_year <- 2000
end_year <- 2021
years = start_year:end_year

# Define coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Load world limits
data(wrld_simpl)

# Define data directories
dir_out = '~/Dropbox/model/fire_database/out_def/AUSTRALIA/'
dir_obs <- "~/Dropbox/model/fire_database/source/"
dir_fire = "~/Documents/dati/fire_climate_data/fire/australia/"

## load region
file_shp = paste0(dir_out, 'misc/northern_australia/northern_australia.shp')
shp_ll_tmp <- readOGR(file_shp)
proj4string(shp_ll_tmp) <- CRS.new
aux=shp_ll_tmp@bbox

## Load grid
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

## get data
# web_name="https://firenorth.org.au/nafi3/downloads/firescars/"
# 
# for (iyear in start_year:end_year) {
#   download.file(
#     paste0(web_name, iyear, "/", iyear, "%20firescar%20shapefiles.zip"),
#     destfile = paste0(dir_fire, "Northern_Territory/firescar", iyear, ".zip"),
#     method = "wget",
#     extra = "-r -p --random-wait"
#   )
#   unzip(paste0(dir_fire, "Northern_Territory/firescar", iyear, ".zip"), exdir=paste0(dir_fire, "Northern_Territory/firescar"))
# }

## read data and interpolate
## fire  to 1x1 grid
BA = array(data = NA, dim = c(length(lon01), length(lat01), length(years) *
                                12))

for (iyear in start_year:end_year) {
  mydata <-
    readOGR(
      paste0(dir_fire, "Northern_Territory/firescar/"),
      paste0("fs", substr(as.character(iyear), 3, 4), "_mths_gda"),
      verbose = FALSE
    )
  mesi = mydata$Month
  if (iyear==2011 ||iyear==2012 ||iyear==2014 ||iyear==2021) {mesi = mydata$month}
  centroids <- getSpPPolygonsLabptSlots(mydata)
  lon_fire = centroids[, 1]
  lat_fire = centroids[, 2]
  ba = area(mydata)
  for (i in 1:length(lon01)) {
    print(paste0('lon ',
                 i,
                 '/',
                 length(lon01)))
    for (j in 1:length(lat01)) {
      idx = which(
        lon_fire >= lon01[i] - 0.5 &
          lon_fire <= lon01[i] + 0.5 &
          lat_fire >= lat01[j] - 0.5 & lat_fire <= lat01[j] + 0.5
      )
      if (length(idx) >= 1) {
        k = (iyear - start_year) * 12
        for (imonth in 1:12) {
          k = k + 1
          idx2 = which(mesi[idx] == imonth)
          if (length(idx2) >= 1) {
            BA[i, j, k] = sum(ba[idx[idx2]], na.rm = TRUE)
          }
        }
      }
    }
  }
}


dim(BA)

image.plot(lon01, lat01, apply(log(BA), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)


lon=lon01
lat=lat01

ii <- !is.na(over(pts, shp_ll_tmp[1,1]))
inout = ii
dim(inout) <- c(length(lon01), length(lat01))
inout[inout == 0] = NA
image.plot(lon01, lat01, inout)
plot(wrld_simpl, add = TRUE)
for (k in 1:dim(BA)[3]) {
  BA[, , k] = BA[, , k]*inout
}

# export as RData
save(BA, file = paste0(dir_out, "BA_",region,"_v1.RData"))
save(lon, file = paste0(dir_out, "lon_BA_",region,"_v1.RData"))
save(lat, file = paste0(dir_out, "lat_BA_",region,"_v1.RData"))
