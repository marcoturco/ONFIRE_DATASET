# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("rgdal","R.matlab","ncdf4","maptools","fields","cshapes","sf","sp","raster")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Define name of the region for the outputs
region="CANADA_NBAC"

# Define time limits
start_year <- 1986
end_year <- 2020

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Load world limits
data(wrld_simpl)

# Load Canada limits
cshp <-
  cshp(as.Date("2019-1-1"),
       dependencies = FALSE)
shp_ll = as_Spatial(cshp[2,1]) #2, is Canada
proj4string(shp_ll) <- CRS.new
aux = shp_ll@bbox

# Define data directories
dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "/diskonfire/ONFIREDATA/canada_nbac/"
dir_grid <- "~/Dropbox/model/fire_database/source/"

# Load fire data from NFDB
mydata<- readOGR(paste0(dir_obs, "nbac_1986_to_2020_20210810.shp"), verbose = FALSE)
mydata <-
  spTransform(mydata,
              CRS(
                "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
              ))
centroids <- getSpPPolygonsLabptSlots(mydata)
lon_fire=centroids[,1]
lat_fire=centroids[,2]

# Get fire data attributes
ba = mydata$ADJ_HA
# plot(mydata$POLY_HA,mydata$ADJ_HA)
anni = format(as.Date(mydata$AFSDATE),"%Y")
anni = as.numeric(as.character(anni))
mesi = format(as.Date(mydata$AFSDATE),"%m")
mesi=as.numeric(as.character(mesi))

# Load grid
fname <- file.path(dir_grid, 'land_sea_mask_1degree.nc4')
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

# Calculate BA for each grid cell
BA = array(data = NA, dim = c(length(lon01), length(lat01), length(start_year:end_year) *
                                12))
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
      k = 0
      for (iyear in start_year:end_year) {
        for (imonth in 1:12) {
          k = k + 1
          idx2 = which(anni[idx] == iyear & mesi[idx] == imonth)
          if (length(idx2) >= 1) {
            BA[i, j, k] = sum(ba[idx[idx2]], na.rm = TRUE)
          }
        }
      }
    }
  }
}

# chech the dataset
summary(as.vector(BA))
image.plot(lon01, lat01, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(shp_ll, add = TRUE)

# Mask points outside the U.S.
ii <- !is.na(over(pts, shp_ll))
inout = ii
dim(inout) <- c(length(lon01), length(lat01))
inout[inout == 0] = NA
# image.plot(lon01, lat01, inout)
# plot(shp_ll, add = TRUE)
for (it in 1:dim(BA)[3]) {
  BA[,,it]=BA[,,it]*inout
}

lonmin=max(lon01)
lonmax=min(lon01)
latmin=max(lat01)
latmax=min(lat01)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    if (length(which(!is.na(BA[i,j,]))) >= 1) {
      idx=which(is.na(BA[i,j,]))
      BA[i,j,idx]=0
      lonmin=min(lonmin,lon01[i])
      lonmax=max(lonmax,lon01[i])
      latmin=min(latmin,lat01[j])
      latmax=max(latmax,lat01[j])
    }
  }
}
ilon=which(lon01>=lonmin & lon01<=lonmax)
ilat=which(lat01>=latmin & lat01<=latmax)
image.plot(lon01[ilon], lat01[ilat], apply(log(BA[ilon,ilat,]+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

lon=lon01[ilon]
lat=lat01[ilat]

BA=BA[ilon,ilat,]
image.plot(lon, lat, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(shp_ll, add = TRUE)

#from hectare to square meter
BA=BA*10000

# export as RData
save(BA, file = paste0(dir_out, "BA_",region,"_v1.RData"))
save(lon, file = paste0(dir_out, "lon_BA_",region,"_v1.RData"))
save(lat, file = paste0(dir_out, "lat_BA_",region,"_v1.RData"))