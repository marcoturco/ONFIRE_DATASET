rm(list = ls())
graphics.off()
gc()

library(ncdf4)
library(fields)
library(maptools)
library(rgdal)        # para leer shapefiles
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data(wrld_simpl)
proj4string(wrld_simpl) <- CRS.new
years=1950:2021

# Define data directories
dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Dropbox/model/fire_database/source/"

# Load grid
fname <- file.path(dir_obs, 'land_sea_mask_1degree.nc4')
obs.nc <- nc_open(fname)
lon = obs.nc$dim$lon$vals 
lat = obs.nc$dim$lat$vals
ilon=which(lon>=110 & lon<=155)
ilat=which(lat>=-45 & lat<=-10)
lon_aus=lon[ilon]
lat_aus=lat[ilat]
rm(lon)
rm(lat)
points <- expand.grid(lon_aus, lat_aus)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new

BA_ALL = array(data = 0, dim = c(length(lon_aus),length(lat_aus),length(years)*12))
# load grid 
region="TASMANIA"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
ilon=match(lon,lon_aus)
ilat=match(lat,lat_aus)
BA_ALL[ilon,ilat,]=BA #[,,1:(dim(BA)[3]-12)]
BA_ALL[ilon,ilat,1:((1960-1950)*12)]=NA
image.plot(lon_aus, lat_aus, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

region="VICTORIA"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_AUSTRALIA_",region,"_v1.RData"))
ilon=match(lon01,lon)
ilat=match(lat01,lat)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    for (k in 1:dim(BA_ALL)[3]) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL[ilon[i], ilat[j], k] = sum(BA_ALL[ilon[i], ilat[j], k], BA[i, j, k], na.rm =
                                            T)
      }
    }
  }
}
image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)


region="AUSTRALIA_NWS"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
ilon=match(lon01,lon)
ilat=match(lat01,lat)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    for (k in 1:dim(BA_ALL)[3]) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL[ilon[i], ilat[j], k] = sum(BA_ALL[ilon[i], ilat[j], k], BA[i, j, k], na.rm =
                                            T)
      }
    }
  }
}
image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

region="SOUTH_AUSTRALIA"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
ilon=match(lon01,lon)
ilat=match(lat01,lat)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    for (k in 1:dim(BA_ALL)[3]) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL[ilon[i], ilat[j], k] = sum(BA_ALL[ilon[i], ilat[j], k], BA[i, j, k], na.rm =
                                            T)
      }
    }
  }
}
image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

region="WESTERN_AUSTRALIA"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
ilon=match(lon01,lon)
ilat=match(lat01,lat)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    for (k in 1:dim(BA_ALL)[3]) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL[ilon[i], ilat[j], k] = sum(BA_ALL[ilon[i], ilat[j], k], BA[i, j, k], na.rm =
                                            T)
      }
    }
  }
}
image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

region="QUEENSLAND"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
image.plot(lon01, lat01, apply((BA), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)
ilon=match(lon01,lon)
ilat=match(lat01,lat)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    for (k in 1:dim(BA_ALL)[3]) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL[ilon[i], ilat[j], k] = sum(BA_ALL[ilon[i], ilat[j], k], BA[i, j, k], na.rm =
                                            T)
      }
    }
  }
}
image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

image.plot(lon, lat, BA_ALL[,,133])
plot(wrld_simpl, add = TRUE)

file_shp = '/Users/marco/Documents/virtualBox/temiav/fire_database/STE_2021_AUST_GDA2020.shp'
shp_ll_tmp <- readOGR(file_shp)
proj4string(shp_ll_tmp) <- CRS.new
# plot(shp_ll_tmp[7, 1])
ii <- is.na(over(pts, shp_ll_tmp[7, 1]))
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image.plot(lon, lat, inout)
plot(wrld_simpl, add = TRUE)
for (k in 1:dim(BA_ALL)[3]) {
    BA_ALL[, , k] = BA_ALL[, , k]*inout
}
# shp_ll_tmp$STE_NAME21


ii <- !is.na(over(pts, wrld_simpl))
inout = ii[,1]
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image.plot(lon, lat, inout)
plot(wrld_simpl, add = TRUE)
for (k in 1:dim(BA_ALL)[3]) {
  BA_ALL[, , k] = BA_ALL[, , k]*inout
}

image.plot(lon, lat, apply((BA_ALL), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

#from hectare to square meter
BA=BA_ALL*10000

## write ntcdf

# define dimensions
londim <- ncdim_def("lon", "degrees_east", as.double(lon))
latdim <- ncdim_def("lat", "degrees_north", as.double(lat))


fechas <- seq(as.Date("1950/1/15"), as.Date("2021/12/15"), by = "month")
anni = 1950:2020
time = seq(0, length(fechas) - 1)
tunits <- "months since 1950-01-15 00:00:00.0 -0:00"
timedim <- ncdim_def("time", tunits, as.double(time))


# define variables
fillvalue <- 1e32
dlname <- "mean"
mean_def <-
  ncvar_def("Monthly_Burned_Area",
            "square meter",
            list(londim, latdim, timedim),
            fillvalue,
            dlname,
            prec = "double")

# create netCDF file and put arrays
ncfname <- paste0(dir_out, "/Monthly_Burned_area_1950_2021_Australia_v1.nc")
ncout <-
  nc_create(ncfname,
            list(mean_def))

# put variables

ncvar_put(ncout, mean_def, BA)

# put additional attributes into dimension and data variables
ncatt_put(ncout, "lon", "axis", "X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout, "lat", "axis", "Y")
ncatt_put(ncout, "time", "axis", "T")

# add global attributes
ncatt_put(
  ncout,
  0,
  "Title",
  'An integrated dataset of observed national inventories of fires - AUSTRALIA'
)
ncatt_put(ncout, 0, "Version", '1.0')
ncatt_put(ncout, 0, "Institution", "University of Murcia (Spain)")
ncatt_put(ncout, 0, "Url", "")
ncatt_put(ncout, 0, "Creator", "Marco Turco <marco.turco@um.es>")
ncatt_put(ncout,
          0,
          "Software",
          "Create in R \n\t https://github.com/marcoturco/")
ncatt_put(ncout, 0, "Date", date())
# ncatt_put(
#   ncout,
#   0,
#   "reference",
#   "TODO"
# )


# close the file, writing data to disk
nc_close(ncout)
