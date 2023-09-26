
rm(list = ls())
graphics.off()
gc()


library(ncdf4)
library(fields)
library(maptools)

data(wrld_simpl)

dir_fire = '/Users/marco/Documents/dati/FireCCI51'
dir_data = '/Users/marco/Dropbox/model/datos/'
dir_obs <- "~/Documents/dati/fire_climate_data/fire/"

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


fname <-
  file.path(dir_fire, '/2001/20010101-ESACCI-L4_FIRE-BA-MODIS-fv5.1.nc')
obs.nc <- nc_open(fname)
ba025 <- ncvar_get(obs.nc, "burned_area")
obs.nc$dim$lon$vals -> lon025
obs.nc$dim$lat$vals -> lat025
lat025 <- rev(lat025)
# ba025=ba025[,ncol(ba025):1]
# image.plot(lon025,lat025,log(ba025+1))
# plot(wrld_simpl,add=TRUE)


# Load grid
fname <- file.path(dir_obs, 'land_sea_mask_1degree.nc4')
obs.nc <- nc_open(fname)
lon = obs.nc$dim$lon$vals 
lat = obs.nc$dim$lat$vals
# Create grid
points <- expand.grid(lon, lat)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new

 
BA = array(data = NA, dim = c(length(lon), length(lat), length(2001:2020) *
                                12))

k = 0
for (iyear in 2001:2020) {
  print(iyear)
  for (imonth in 1:12) {
    print(imonth)
    k = k + 1
    
    fname <-
      paste0(
        dir_fire,
        '/',
        iyear,
        '/',
        iyear,
        sprintf("%02d", as.numeric(imonth)),
        '01-ESACCI-L4_FIRE-BA-MODIS-fv5.1.nc'
      )
    obs.nc <- nc_open(fname)
    # ba025 <- ncvar_get(obs.nc, "burned_area")
    aux <- ncvar_get(obs.nc, "burned_area_in_vegetation_class")
    aux=aux[,,4:dim(aux)[3]] #delete non-nat fires
    ba025=apply(aux,c(1,2),sum,na.rm=TRUE)
    ba025=ba025[,ncol(ba025):1]
    # image.plot(lon025, lat025, apply(log(ba025+1),c(1,2),mean,na.rm=TRUE))
    # plot(wrld_simpl, add = TRUE)
    for (i in 1:length(lon)) {
      for (j in 1:length(lat)) {
        idx1 = which(lon025 >= lon[i] - 0.5 &
                       lon025 <= lon[i] + 0.5)
        idx2 = which(lat025 >= lat[j] - 0.5 &
                       lat025 <= lat[j] + 0.5)
        
        if (length(idx1) & length(idx2) >= 1) {
          BA[i, j, k] = sum(ba025[idx1, idx2])
        }
      }
    }
  }
}



image.plot(lon, lat, apply(log(BA+1),c(1,2),mean,na.rm=TRUE))
plot(wrld_simpl, add = TRUE)


save(BA, file = paste0(dir_obs,"/BA_200101_202012_1degree_nat.RData"))

