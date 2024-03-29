#===============================================================================
# Description: Prepares region-level fire data
#===============================================================================

#===============================================================================
# 1). Preliminary -----
#===============================================================================

# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("raster","fields","ncdf4","rgdal", "rgeos","R.matlab","cshapes","sf")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Directories
dir <- list()
dir$out = '~/Dropbox/model/fire_database/out_def/US/'
dir$fire = '~/Documents/dati/fire_us/mtbs_perimeter_data_2021/'
dir$obs <- "~/Documents/dati/fire_climate_data/fire/"

# Define time limits
start_year <- 1984
end_year <- 2020


# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

region="US_MTBS"
#===============================================================================
# 2). MTBS data -----
#===============================================================================

# The input file geodatabase
mtbs <- readOGR(paste0(dir$fire,"/mtbs_perims_DD.shp"), verbose = FALSE)
# summary(mtbs)
mtbs = spTransform(mtbs,CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
all_ba=as.numeric(mtbs$BurnBndAc) #Number of acres mapped
all_ba=all_ba*4046.86 #acres a m2
lon_fire=as.numeric(mtbs$BurnBndLon) #Longitude of the mapped centroid of fire perimeter
lat_fire=as.numeric(mtbs$BurnBndLat) #Latitude of the mapped centroid of fire perimeter
anni=as.numeric(format(as.Date(mtbs$Ig_Date),"%Y"))
mesi=as.numeric(format(as.Date(mtbs$Ig_Date),"%m"))
aux = mtbs@bbox


cc
# Load grid
fname <- file.path(dir$obs, 'land_sea_mask_1degree.nc4')
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
            BA[i, j, k] = sum(all_ba[idx[idx2]], na.rm = TRUE)
          }
        }
      }
    }
  }
}

# Load US limits
cshp <-
  cshp(as.Date("2019-1-1"),
       dependencies = FALSE)
shp_ll = as_Spatial(cshp[1,1]) #1, is U.S.
proj4string(shp_ll) <- CRS.new


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
plot(shp_ll, add = TRUE)

lon=lon01[ilon]
lat=lat01[ilat]

BA=BA[ilon,ilat,]
image.plot(lon, lat, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(shp_ll, add = TRUE)



# export as RData
save(BA, file = paste0(dir$out, "BA_",region,"_v1.RData"))
save(lon, file = paste0(dir$out, "lon_BA_",region,"_v1.RData"))
save(lat, file = paste0(dir$out, "lat_BA_",region,"_v1.RData"))

# export as netcdf
# define dimensions
londim <- ncdim_def("lon", "degrees_east", as.double(lon))
latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
fechas <- seq(as.Date(paste0(start_year,"/1/15")), as.Date(paste0(end_year,"/12/15")), by = "month")
time = seq(0, length(fechas)-1)
tunits <- paste0("months since ",start_year,"-01-15")
timedim <- ncdim_def("time", tunits, as.double(time))
# define variables
mean_def <-
  ncvar_def("Monthly_Burned_Area",
            "square meter",
            list(londim, latdim, timedim),1e32,prec = "double")
# create netCDF file and put arrays
ncfname <- paste0(dir$out, "/Monthly_Burned_area_",start_year,"_",end_year,"_",region,"_v1.nc")
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
  'An integrated dataset of observed national inventories of fires - US MTBS data'
)
ncatt_put(ncout, 0, "Version", '1.0')
ncatt_put(ncout, 0, "Institution", "University of Murcia (Spain)")
ncatt_put(ncout, 0, "Url", "")
ncatt_put(ncout, 0, "Creator", "Marco Turco <marco.turco@um.es")
ncatt_put(ncout,
          0,
          "Software",
          "Create in R \n\t https://github.com/marcoturco/")
ncatt_put(ncout, 0, "Date", date())
# ncatt_put(
#   ncout,
#   0,
#   "reference",
#   "An integrated dataset of observed national inventories of fires"
# )

# close the file, writing data to disk
nc_close(ncout)

## Export as geoTIFF
BA <- brick(ncfname,  varname = "Monthly_Burned_Area")
writeRaster(BA,paste0(substr(ncfname,1,nchar(ncfname)-3),'.tif'),options=c('TFW=YES'),overwrite=TRUE)

# Save the data frame as a csv file
df <- as.data.frame(BA, xy=T)
# Define the start and end years
# Generate a sequence of years and months
years <- seq(from = start_year, to = end_year, by = 1)
months <- rep(seq(from = 1, to = 12, by = 1), length(years))
# Create the column names by combining the year and month
col_names <- paste(years, formatC(months, width = 2, flag = "0"), sep = "-")
# Set the column names of the data frame
colnames(df)[3:dim(df)[2]] <- col_names
colnames(df)[1] <- "lon"
colnames(df)[2] <- "lat"
write.csv(df, file = paste0(substr(ncfname,1,nchar(ncfname)-3),'.csv'), row.names = F)