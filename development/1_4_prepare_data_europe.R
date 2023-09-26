# Clean up
rm(list = ls())
graphics.off()
gc()

# Packages
wants <- c("rgdal","R.matlab","ncdf4","maptools","fields","cshapes","sf","sp","raster","rio","stringr")
needs <- wants[!(wants %in% installed.packages()[, "Package"])]
if (length(needs))
  install.packages(needs)
lapply(wants, function(i)
  require(i, character.only = TRUE))
rm(needs, wants)

# Define name of the region for the outputs
region="EUROPE"

# Define time limits
start_year <- 1980
end_year <- 2020
years = start_year:end_year

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Load world limits
data(wrld_simpl)
proj4string(wrld_simpl) <- CRS.new

# Define data directories
dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Dropbox/model/fire_database/source/"

# Load fire data
my_data <-
  import_list(
    paste0(dir_obs,'Request230222_turco1.xlsx'
    ),
    which = c(1, 2)
  )
EFFIS = my_data$Data

sum(EFFIS$SumOfBATOT)
sum(EFFIS$CountOfFIREID_EU)

dim(my_data$Data)
nuts3_effis=EFFIS$NUTS3_2016
nuts3_effis=unique(nuts3_effis)

# load nuts3 regions
# Fichero .shp
file_shp = paste0(dir_obs,'nuts3_europe/NUTS_RG_01M_2016_3035.shp')
shp <- readOGR(file_shp)
shp_ll <-
  spTransform(shp,
              CRS(
                "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
              ))
# Filter NUTS polygons based on EFFIS NUTS3 codes
nuts3_filtered <- shp_ll[shp_ll$NUTS_ID %in% nuts3_effis, ]
aux = nuts3_filtered@bbox

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
# points <- expand.grid(lon01, lat01)
# pts = SpatialPoints(points)
# proj4string(pts) <- CRS.new
# Create a 0.1º grid extent based on the original 1º grid
x_min <- min(lon01)
x_max <- max(lon01)
y_min <- min(lat01)
y_max <- max(lat01)
# Define the desired resolution for the new 0.1º grid
resolution <- 0.1
# Create the 0.1º grid using the desired resolution
lon001=seq(x_min, x_max, resolution)
lat001=seq(y_min, y_max, resolution)
points <- expand.grid(lon001, lat001)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new

# Calculate BA for each grid cell
BA = array(data = NA, dim = c(length(lon001), length(lat001), length(years) *
                                12))

for (ireg in 1:length((EFFIS$NUTS3_2016))) {
  print(paste0('ireg ',
               ireg,
               '/',
               length((EFFIS$NUTS3_2016))))
  inut = which(str_detect(
    as.character(shp_ll$NUTS_ID),
    as.character(EFFIS$NUTS3_2016[ireg])
  ))
  if (length(inut)==0) {next}
  ii <- !is.na(over(pts, shp_ll[inut, ]))
  inout = ii[, 1]
  dim(inout) <- c(length(lon001), length(lat001))
  inout[inout == 0] = NA
  if (sum(inout, na.rm = TRUE) > 0) {
    # image.plot(lon001, lat001, inout)
    # plot(shp_ll[inut, ], add = TRUE)
    # plot(wrld_simpl, add = TRUE)
    for (i in 1:length(lon001)) {
      for (j in 1:length(lat001)) {
        if (!is.na(inout[i, j])) {
          BA[i, j , (match(EFFIS$YEAR[ireg], years) - 1) * 12 + match(EFFIS$MONTH[ireg], seq(1, 12))] =
            (EFFIS$SumOfBATOT[ireg] / sum(inout, na.rm = TRUE))
        }
      }
    }
  }
}





# chech the dataset
summary(as.vector(BA))
image.plot(lon001, lat001, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)


# Calculate BA for each 1º grid cell
BA2 = array(data = NA, dim = c(length(lon01), length(lat01), length(years) *
                                 12))
for (k in 1:dim(BA)[3]) {
  print(paste0('month ',k,' of ',dim(BA)[3]))
  for (i in 1:length(lon01)) {
    
    for (j in 1:length(lat01)) {
      idx1 = which(lon001 >= lon01[i] - 0.5 &
                     lon001 <= lon01[i] + 0.5)
      idx2 = which(lat001 >= lat01[j] - 0.5 &
                     lat001 <= lat01[j] + 0.5)
      
      if (length(idx1) >= 1 & length(idx2) >= 1 & length(which(is.na(BA[idx1, idx2,k]))) < length(idx1)*length(idx2) ) {
        BA2[i, j, k] = sum(BA[idx1, idx2,k],na.rm=T)
      }
    }
  }
}


summary(as.vector(BA2))
image.plot(lon01, lat01, apply(log(BA2+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)


rm(BA)
BA=BA2
rm(BA2)

lonmin=max(lon01)
lonmax=min(lon01)
latmin=max(lat01)
latmax=min(lat01)
for (i in 1:length(lon01)) {
  for (j in 1:length(lat01)) {
    if (length(which(!is.na(BA[i,j,]))) >= 1) {
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
plot(wrld_simpl, add = TRUE)

## Find which points fall over land
points <- expand.grid(lon, lat)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new
ii <- !is.na(over(pts, wrld_simpl)$FIPS)
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image.plot(lon, lat, inout)
# plot(shp_ll, add = TRUE)
for (it in 1:dim(BA)[3]) {
  BA[,,it]=BA[,,it]*inout
}


image.plot(lon, lat, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

#from hectare to square meter
BA=BA*10000

#set entire period                                                                                                                                                                                           

# export as RData
save(BA, file = paste0(dir_out, "BA_",region,"_v1.RData"))
save(lon, file = paste0(dir_out, "lon_BA_",region,"_v1.RData"))
save(lat, file = paste0(dir_out, "lat_BA_",region,"_v1.RData"))

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
ncfname <- paste0(dir_out, "/Monthly_Burned_area_",start_year,"_",end_year,"_",region,"_v1.nc")
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
  'An integrated dataset of observed national inventories of fires - EUROPE'
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
#   "An integrated dataset of observed national inventories of fires"
# )

# close the file, writing data to disk
nc_close(ncout)


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

