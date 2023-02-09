rm(list = ls())
graphics.off()
gc()

library(rgdal)
library(R.matlab)
library(ncdf4)
library(maptools)
library(fields)
library(cshapes)
library(sf)
library(sp)
library(raster)
library(rio)        
library(stringr)


# Define name of the region for the outputs
region="CHILE"

# Define time limits
start_year <- 1985
end_year <- 2021

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Load world limits
data(wrld_simpl)

# Load country limits
cshp <-
  cshp(as.Date("2019-1-1"),
       dependencies = FALSE)
ic=which(cshp$gwcode==155) # 155 is the gwcode of Chile
shp_ll = as_Spatial(cshp[ic,1]) 
proj4string(shp_ll) <- CRS.new

# Define data directories
dir_out = '/Users/marco/Dropbox/model/fire_database/out_def/'
dir_obs <- "/Users/marco/Documents/dati/fire_climate_data/fire/"

# Load fire data
my_data <-
  import_list(
    paste0(dir_obs,
      '/chile/1671211598TABLA6_TEMPORADA2021_06b_12.09.2022_version2022_DICIEMBRE.xls'
    )
  )


num_reg_max=16
names_reg=c("XV","I","II","III","IV","V","RM","VI","VII","XVI","VIII","IX","XIV","X","XI","XII")
years=1984:2022
CHILE = array(data = NA, dim = c(length(years) *
                                  12,num_reg_max))

k=0
for (iyear in years[1]:(years[length(years)]-1)) {
  k=k+1
  i1 = (k - 1) * 12 + 7
  i2 = (k - 1) * 12 + 12 + 6
  dum=dim(my_data[[paste0('', iyear+1, '')]])
  num_reg=dum[2]-3
  for (ireg in 1:num_reg) {
    aux = my_data[[paste0('', iyear+1, '')]][[paste0('', '...', ireg+1, '')]]
    name_reg = aux[8]
    iok_reg=which(names_reg==name_reg)
    if (length(which(is.na(as.numeric(aux[9:20]))))<12) { #if a entire 12-month period is without records, means that is NA, and not zero
      aux2=as.numeric(aux[9:20])
      aux2[is.na(aux2)]=0
      CHILE[i1:i2, iok_reg] = aux2
    }
  }
}

#from regions to 1ºx1º

# chile regions
# Fichero .shp
file_shp = paste0(dir_obs,'/chile/Regiones/Regional.shp')
shp <- readOGR(file_shp)
shp_ll <-
  spTransform(shp,
              CRS(
                "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
              ))
aux = shp_ll@bbox
joint_reg=c(1,2,3,16,15,4,5,6,7,14,13,12,11,10,9,8)

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

# Calculate BA for each grid cell
BA = array(data = NA, dim = c(length(lon01), length(lat01), length(years) *
                                12))

for (ireg in 1:length(joint_reg)) {
  print(paste0('ireg ',
               ireg,
               '/',
               length(joint_reg)))
  ii <- !is.na(over(pts, shp_ll[ireg, ]))
  inout = ii[, 1]
  dim(inout) <- c(length(lon01), length(lat01))
  inout[inout == 0] = NA
  # image.plot(lon01, lat01, inout)
  # plot(shp_ll[ireg, ], add = TRUE)
  # plot(wrld_simpl, add = TRUE)
  
  if (sum(inout, na.rm = TRUE) > 0) {
    print(ireg)
    for (i in 1:length(lon01)) {
      for (j in 1:length(lat01)) {
        if (!is.na(inout[i, j])) {
          BA[i, j , ] =
            (CHILE[,joint_reg[ireg]]/ sum(inout, na.rm = TRUE))
        }
      }
    }
  }
}

# chech the dataset
summary(as.vector(BA))
image.plot(lon01, lat01, apply(log(BA+1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

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


#from hectare to square meter
BA=BA*10000

#set entire period
BA=BA[,,-(1:12)] #delete 1984 as it is not complete
BA=BA[,,-((dim(BA)[3]-11):dim(BA)[3])] #delete 2022 as it is not complete

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
  'An integrated dataset of observed national inventories of fires - CANADA'
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