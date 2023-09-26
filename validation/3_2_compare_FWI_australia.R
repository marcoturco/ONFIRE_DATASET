rm(list = ls())
graphics.off()
gc()


# Package ‘boot.pval’
# https://online.stat.psu.edu/stat462/node/135/

library(ncdf4)
library(fields)
library(maptools)
library(rgdal)        
library(RColorBrewer)
source('~/Dropbox/model/fire_database/script_def/validation/image_mask.R')

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


data(wrld_simpl)

## fix parameters
# zlim <- c(0, 1)
# my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
# my_col <-
#   (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
# my_col = my_col[1:length(my_col) - 1]
# def_breaks = seq(0,zlim[2],length.out=length(my_breaks))
brk_prob <- seq(0, 1, length.out = 5)
col_prob <-
  (colorRampPalette(brewer.pal(length(brk_prob), "YlGn"))(length(brk_prob) -1))


dir_fwi='~/Documents/dati/obs/ERA5/FWI/'
dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Documents/dati/fire_climate_data/fire/"

fname<-file.path(dir_fwi, 'ERA5-FWI-1979-2021-IPCC-GRID-MONTHLY.nc')
obs.nc <- nc_open(fname)
obs.nc$dim$lon$vals -> lon_com
obs.nc$dim$lat$vals -> lat_com
FWI <- ncvar_get(obs.nc,"fwi") 
points <- expand.grid(lon_com, lat_com)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new
image.plot(lon_com,lat_com,apply(FWI,c(1,2),mean,na.rm=TRUE))
plot(wrld_simpl,add=TRUE)
years_fwi=1979:2021

#common period 1
years=1986:2020 #largest common period
icommon <- seq(as.Date("1986-01-01"), as.Date("2020-12-01"), by = "month")
iFWI <- seq(as.Date("1979-01-01"), as.Date("2021-12-01"), by = "month")
common_periods <- intersect(iFWI, icommon)
indices_FWI <- match(common_periods, iFWI)
FWI=FWI[,,indices_FWI]

## AUSTRALIA
# load grid 
load(paste0(dir_out, "AUSTRALIA/BA_AUSTRALIA_v1.RData"))
load(paste0(dir_out, "AUSTRALIA/lon_BA_AUSTRALIA_v1.RData"))
load(paste0(dir_out, "AUSTRALIA/lat_BA_AUSTRALIA_v1.RData"))
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
lon_com=lon_com[ilon]
lat_com=lat_com[ilat]
FWI=FWI[ilon,ilat,]
iok <- seq(as.Date("1950-01-01"), as.Date("2021-12-01"), by = "month")
common_periods <- intersect(iok, icommon)
indices_ok <- match(common_periods, iok)
BA=BA[,,indices_ok]*1e-6 #from m2 to km 2



## correlation
zlim <- c(0, 1)
my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
my_col <-
  (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
my_col = my_col[1:length(my_col) - 1]
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))


corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    if (sum(BA[i, j,],na.rm=T) >0 & sum(FWI[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA[i, j,],
        FWI[i, j,],
        use = "pairwise.complete.obs",
        alternative = "greater",
        method = "spearman"
      )
      corre[i, j] = dum$estimate
      pvalue[i, j] <- dum$p.value
    }
  }
}
corre2=corre
corre2[pvalue>0.05]=NA
mean(corre2,na.rm = TRUE)
image.plot(lon_com, lat_com, corre2)
plot(wrld_simpl, add = TRUE)



# toplot=corre2
# toplot[toplot<=my_breaks[1]]=my_breaks[1]
# iok=which(pvalue<=0.05)
# points <- expand.grid(lon, lat)
# setEPS()
# postscript(paste0(dir_out,"figures/corr_fwi_australia.eps"), width=11,height=8.5)
# image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_fwi_australia.eps"), width=11,height=8.5)
image_mask(
  lon_com,
  lat_com,
  aux,
  col = col_prob,
  zlim = c(0, 1),
  crop.color = 'gray',
  axis.args = list(
    at = c(-0.2, 0.15, 0.47, 0.78, 1.09, 1.25),
    labels = c(brk_prob, 'Not sig.')
  )
)
dev.off()
summary(as.vector(corre2))

save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig7a.RData"))


########## con firecci51
#firecci51

fname<-file.path(dir_fwi, 'ERA5-FWI-1979-2021-IPCC-GRID-MONTHLY.nc')
obs.nc <- nc_open(fname)
FWI <- ncvar_get(obs.nc,"fwi") 
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
FWI=FWI[ilon,ilat,]
rm(lon)
rm(lat)


#common period 1
years=2001:2020 #largest common period
icommon <- seq(as.Date("2001-01-01"), as.Date("2020-12-01"), by = "month")
iFWI <- seq(as.Date("1979-01-01"), as.Date("2021-12-01"), by = "month")
common_periods <- intersect(iFWI, icommon)
indices_FWI <- match(common_periods, iFWI)
FWI=FWI[,,indices_FWI]

load(paste0(dir_obs,"/BA_200101_202012_1degree_nat.RData"))
F51=BA
F51=F51[ilon,ilat,]

corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    if (sum(F51[i, j,],na.rm=T) >0 & sum(FWI[i, j,],na.rm=T) >0) {
      dum = cor.test(
        F51[i, j,],
        FWI[i, j,],
        use = "pairwise.complete.obs",
        alternative = "greater",
        method = "spearman"
      )
      corre[i, j] = dum$estimate
      pvalue[i, j] <- dum$p.value
    }
  }
}
corre2=corre
corre2[pvalue>0.05]=NA
mean(corre2,na.rm = TRUE)
image.plot(lon_com, lat_com, corre2)
plot(wrld_simpl, add = TRUE)

# toplot=corre2
# toplot[toplot<=my_breaks[1]]=my_breaks[1]
# iok=which(pvalue<=0.05)
# setEPS()
# postscript(paste0(dir_out,"figures/corr_fwi_aus_f51.eps"), width=11,height=8.5)
# image(lon_com,lat_com,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()
# summary(as.vector(corre2))

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_fwi_aus_f51.eps"), width=11,height=8.5)
image_mask(
  lon_com,
  lat_com,
  aux,
  col = col_prob,
  zlim = c(0, 1),
  crop.color = 'gray',
  axis.args = list(
    at = c(-0.2, 0.15, 0.47, 0.78, 1.09, 1.25),
    labels = c(brk_prob, 'Not sig.')
  )
)
dev.off()
summary(as.vector(corre2))


lon=lon_com
lat=lat_com
save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig7d.RData"))
