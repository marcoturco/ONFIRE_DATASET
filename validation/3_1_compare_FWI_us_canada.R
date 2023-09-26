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
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
FWI <- ncvar_get(obs.nc,"fwi") 
ilon=which(lon>=-171.5 & lon<=-53.5)
ilat=which(lat>=19.5 & lat<=70.5)
lon_com=lon[ilon]
lat_com=lat[ilat]
FWI=FWI[ilon,ilat,]
rm(lon)
rm(lat)
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

BA_ALL_1 = array(data = NA, dim = c(length(lon_com),length(lat_com),length(years)*12)) #NBAC+FPA_FOD
BA_ALL_2 = array(data = NA, dim = c(length(lon_com),length(lat_com),length(years)*12)) #NFDB+MTBS

# load grid 
region="CANADA_NBAC"
iNBAC <- seq(as.Date("1986-01-01"), as.Date("2020-12-01"), by = "month")
common_periods <- intersect(iNBAC, icommon)
indices_NBAC <- match(common_periods, iNBAC)
load(paste0(dir_out, "CANADA/BA_",region,"_v1.RData"))
load(paste0(dir_out, "CANADA/lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "CANADA/lat_BA_",region,"_v1.RData"))
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
BA_ALL_1[ilon,ilat,]=BA[,,indices_NBAC]*1e-6 #from m2 to km 2

region="US_FPA_FOD"
iFPA_FOD <- seq(as.Date("1992-01-01"), as.Date("2020-12-01"), by = "month")
common_periods <- intersect(iFPA_FOD, icommon)
indices_FPA_FOD <- match(common_periods, iFPA_FOD)
indices_ALL <- match(common_periods, iNBAC)
load(paste0(dir_out, "US/BA_",region,"_v1.RData"))
load(paste0(dir_out, "US/lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "US/lat_BA_",region,"_v1.RData"))
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    for (k in 1:length(indices_FPA_FOD)) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL_1[ilon[i], ilat[j], indices_ALL[k]] = sum(BA_ALL_1[ilon[i], ilat[j], k], BA[i, j, indices_FPA_FOD[k]], na.rm =
                                                           T)*1e-6 #from m2 to km 2
      }
    }
  }
}
image.plot(lon_com, lat_com, apply((BA_ALL_1), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)

# second merge
region="CANADA_NFDB"
iNFDB <- seq(as.Date("1959-01-01"), as.Date("2021-12-01"), by = "month")
common_periods <- intersect(iNFDB, icommon)
indices_NFDB <- match(common_periods, iNFDB)
load(paste0(dir_out, "CANADA/BA_",region,"_v1.RData"))
load(paste0(dir_out, "CANADA/lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "CANADA/lat_BA_",region,"_v1.RData"))
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
BA_ALL_2[ilon,ilat,]=BA[,,indices_NFDB]*1e-6 #from m2 to km 2

region="US_MTBS"
# iMTBS <- seq(as.Date("1984-01-01"), as.Date("2020-12-01"), by = "month") #but i only choose from 1992
iMTBS <- seq(as.Date("1992-01-01"), as.Date("2020-12-01"), by = "month") #but i only choose from 1992
common_periods <- intersect(iMTBS, icommon)
indices_MTBS <- match(common_periods, iMTBS)
load(paste0(dir_out, "US/BA_",region,"_v1.RData"))
load(paste0(dir_out, "US/lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "US/lat_BA_",region,"_v1.RData"))
# del period <1992
BA=BA[,,-(1:96)]
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    for (k in 1:length(indices_MTBS)) {
      if (!is.na(BA[i, j, k])) {
        BA_ALL_2[ilon[i], ilat[j], indices_ALL[k]] = sum(BA_ALL_2[ilon[i], ilat[j], k], BA[i, j, indices_MTBS[k]], na.rm =
                                                           T)*1e-6 #from m2 to km 2
      }
    }
  }
}
image.plot(lon_com, lat_com, apply((BA_ALL_2), c(1, 2), mean, na.rm = TRUE))
plot(wrld_simpl, add = TRUE)



## correlation
# zlim <- c(0, 1)
# my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
# my_col <-
#   (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
# my_col = my_col[1:length(my_col) - 1]
# def_breaks = seq(0,zlim[2],length.out=length(my_breaks))


corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    if (sum(BA_ALL_1[i, j,],na.rm=T) >0 & sum(FWI[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_1[i, j,],
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
# postscript(paste0(dir_out,"figures/corr_fwi_can_usa_1.eps"), width=11,height=8.5)
# image(lon_com,lat_com,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()
# summary(as.vector(corre2))

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_fwi_can_usa_1.eps"), width=11,height=8.5)
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

save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig6a.RData"))

corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    if (sum(BA_ALL_2[i, j,],na.rm=T) >0 & sum(FWI[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_2[i, j,],
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
# postscript(paste0(dir_out,"figures/corr_fwi_can_usa_2.eps"), width=11,height=8.5)
# image(lon_com,lat_com,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()
# summary(as.vector(corre2))

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_fwi_can_usa_2.eps"), width=11,height=8.5)
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

save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig6b.RData"))

########## con firecci51
#firecci51

fname<-file.path(dir_fwi, 'ERA5-FWI-1979-2021-IPCC-GRID-MONTHLY.nc')
obs.nc <- nc_open(fname)
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
FWI <- ncvar_get(obs.nc,"fwi") 
ilon=which(lon>=-171.5 & lon<=-53.5)
ilat=which(lat>=19.5 & lat<=70.5)
lon_com=lon[ilon]
lat_com=lat[ilat]
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
    if (sum(BA_ALL_2[i, j,],na.rm=T) >0 & sum(F51[i, j,],na.rm=T) >0 & sum(FWI[i, j,],na.rm=T) >0) {
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
# postscript(paste0(dir_out,"figures/corr_fwi_can_usa_f51.eps"), width=11,height=8.5)
# image(lon_com,lat_com,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()
# summary(as.vector(corre2))

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_fwi_can_usa_f51.eps"), width=11,height=8.5)
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
save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig6c.RData"))
