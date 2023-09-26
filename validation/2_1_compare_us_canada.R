rm(list = ls())
graphics.off()
gc()

library(ncdf4)
library(fields)
library(maptools)
library(RColorBrewer)
library(viridis)
source('~/Dropbox/model/fire_database/script_def/validation/image_mask.R')

## fix parameters
zlim <- c(0, 1000)
my_breaks <- c(0,1,2.5,5,10,25,50,100,250,500,1000)
my_col <-rev(inferno(length(my_breaks)-1))
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

brk_prob <- seq(0, 1, length.out = 5)
col_prob <-
  (colorRampPalette(brewer.pal(length(brk_prob), "YlGn"))(length(brk_prob) -1))

# Defire coordinate system
CRS.new <-
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data(wrld_simpl)
wrld_simpl=spTransform(wrld_simpl,
            CRS(
              "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
            ))

dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Documents/dati/fire_climate_data/fire/"

# Create a function to multiply each year of BA with the inout mask
multiply_years <- function(year) {
  year * inout
}


# Load grid
fname <- file.path(dir_out, 'misc/land_sea_mask_1degree.nc4')
obs.nc <- nc_open(fname)
lon = obs.nc$dim$lon$vals 
lat = obs.nc$dim$lat$vals
ilon=which(lon>=-171.5 & lon<=-53.5)
ilat=which(lat>=19.5 & lat<=70.5)
lon_com=lon[ilon]
lat_com=lat[ilat]
rm(lon)
rm(lat)


#firecci51
load(paste0(dir_obs,"/BA_200101_202012_1degree_nat.RData"))
F51=BA[ilon,ilat,]
years_f51=2001:2020

#common period 1
years=1986:2020 #largest common period
icommon <- seq(as.Date("1986-01-01"), as.Date("2020-12-01"), by = "month")
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

## monthly to total annual
BA_ALL_1_y = array(data = NA, dim = c(length(lon_com), length(lat_com), length(years) ))
BA_ALL_2_y = array(data = NA, dim = c(length(lon_com), length(lat_com), length(years) ))
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    idx = which(!is.na(BA_ALL_1[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        BA_ALL_1_y[i,j , iyear] = sum(BA_ALL_1[i, j, i1:i2], na.rm = TRUE) 
      }
    }
    idx = which(!is.na(BA_ALL_2[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        BA_ALL_2_y[i,j , iyear] = sum(BA_ALL_2[i, j, i1:i2], na.rm = TRUE) 
      }
    }
  }
}

F51_y = array(data = NA, dim = c(length(lon_com), length(lat_com), length(years_f51) ))
for (i in 1:length(lon_com)) {
  for (j in 1:length(lat_com)) {
    idx = which(!is.na(BA_ALL_1[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years_f51)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        F51_y[i,j , iyear] = sum(F51[i, j, i1:i2], na.rm = TRUE) 
      }
    }
  }
}

lon=lon_com
lat=lat_com
points <- expand.grid(lon, lat)
pts = SpatialPoints(points)
proj4string(pts) <- CRS.new

#meanBA
toplot=apply(BA_ALL_1_y,c(1,2),mean,na.rm=T)
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"/figures/colorbar_mean_BA.eps"), width=11,height=8.5)
image.plot(lon,lat,(toplot), zlim = zlim, 
           smallplot = c(0.85, 0.9, .31, .75),
           legend.only = TRUE, 
           legend.width = 0.5,
           horizontal = F,
           col = my_col,
           lab.breaks= my_breaks,
           # axis.args = list(at = def_breaks, labels = brks),
           axis.args = list(cex.axis = 1, font = 1),
           legend.mar = 4.7#,
           # axis.args = list(at = def_breaks, labels ,cex.axis = 1, font = 1))
)
dev.off()

setEPS()
# BA_ALL_1 = array(data = NA, dim = c(length(lon_com),length(lat_com),length(years)*12)) #NBAC+FPA_FOD
# BA_ALL_2 = array(data = NA, dim = c(length(lon_com),length(lat_com),length(years)*12)) #NFDB+MTBS
postscript(paste0(dir_out,"figures/mean_BA_NBAC_FPA_FOD.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()


toplot = apply(BA_ALL_1_y,c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig2a.RData"))

# rm(list=ls())
# ls()
# dir_out = '~/Dropbox/model/fire_database/out_def/'
# load(paste0(dir_out,"ForFigures/ForFig2a.RData"))
# ls()

toplot=apply(BA_ALL_2_y,c(1,2),mean,na.rm=T)
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"figures/mean_BA_NFDB_MTBS.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()

toplot = apply(BA_ALL_2_y,c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig2b.RData"))

cor.test(
  as.vector(apply(BA_ALL_1_y, c(1, 2), mean, na.rm = T)),
  as.vector(apply(BA_ALL_2_y, c(1, 2), mean, na.rm = T)),
  use = "pairwise.complete.obs",
  method = "spearman"
)


#comparison only canada 
i_c=which(wrld_simpl$NAME=="Canada")
ii <- !is.na(over(pts, wrld_simpl[i_c,1]))
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image(lon,lat,inout)
plot(wrld_simpl, add = TRUE)

image(lon,lat,inout*(apply(BA_ALL_1_y,c(1,2),mean,na.rm=T)))
plot(wrld_simpl, add = TRUE)

# mean1=mean(inout*(apply(BA_ALL_1_y,c(1,2),mean,na.rm=T)),na.rm=T)
aux1 <- apply(BA_ALL_1_y, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean1=mean(aux2)

aux1 <- apply(BA_ALL_2_y, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean2=mean(aux2)

100*(round(mean2)-round(mean1))/round(mean2)


#comparison only us 
i_c=which(wrld_simpl$ISO2=="US")
plot(wrld_simpl[i_c,])
ii <- !is.na(over(pts, wrld_simpl[i_c,1]))
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image(lon,lat,inout)
plot(wrld_simpl, add = TRUE)

image(lon,lat,inout*(apply(BA_ALL_1_y,c(1,2),mean,na.rm=T)))
plot(wrld_simpl, add = TRUE)
image(lon,lat,inout*(apply(BA_ALL_2_y,c(1,2),mean,na.rm=T)))
plot(wrld_simpl, add = TRUE)

aux1 <- apply(BA_ALL_1_y, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean1=mean(aux2[years>=1992])

aux1 <- apply(BA_ALL_2_y, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean2=mean(aux2[years>=1992])

100*(round(mean2)-round(mean1))/round(mean2)




## with f51
toplot=apply(F51_y/1000000,c(1,2),mean,na.rm=T)
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"figures/mean_BA_F51.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()
summary(F51_y/1000000)

toplot = apply(F51_y/1000000,c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig2e.RData"))


toplot=apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)
summary(as.vector(toplot))
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"figures/mean_BA_NFDB_MTBS_F51.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()

toplot = apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig2d.RData"))


toplot=apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)
summary(as.vector(toplot))
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"figures/mean_BA_NBAC_FPA_FOD_F51.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()

toplot = apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig2c.RData"))



cor.test(
  as.vector(apply(F51_y/1000000,c(1,2),mean,na.rm=T)),
  as.vector(apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor.test(
  as.vector(apply(F51_y/1000000,c(1,2),mean,na.rm=T)),
  as.vector(apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor.test(
  as.vector(apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),
  as.vector(apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

#comparison only canada 
i_c=which(wrld_simpl$NAME=="Canada")
ii <- !is.na(over(pts, wrld_simpl[i_c,1]))
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image(lon,lat,inout)
plot(wrld_simpl, add = TRUE)



aux1 <- apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]], MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean1=mean(aux2)

aux1 <- apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]], MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean2=mean(aux2)

aux1 <- apply(F51_y/1000000, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)
mean2=mean(aux2)


#comparison only US 
i_c=which(wrld_simpl$ISO2=="US")
ii <- !is.na(over(pts, wrld_simpl[i_c,1]))
inout = ii
dim(inout) <- c(length(lon), length(lat))
inout[inout == 0] = NA
image(lon,lat,inout)
plot(wrld_simpl, add = TRUE)


# Apply the function to each year of BA_ALL_1_y
aux1 <- apply(BA_ALL_1_y, MARGIN = 3, multiply_years)
aux2=apply(aux1,c(2),sum,na.rm=T)

aux3 <- apply(BA_ALL_2_y, MARGIN = 3, multiply_years)
aux4=apply(aux3,c(2),sum,na.rm=T)

aux5 <- apply(F51_y, MARGIN = 3, multiply_years)
aux6=apply(aux5,c(2),sum,na.rm=T)/1000000

# Create a blank plot
plot(NULL, xlim = c(1992, 2020), ylim = c(0, max(aux2, aux4, aux6)), 
     xlab = "Years", ylab = "Values", type = "n")

# Add lines for each variable
lines(years[years >= 1992], aux2[years >= 1992], col = "red", type = "l")
lines(years[years >= 1992], aux4[years >= 1992], col = "blue", type = "l")
lines(years[years >= 2001], aux6, col = "green", type = "l")
aux2[years >= 1992]/aux4[years >= 1992]
# Add a legend
legend("topleft", legend = c("aux2", "aux4", "aux6"), col = c("red", "blue", "green"), lty = 1)

mean(aux2[years >= 1992])
mean(aux4[years >= 1992])

# mean1=mean(inout*(apply(BA_ALL_1_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),na.rm=T)
# mean2=mean(inout*(apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)),na.rm=T)
# mean3=mean(inout*(apply(F51_y/1000000,c(1,2),mean,na.rm=T)),na.rm=T)
# 100*(round(mean2)-round(mean1))/round(mean2)

mean(aux2[years >= 2001])
mean(aux4[years >= 2001])
mean(aux6)


## correlation
zlim <- c(0, 1)
my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
my_col <-
  (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
my_col = my_col[1:length(my_col) - 1]
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))


# toplot=apply(BA_ALL_2_y[,,(dim(BA_ALL_2_y)[3]-19):dim(BA_ALL_2_y)[3]],c(1,2),mean,na.rm=T)
# summary(as.vector(toplot))
# toplot[toplot>zlim[2]]=zlim[2]
# image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# plot(wrld_simpl, add = TRUE)




corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    # if (!is.na(toplot[i,j])) {
    # if (sum(BA_ALL_1_y[i, j,],na.rm=T) >0 & sum(BA_ALL_2_y[i, j,],na.rm=T) >0) {
    if (sum(BA_ALL_1[i, j,],na.rm=T) >0 & sum(BA_ALL_2[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_1[i, j,],
        BA_ALL_2[i, j,],
        use = "pairwise.complete.obs",
        alternative = "greater",
        method = "spearman"
        # method = "pearson"
      )
      corre[i, j] = dum$estimate
      pvalue[i, j] <- dum$p.value
    }
  }
}


corre2=corre
corre2[pvalue>0.05]=NA
mean(corre2,na.rm = TRUE)
image.plot(lon, lat, corre2)
plot(wrld_simpl, add = TRUE)

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_can_usa.eps"), width=11,height=8.5)
image_mask(
  lon,
  lat,
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


save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig4a.RData"))


# 1 with 2 short period
BA_ALL_1=BA_ALL_1[,,-(1:180)]
BA_ALL_2=BA_ALL_2[,,-(1:180)]
corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre





for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    if (sum(BA_ALL_1[i, j,],na.rm=T) >0 & sum(BA_ALL_2[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_1[i, j,],
        BA_ALL_2[i, j,],
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
image.plot(lon, lat, corre2)
plot(wrld_simpl, add = TRUE)

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_can_usa_sort_period.eps"), width=11,height=8.5)
image_mask(
  lon,
  lat,
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


save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig4b.RData"))


corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    if (sum(BA_ALL_1[i, j,],na.rm=T) >0 & sum(F51[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_1[i, j,],
        F51[i, j,],
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
image.plot(lon, lat, corre2)
plot(wrld_simpl, add = TRUE)

# toplot=corre2
# toplot[toplot<=my_breaks[1]]=my_breaks[1]
# iok=which(pvalue<=0.05)
# setEPS()
# postscript(paste0(dir_out,"figures/corr_can_usa_NBAC_FPA_FOD_f51.eps"), width=11,height=8.5)
# image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()
# summary(as.vector(corre2))

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_can_usa_NBAC_FPA_FOD_f51.eps"), width=11,height=8.5)
image_mask(
  lon,
  lat,
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

save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig4c.RData"))

corre = array(data = NA, dim = c(length(lon_com), length(lat_com)))
pvalue=corre
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    if (sum(BA_ALL_2[i, j,],na.rm=T) >0 & sum(F51[i, j,],na.rm=T) >0) {
      dum = cor.test(
        BA_ALL_2[i, j,],
        F51[i, j,],
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
image.plot(lon, lat, corre2)
plot(wrld_simpl, add = TRUE)

# toplot=corre2
# toplot[toplot<=my_breaks[1]]=my_breaks[1]
# iok=which(pvalue<=0.05)
# setEPS()
# postscript(paste0(dir_out,"figures/corr_can_usa_NFDB_MTBS_f51.eps"), width=11,height=8.5)
# image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()

aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_can_usa_NFDB_MTBS_f51.eps"), width=11,height=8.5)
image_mask(
  lon,
  lat,
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

save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig4d.RData"))
