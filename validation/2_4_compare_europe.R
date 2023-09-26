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

dir_out = '~/Dropbox/model/fire_database/out_def/'
dir_obs <- "~/Documents/dati/fire_climate_data/fire/"


# Load grid
fname <- file.path(dir_out, 'misc/land_sea_mask_1degree.nc4')
obs.nc <- nc_open(fname)
lon = obs.nc$dim$lon$vals 
lat = obs.nc$dim$lat$vals
# ilon=which(lon>=-171.5 & lon<=-53.5)
# ilat=which(lat>=19.5 & lat<=70.5)
lon_com=lon
lat_com=lat
rm(lon)
rm(lat)



#firecci51
load(paste0(dir_obs,"/BA_200101_202012_1degree_nat.RData"))
F51=BA
years_f51=2001:2020
icommon <- seq(as.Date("2001-01-01"), as.Date("2020-12-01"), by = "month")

## EUROPE
# load grid 
load(paste0(dir_out, "EUROPE/BA_EUROPE_v1.RData"))
load(paste0(dir_out, "EUROPE/lon_BA_EUROPE_v1.RData"))
load(paste0(dir_out, "EUROPE/lat_BA_EUROPE_v1.RData"))
ilon=match(lon,lon_com)
ilat=match(lat,lat_com)
iok <- seq(as.Date("1980-01-01"), as.Date("2020-12-01"), by = "month")
common_periods <- intersect(iok, icommon)
indices_ok <- match(common_periods, iok)
BA=BA[,,indices_ok]*1e-6 #from m2 to km 2

## monthly to total annual
BA_y = array(data = NA, dim = c(length(lon), length(lat), length(years_f51) ))
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    idx = which(!is.na(BA[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years_f51)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        BA_y[i,j , iyear] = sum(BA[i, j, i1:i2], na.rm = TRUE) 
      }
    }
  }
}


F51_y = array(data = NA, dim = c(length(lon), length(lat), length(years_f51) ))
ino=array(data = NA, dim = c(length(lon), length(lat) ))
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    idx = which(!is.na(BA[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years_f51)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        F51_y[i,j ,iyear] = sum(F51[ilon[i], ilat[j], i1:i2], na.rm = TRUE) 
      }
      if (length(which(is.na(BA_y[i, j, ]) | BA_y[i, j, ]==0)) >= 1) {
        ino[i,j]=1
        F51_y[i,j ,which(is.na(BA_y[i, j, ]) | BA_y[i, j, ]==0)]=NA
      }
    }
  }
}

lonmin=max(lon)
lonmax=min(lon)
latmin=max(lat)
latmax=min(lat)
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    if (length(which(!is.na(BA[i,j,]))) >= 1) {
      lonmin=min(lonmin,lon[i])
      lonmax=max(lonmax,lon[i])
      latmin=min(latmin,lat[j])
      latmax=max(latmax,lat[j])
    }
  }
}


xlim <- c(min(lonmin)-1, max(lonmax)+1)
ylim <- c(min(latmin)-1, max(latmax)+1)


#meanBA
toplot=apply(BA_y,c(1,2),mean,na.rm=T)
summary(as.vector(apply(BA_y,c(1,2),mean,na.rm=T)))
toplot[toplot>zlim[2]]=zlim[2]

setEPS()
postscript(paste0(dir_out,"figures/mean_BA_EUROPE_ONFIRE.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, xlim=xlim, ylim=ylim,zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()


toplot = apply(BA_y,c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig3e.RData"))



## with f51
toplot=apply(F51_y/1000000,c(1,2),mean,na.rm=T)
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"figures/mean_BA_EUROPE_F51.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, xlim=xlim, ylim=ylim,zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()
summary(F51_y/1000000)


toplot = apply(F51_y/1000000,c(1,2),mean,na.rm=T)
save(list=c("lon", "lat", "toplot"), file=paste0(dir_out,"ForFigures/ForFig3f.RData"))


100*(mean((apply(BA_y,c(1,2),mean,na.rm=T)),na.rm=T)-mean((apply(F51_y/1000000,c(1,2),mean,na.rm=T)),na.rm=T))/mean((apply(BA_y,c(1,2),mean,na.rm=T)),na.rm=T)

plot(apply(BA_y,c(1,2),mean,na.rm=T),apply(F51_y/1000000,c(1,2),mean))


cor.test(
  as.vector(apply(BA_y,c(1,2),mean,na.rm=T)),
  as.vector(apply(F51_y/1000000,c(1,2),mean,na.rm=T)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

# mean1=mean((apply(BA_y,c(1,2),mean,na.rm=T)),na.rm=T)
# mean2=mean((apply(F51_y/1000000,c(1,2),mean,na.rm=T)),na.rm=T)
# 100*(round(mean2)-round(mean1))/round(mean2)
# 


aux2=apply(BA_y,c(3),sum,na.rm=T)
mean1=mean(aux2)


aux2=apply(F51_y/1000000,c(3),sum,na.rm=T)
mean2=mean(aux2)

plot.ts(years_f51,apply(BA_y,c(3),sum,na.rm=T))
lines(years_f51,apply(BA_y,c(3),sum,na.rm=T),col="black")
lines(years_f51,apply(F51_y/1000000,c(3),sum,na.rm=T),col="red")



## correlation


# source('~/Dropbox/model/script_climate_fire_v2/image_mask.R')
# 
# aux = best_mod * NA
# aux[mask == 1] = -999
# aux[!is.na(best_mod)] = best_mod[!is.na(best_mod)]
# 
# 
# 
# image_mask(
#   lon[ilon],
#   lat[ilat],
#   aux[ilon, ilat],
#   col = col_prob,
#   # zlim=c(0,1),crop.color='gray',axis.args = list(at = c(brk_prob,brk_prob[length(brk_prob)]+diff(brk_prob)[1]), labels=c(brk_prob,'No effect')))
#   zlim = c(1, 5),
#   crop.color = 'gray',
#   axis.args = list(
#     at = c(1, 2, 3, 4, 5, 6),
#     # at = c(0.5, 2.15, 3.47, 4.78, 5.09, 5.25),
#     labels = c(list_models, 'No effect')
#   )
# )<<<



zlim <- c(0, 1)
my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
my_col <-
  (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
my_col = my_col[1:length(my_col) - 1]
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

# 1 with 2 short period

corre = array(data = NA, dim = c(length(lon), length(lat)))
pvalue=corre
for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    
    BA_clean <- BA[i, j, ][complete.cases(BA[i, j, ])]
    F51_clean <- F51[ilon[i], ilat[j], ][complete.cases(F51[ilon[i], ilat[j], ])]
    n_obs <- length(BA_clean)
    
    # if (sum(BA[i, j,],na.rm=T) >0 & sum(F51[ilon[i], ilat[j],],na.rm=T) >0) {
      
    if (n_obs > 1) {
      dum = cor.test(
        BA[i, j,],
        F51[ilon[i], ilat[j],],
        use = "pairwise.complete.obs",
        alternative = "greater",
        # method = "pearson"
        method = "spearman"
      )
      corre[i, j] = dum$estimate
      pvalue[i, j] <- dum$p.value
    }
    else {
      # Handle the case where there are not enough observations
      print("Not enough finite observations for the correlation test.")
    }
  }
}
corre2=corre
corre2[pvalue>0.05]=NA
mean(corre2,na.rm = TRUE)
image.plot(lon, lat, corre)
plot(wrld_simpl, add = TRUE)

# toplot=corre2
# toplot[toplot<=my_breaks[1]]=my_breaks[1]
# iok=which(pvalue<=0.05)
# points <- expand.grid(lon, lat)
# setEPS()
# postscript(paste0(dir_out,"figures/corr_EUROPE_sort_period.eps"), width=11,height=8.5)
# image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
# # points(points[iok,],pch = 20,cex=1)
# plot(wrld_simpl, add = TRUE)
# dev.off()


# Set the desired limits for x and y axes


aux = corre2 * NA
aux[pvalue>0.05] = -999
aux[!is.na(corre2)] = corre2[!is.na(corre2)]
postscript(paste0(dir_out,"figures/corr_mask_EUROPE_sort_period.eps"), width=11,height=8.5)
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
  ),
  xlim = xlim, ylim = ylim
)
dev.off()


save(list=c("lon", "lat", "corre","pvalue"), file=paste0(dir_out,"ForFigures/ForFig5c.RData"))


