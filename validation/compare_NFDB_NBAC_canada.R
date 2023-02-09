rm(list = ls())
graphics.off()
gc()

library(ncdf4)
library(fields)
library(maptools)
library(RColorBrewer)
library(viridis)
library(RobustLinearReg)
library(wql) #The following objects are masked from ‘package:TTAinterfaceTrendAnalysis’:

## fix parameters
years=1986:2020
zlim <- c(0, 1)
my_breaks <- seq(zlim[1],zlim[2], length.out = 5)
my_col <-
  (colorRampPalette(brewer.pal(length(my_breaks), "YlGnBu"))(length(my_breaks)))
my_col = my_col[1:length(my_col) - 1]
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

data(wrld_simpl)

dir_out = '~/Dropbox/model/fire_database/out_def/'

## fire from nfdb
region="CANADA"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
NFDB=BA #1959:2021
lat_NFDB=lat
lon_NFDB=lon

region="CANADA_NBAC"
load(paste0(dir_out, "BA_",region,"_v1.RData"))
load(paste0(dir_out, "lon_BA_",region,"_v1.RData"))
load(paste0(dir_out, "lat_BA_",region,"_v1.RData"))
NBAC=BA #1986:2020
lat_NBAC=lat
lon_NBAC=lon

# iok=match(lon_NFDB,lon_NBAC)

#set same temporal period
NFDB=NFDB[,,-((dim(NFDB)[3]-11):dim(NFDB)[3])] #delete 2021
NFDB=NFDB[,,(dim(NFDB)[3]-dim(NBAC)[3]+1):dim(NFDB)[3]] #delete years<1986

## monthly to total annual
NFDBY = array(data = NA, dim = c(length(lon_NFDB), length(lat_NFDB), length(years) ))
NBACY=NFDBY

for (i in 1:length(lon_NFDB)) {
  for (j in 1:length(lat_NFDB)) {
    idx = which(!is.na(NFDB[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        NFDBY[i,j , iyear] = sum(NFDB[i, j, i1:i2], na.rm = TRUE) 
      }
    }
  }
}
for (i in 1:length(lon_NFDB)) {
  for (j in 1:length(lat_NFDB)) {
    idx = which(!is.na(NBAC[i, j, ]))
    if (length(idx) >= 1) {
      for (iyear in 1:length(years)) {
        i1 = (iyear - 1) * 12 + 1
        i2 = (iyear - 1) * 12 + 12
        NBACY[i,j , iyear] = sum(NBAC[i, j, i1:i2], na.rm = TRUE) 
      }
    }
  }
}


#correlation
corre = array(data = NA, dim = c(length(lon), length(lat)))
pvalue=corre
for (i in 1:length(lon_NFDB)) {
  for (j in 1:length(lat_NFDB)) {
    ilon = which(lon_NFDB[i] == lon_NBAC)
    ilat = which(lat_NFDB[j] == lat_NBAC)
    if (sum(NFDBY[i, j,],na.rm=T) >0 & sum(NBACY[ilon, ilat,],na.rm=T) >0) {
      dum = cor.test(
        NFDBY[i, j,],
        NBACY[ilon, ilat,],
        use = "pairwise.complete.obs",
        alternative = "greater",
        method = "pearson"
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


#corre[pvalue>0.05]=NA

setEPS()
postscript(paste0(dir_out,"colorbar_corr_nfdb_nbac.eps"), width=11,height=8.5)
image.plot(lon,lat,corre2, zlim = zlim, 
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
postscript(paste0(dir_out,"corr_nfdb_nbac.eps"), width=11,height=8.5)
image(lon,lat,(corre2), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()

#BIAS
zlim <- c(0, 1000)
my_breaks <- c(0,1,2.5,5,10,25,50,100,250,500,1000)
my_col <-rev(inferno(length(my_breaks)-1))
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

toplot=apply(NFDBY,c(1,2),mean,na.rm=T)*1e-6
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"mean_BA_NFDBY_canada.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()

toplot=apply(NBACY,c(1,2),mean,na.rm=T)*1e-6
toplot[toplot>zlim[2]]=zlim[2]
setEPS()
postscript(paste0(dir_out,"mean_BA_NBACY_canada.eps"), width=11,height=8.5)
image(lon,lat,(toplot), breaks = my_breaks, col = my_col, legend = FALSE, zlim=zlim)
plot(wrld_simpl, add = TRUE)
dev.off()


image.plot(lon, lat, (apply(NFDBY,c(1,2),mean,na.rm=T)*1e-6))
plot(wrld_simpl, add = TRUE)
image.plot(lon, lat, apply(NBACY,c(1,2),mean,na.rm=T)*1e-6)
plot(wrld_simpl, add = TRUE)

100*(mean((apply(NFDBY,c(1,2),mean,na.rm=T)*1e-6),na.rm=T)-mean((apply(NBACY,c(1,2),mean,na.rm=T)*1e-6),na.rm=T))/mean((apply(NFDBY,c(1,2),mean,na.rm=T)*1e-6),na.rm=T)
# [1] 26.6068

# trend

#map

sid1=apply(NFDBY,c(3),sum,na.rm=T)
sid2=apply(NBACY,c(3),sum,na.rm=T)

plot(sid1)
lines(sid2,col="red")
lines(sid1,col="blue")
dum = mannKen(sid2)

changes <- array(data = NA, dim = c(length(lon), length(lat)))
pval   <- array(data = NA, dim = c(length(lon), length(lat)))

for (i in 1:length(lon_NFDB)) {
  for (j in 1:length(lat_NFDB)) {
    if (sum(NFDBY[i, j,],na.rm=T) >0 ) {
      
      dum = mannKen(NFDBY[i,j,]*1e-6)
      # changes[i,j] = (100*(dum$sen.slope*length(years))/mean(NFDBY[i,j,],na.rm=T))
      changes[i,j] = dum$sen.slope*10
      pval[i,j] = dum$p.value
      rm(dum)
    }
  }
}


image.plot(lon,lat,changes)
plot(wrld_simpl,add = TRUE)


maxFWI=10
data2plot=changes
data2plot[data2plot>maxFWI]=maxFWI
data2plot[data2plot< -maxFWI]=-maxFWI

# pdf(file = paste0(dir_out,"MyPlot.pdf"),   # The directory you want to save the file in
#     width = 4, # The width of the plot in inches
#     height = 4) # The height of the plot in inches

zlim <- c(-maxFWI, maxFWI)
my_breaks <- seq(zlim[1],zlim[2], length.out = 11)
my_col <- rev(colorRampPalette(brewer.pal(11, "BrBG"))(10))
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

# pdf(paste0(dir_out,"mean_fwi_europe.pdf"), width=11,height=8.5)
image.plot(lon,lat,(data2plot), zlim = zlim, 
           # smallplot = c(0.85, 0.9, .31, .75),
           legend.only = F, 
           legend.width = 0.5,
           horizontal = F,
           col = my_col,
           lab.breaks= my_breaks,
           # axis.args = list(at = def_breaks, labels = brks),
           axis.args = list(cex.axis = 1, font = 1),
           legend.mar = 4.7#,
           # axis.args = list(at = def_breaks, labels ,cex.axis = 1, font = 1))
)
plot(wrld_simpl, add = TRUE)

iok=which(pval<=0.05)
points(ptn[iok,],pch = 20,cex=0.1)

