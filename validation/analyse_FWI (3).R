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

data(wrld_simpl)

## fix parameters
# dir_fwi='/Users/marco/Documents/dati/obs/ERA5/FWI/'
dir_data = '/Users/marco/Dropbox/estcena/scripts/alberto_tfg/datos/'
dir_out = '/Users/marco/Dropbox/estcena/scripts/alberto_tfg/out/'
year_fwi=1979:2021



#load data
fname<-file.path(dir_data, 'ERA5_1979-2021-1-yearly-europe.nc')
obs.nc <- nc_open(fname)
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
ptn <- expand.grid(lon, lat)
FWI <- ncvar_get(obs.nc,"fwi") 
ni = dim(FWI)[1]
nj = dim(FWI)[2]

# obs[obs=="-9999"] <- NA
# lon=lon-180
# # dum=obs
# # dum[1:(nrow(obs)/2),,]=obs[(nrow(obs)/2+1):nrow(obs),,]
# # dum[(nrow(obs)/2+1):nrow(obs),,]=obs[1:(nrow(obs)/2),,]
# # rm(obs)
# FWI=dum
# rm(dum)

image.plot(lon,lat,apply(FWI,c(1,2),mean,na.rm=TRUE))
plot(wrld_simpl,add=TRUE)

dim(FWI)



## annual mean FWIy
## to plot
maxFWI=100
data2plot=FWI
data2plot[data2plot>maxFWI]=maxFWI

# pdf(file = paste0(dir_out,"MyPlot.pdf"),   # The directory you want to save the file in
#     width = 4, # The width of the plot in inches
#     height = 4) # The height of the plot in inches

zlim <- c(0, maxFWI)
my_breaks <- seq(zlim[1],zlim[2], length.out = 6)
my_col <-rev(inferno(length(my_breaks)-1))
def_breaks = seq(0,zlim[2],length.out=length(my_breaks))

data2plot=apply(FWI,c(1,2),mean,na.rm=T)
data2plot[data2plot>maxFWI]=maxFWI

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
# dev.off()

## tendencia annual mean regional averaged

## compute global annual (or seasonal) averages
### show cosine of the latitude of each grid cell (which is used to createregional  mean FWI)
#image.plot(lon,lat,ww); plot(wrld_simpl,add = TRUE)
# main="Cosine of the latitude (used to areally weight cells)")
ww <- matrix(data = NA,nrow = ni,ncol = nj)
for (i in 1:ni) {
  for (j in 1:nj) {
    ww[i,j] = cos(lat[j] * pi / 180)
  }
}
image.plot(lon,lat,ww)
plot(wrld_simpl,add = TRUE)

dum <- matrix(FWI,ni * nj,length(year_fwi))
FWI_reg <- apply(dum,2,weighted.mean,w = c(ww / sum(ww)),na.rm = TRUE)

plot.ts(year_fwi,FWI_reg)



model <- theil_sen_regression(FWI_reg~year_fwi)
sid=summary(model)


#Plotting the results on the project step-by-spet
# plot(
#   FWI_reg~year_fwi,
#   col = "gray",
#   xlab = "x",
#   ylab = "y",
#   main = "Compare regressions")
# apply(coefs, 2, abline, col = rgb(1, 0, 0, 0.03))

# x=year_fwi
# new.data = seq(min(x), max(x), by = 1)
# conf_interval <-
#   predict(
#     model,
#     newdata = data.frame(x = new.data),
#     interval = "confidence",
#     level = 0.95)
# 
# lines(new.data, conf_interval[, 2], col = "black", lty = 3, lwd=3)
# lines(new.data, conf_interval[, 3], col = "black", lty = 3, lwd=3)


# abline(model,col='blue')

# lines(FWI_reg, conf_interval[, 2], col = "black", lty = 3, lwd=3)
# lines(FWI_reg, conf_interval[, 3], col = "black", lty = 3, lwd=3)
# legend("topleft",
#        legend = c("Bootstrap", "Population", 'Sample'),
#        col = c("red", "blue", 'green'),
#        lty = 1:3,
#        cex = 0.8)
# 







aux=mannKen(FWI_reg)

pval=character()
if (aux$p.value<0.05) {pval="*"}
perc_change=round(100*(aux$sen.slope*length(year_fwi))/mean(FWI_reg))

plot(
  FWI_reg~year_fwi,
  col = "gray",
  xlab = "years",
  ylab = "FWI",
  main = "Annual FWI")
lines(FWI_reg~year_fwi,col = "gray")
abline(model,col='red')
text(min(year_fwi)+3, max(FWI_reg)-0.01*max(FWI_reg), paste0(perc_change,"%",pval))

#map
changes <- matrix(data = NA,nrow = ni,ncol = nj)
pval   <- matrix(data = NA,nrow = ni,ncol = nj)

for (i in 1:ni) {
  for (j in 1:nj) {
    if (!is.na(sum(FWI[i,j,]))) {
      dum = mannKen(FWI[i,j,])
      changes[i,j] = round(100*(dum$sen.slope*length(year_fwi))/mean(FWI[i,j,]))
      pval[i,j] = dum$p.value
      rm(dum)
    }
  }
}


image.plot(lon,lat,changes)
plot(wrld_simpl,add = TRUE)

maxFWI=100
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