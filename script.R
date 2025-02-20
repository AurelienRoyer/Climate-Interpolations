################# Supplementary material for the work Late Pleistocene temperature patterns 
##in the Western Palearctic: insights from rodent associations compared with General Circulation Models
## by Aurélien Royer, Julien Crétat, Rémi Laffont, Sara Gamboa, Belén Luna, Iris Menéndez, 
##Benjamin Pohl, Sophie Montuire, Manuel Hernandez Fernandez 
## version 20/02/2025

## settings the file with the data ----
setwd(".../suppl.R.script")

## load of the packages ----
library(IsoriX)
library(raster)
library(measurements)
library(MASS) 
library(rgdal)
library (tidyverse)
library(terra)

## load of the required functions ----
source("Utils_Iso.R")
source("PalBerV2.R")

## load bathimetrie files ----
test2<-load("finalMap_bathi_YD.rda")
finalMap_bathi_YD<-eval(as.symbol(test2))
test2<-load("finalMap_bathi_ALRO.rda")
finalMap_bathi_ALRO<-eval(as.symbol(test2))
test2<-load("finalMap_bathi_BOL.rda")
finalMap_bathi_BOL<-eval(as.symbol(test2))
test2<-load("finalMap_bathi_HE.rda")
finalMap_bathi_HE<-eval(as.symbol(test2))
test2<-load("finalMap_bathi_LGM.rda")
finalMap_bathi_LGM<-eval(as.symbol(test2))
liste.carte.bathi<-list(finalMap_bathi_YD,finalMap_bathi_ALRO,finalMap_bathi_BOL,finalMap_bathi_HE,finalMap_bathi_LGM)


## load ice sheet files ----
test2<-load("finalMap_ice_YD.rda")
finalMap_ice_YD<-eval(as.symbol(test2))
test2<-load("finalMap_ice_ALRO.rda")
finalMap_ice_ALRO<-eval(as.symbol(test2))
test2<-load("finalMap_ice_BOL.rda")
finalMap_ice_BOL<-eval(as.symbol(test2))
test2<-load("finalMap_ice_HE.rda")
finalMap_ice_HE<-eval(as.symbol(test2))
test2<-load("finalMap_ice_LGM.rda")
finalMap_ice_LGM<-eval(as.symbol(test2))
liste.carte.ice<-list(finalMap_ice_YD,finalMap_ice_ALRO,finalMap_ice_BOL,finalMap_ice_HE,finalMap_ice_LGM)


## load altitude maps ----
ElevEurope3<-rast("ElevEurope.tif")
ElevEurope.modern3<-rast("ElevEurope.modern.tif")


## To Calculate BCI from faunal list and estimate climate interpolation: example of Peyrazet level 5 ----
EUL="FALSE"
Peyrazet_5<-list(c("Arvicola_amphibius","Microtus_arvalis","Microtus_agrestis","Alexandromys_oeconomus","Apodemus_sylvaticus","Eliomys_quercinus","Sicista_betulina","Arvicola_sapidus","Chionomys_nivalis","Spermophilus_sp."))
Peyrazet_5_BCI<-Func_BCI_Calcul(Peyrazet_5, EUL = TRUE , verif = TRUE)
Peyrazet_5_interpo<-func_LDA(Peyrazet_5_BCI)

## Load csv files with climate inference to generate Spatial interpolation of past periods ----
LGM.iso.csv<-read.csv("LGM.iso.csv",sep=";")
HE.iso.csv<-read.csv("HE.iso.csv",sep=";")
BOL.iso.csv<-read.csv("BOL.iso.csv",sep=";")
ALRO.iso.csv<-read.csv("ALRO.iso.csv",sep=";")
YD.iso.csv<-read.csv("YD2.iso.csv",sep=";")

period.names<-c("YD","Allerod","Bolling","HS1","LGM")


## Spatial interpolation: example with LGM ----
liste.essai<-c("LGM.iso.csv")

period<-get(liste.essai[1])
for (i in 1: nrow(period)) {
  
  alt<-period[i,"alt"]
  cor_alt<-alt/100*-0.5
  period[i,"MAT.moy"]<-period[i,"MAT.moy"]-cor_alt
  period[i,"MAT.min"]<-period[i,"MAT.min"]-cor_alt
  period[i,"MAT.max"]<-period[i,"MAT.max"]-cor_alt
  period[i,"Tmax.moy"]<- period[i,"Tmax.moy"]-cor_alt
  period[i,"Tmax.min"]<-period[i,"Tmax.min"]-cor_alt
  period[i,"Tmax.max"]<-period[i,"Tmax.max"]-cor_alt
  period[i,"Tmin.moy"]<-period[i,"Tmin.moy"]-cor_alt
  period[i,"Tmin.min"]<-period[i,"Tmin.min"]-cor_alt
  period[i,"Tmin.max"]<-period[i,"Tmin.max"]-cor_alt
}
#write.csv(period,"LGM.iso_altCOR.csv") 

#calculate the GMML with the file having altitude effect removed
liste.essai<-c("LGM.iso_altCOR.csv")
iso.mat <- Compute_Iso2(liste.essai, ElevEurope3, par2use = "MAT")
iso.max <- Compute_Iso2(liste.essai, ElevEurope3, par2use = "Tmax")
iso.min <- Compute_Iso2(liste.essai, ElevEurope3, par2use = "Tmin")

liste.iso<-c("iso.mat")

# reinjection the altitude effect

  iso<-get(liste.iso[[1]])
  for (ii in 1: length(liste.essai)){
   iso[[ii]][[1]]$isoscapes$mean_corALT<-iso[[ii]][[1]]$isoscapes$mean+(ElevEurope3/100*-0.5)

  }
  
# to save the file
  #period1<-gsub(".csv","",paste0(liste.essai[ii]))
  #saveRDS_IsoriX(iso[[ii]][[1]], file = paste0("Europe_spatial_interp_",liste.iso[[1]],"_",period1,".rds"), compress = "xz") 
  #iso.mat<-readRDS(paste0("Europe_spatial_interp_",liste.iso[[1]],"_",period1,".rds"))
 
  #line to use only if you do not load the file with the previous line. 
  iso<- iso[[1]][[1]]
  #If you use readRDS, no need to use the previous line, and use this one : 
  iso<-iso.mat
 
  
## plot the spatial interpolation by generating a png file ----
  
 #preparation of symbols
  biome.proba.cex<-LGM.iso.csv$proba_value
  biome.proba.cex[biome.proba.cex > 0.95]<-2.2
  biome.proba.cex[biome.proba.cex < 0.95]<-1.5
  
  #preparation of the parameter range
  inf.range<-c(-50)
  sup.range<-c(40)
  nlabel<-c(19)
  nstep<-c(1)
  
  ## prepare the colors for the map
  isopal.mat.all<-c("#010782","#001C98", "#003D98", "#006198","#006D9E","#007DA7","#008FB3","#59A5C1","#97BFD2",
                    "#DEE6EA","#C8977F","#B77A55","#A66025","#964900","#812900","#661502","#592b27")
  
  ii=5 ## ii=1: YD, 2= Allerod, 3=Bolling, 4=HS1, 5=LGM ### to use period.names, liste.carte.bathi and liste.carte.ice. 
  title=paste0("Interpolation_Alt.cor", "_MAT", " - " ,period.names[[ii]])
  png(file=paste(title,".png"),width =1139, height = 742)
  plot_ISOSCAPE(iso, which = "mean_corALT",
                y_title = list(which = FALSE, title=title),
                sources = list(pch=17, cex=biome.proba.cex, col=attr(iso,which='col'),lwd=2, draw=TRUE),
                borders = list(borders = NA, lwd = 0.5, col = "white"),
                mask    = list(list(mask = vect(liste.carte.bathi[[ii]]),lwd=0.5, 
                                    col = "black", fill = "darkgrey"),
                               list(mask = vect(liste.carte.ice[[ii]]),lwd=0.5, 
                                    col = "black",  fill = "paleturquoise1")),
                palette = list(range = c(inf.range,
                                         sup.range), 
                               step = nstep,
                               n_labels = nlabel, 
                               fn =grDevices::colorRampPalette(isopal.mat.all, bias = 1), digits = 2)
                )
  dev.off()
  
  
## Extraction of LGM from a point ----

  r <- rast()
  crs(r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  data_loca<-data.frame(c(5.016667,-2.48794),c(47.316667,39.896027)) ## longitude and latitude of Dijon and Madrid
  colnames(data_loca)<-c("Longitude","Latitude")
  points <- vect(data_loca, geom=c("Longitude", "Latitude"), crs=crs(r))
  values_points<- terra::extract(iso$isoscapes$mean_corALT,points)
  
  
  
  
  
## Load csv files with climate inference to generate Spatial interpolation of modern period ----
  modern.iso.csv<-read.csv("Modern157.estimation.iso.csv",sep=";")
  
## correction altitude
  
  period<-modern.iso.csv
  for (i in 1: nrow(period)) {
    
    alt<-period[i,"alt"]
    cor_alt<-alt/100*-0.5
    period[i,"MAT.moy"]<-period[i,"MAT.moy"]-cor_alt
    period[i,"MAT.min"]<-period[i,"MAT.min"]-cor_alt
    period[i,"MAT.max"]<-period[i,"MAT.max"]-cor_alt
    period[i,"Tmax.moy"]<- period[i,"Tmax.moy"]-cor_alt
    period[i,"Tmax.min"]<-period[i,"Tmax.min"]-cor_alt
    period[i,"Tmax.max"]<-period[i,"Tmax.max"]-cor_alt
    period[i,"Tmin.moy"]<-period[i,"Tmin.moy"]-cor_alt
    period[i,"Tmin.min"]<-period[i,"Tmin.min"]-cor_alt
    period[i,"Tmin.max"]<-period[i,"Tmin.max"]-cor_alt
  }
  #write.csv(period,"Modern157.estimation.iso_altCOR.csv")
  
  #load file without the three outliers :  birjan, bouarfa, asabad
  liste.essai<-c("Modern157.estimation.iso_altCOR-outlier2.csv")
  
  iso.modern.mat <- Compute_Iso2(liste.essai, ElevEurope.modern3, par2use = "MAT")
  iso.modern.max <- Compute_Iso2(liste.essai, ElevEurope.modern3, par2use = "Tmax")
  iso.modern.min <- Compute_Iso2(liste.essai, ElevEurope.modern3, par2use = "Tmin")
  
  ## correctif altitude
  liste.iso <-c("iso.modern.mat")
  for (jj in 1:length(liste.iso)) {
    iso<-get(liste.iso[[jj]])
    iso$`Modern157.estimation.iso_altCOR-outlier2.csv`$MAT$isoscapes$mean_corALT<-iso$`Modern157.estimation.iso_altCOR-outlier2.csv`$MAT$isoscapes$mean+(as(ElevEurope.modern3, "SpatRaster") /100*-0.5)
  }
  
  #saveRDS_IsoriX(iso[[1]][[1]], file = paste0("Europe_spatial_interp_","iso.modern.mat",".rds"), compress = "xz") 
  #line to use only if you do not load the file with the previous line. 
  iso<- iso[[1]][[1]]
  
  ## prepa symbols
  bio.clim<- read.csv("Modern157.estimation.iso_altCOR-outlier2.csv",sep=",")
  biome.proba.cex<-bio.clim$proba_value
  biome.proba.cex[biome.proba.cex > 0.95]<-1.7
  biome.proba.cex[biome.proba.cex < 0.95]<-1

  biome.reel<-modern.iso.csv$Modern_climatic_zone
  biome.compa<-NULL
  for (i in 1:length(biome.reel)){
    biome.compa[i]<-identical(biome.reel[i],bio.clim$First_predite[i])
    
  }
  sum(biome.compa, na.rm = TRUE)
  biome.compa[biome.compa == TRUE]<-17
  biome.compa[biome.compa == FALSE]<-15
  
  title=paste0("Interpolation_Alt.cor", "_MAT", " - " ,"Modern")
  png(file=paste(title,".png"),width =1139, height = 742)
  plot_ISOSCAPE(iso, which = "mean_corALT",
                y_title = list(which = FALSE, title=title),
                sources = list(pch=biome.compa, cex=biome.proba.cex, col=attr(iso,which='col'),lwd=2, draw=TRUE),
                borders = list(borders = NA, lwd = 0.5, col = "white"),
                mask=list(mask = NA,col = "black",fill='darkgrey'),
                palette = list(range = c(inf.range,sup.range), step = nstep,n_labels = nlabel, 
                               fn =grDevices::colorRampPalette(IsoriX::isopalette1, bias = 1), digits = 2))
  dev.off()
  
  
  
  
  
  ############# Anomaly Modern-LGM
  iso.mat.LGM<-readRDS("Europe_spatial_interp_iso.mat_LGM.iso_altCOR.rds")
  iso.mat.modern<-readRDS(paste0("Europe_spatial_interp_","iso.modern.mat",".rds"))
  iso.mat.modern$isoscapes<-terra::resample(iso.mat.modern$isoscapes,iso.mat.LGM$isoscapes)
  
  diff.carte.iso<-iso.mat.LGM
  diff.carte.iso.2<-iso.mat.LGM$isoscapes$mean_corALT-iso.mat.modern$isoscapes$mean_corALT
  diff.carte.iso$isoscapes$mean<-diff.carte.iso.2$mean_corALT
  
  
  inf.range<--25
  sup.range<-25
  nstep<-0.5
  nlabel<-51
  ii=5
  
  title=paste("Difference - MAT - Present-Day"," - ",period.names[[ii]])
  png(file=paste("Difference - MAT - Present-Day"," - ",period.names[[ii]],".png"),width =1139, height = 742)
  plot_ISOSCAPE(diff.carte.iso,which = "mean",
                        y_title = list(which = FALSE, title=title),
                        sources = list(pch=17, cex=1.2, col="gray32",lwd=2, draw=TRUE),
                        borders = list(borders = NA, lwd = 0.5, col = "white"),
                        mask    = list(
                          list(mask = vect(liste.carte.bathi[[ii]]),col = "black", fill = "darkgrey"),
                          list(mask = vect(liste.carte.ice[[ii]]),col = "black", fill = "paleturquoise1")
                        ),
                        palette = list(range = c(inf.range,sup.range), 
                                       step = nstep,n_labels = nlabel, 
                                       fn = grDevices::colorRampPalette(c("blue3","blue1","deepskyblue1", "white","orange", "red","red4"), bias=1), digits = 2))
  
  dev.off()
  
  
