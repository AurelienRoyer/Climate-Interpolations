################# Supplementary material for the work Late Pleistocene temperature patterns 
##in the Western Palearctic: insights from rodent associations compared with General Circulation Models
## by Aurélien Royer, Julien Crétat, Rémi Laffont, Sara Gamboa, Belén Luna, Iris Menéndez, 
##Benjamin Pohl, Sophie Montuire, Manuel Hernandez Fernandez 
## version 20/02/2025

#### LOOCV script

csv<-list.files(pattern="iso_altCOR.csv", recursive = F)
## example with LGM
iso3<-read.csv(csv[4])

liste.essai<-NULL
# creation a file without one site 
  for (HH in 1: nrow(iso3)) {
    nom2<-iso3[HH,"couche"]
    LGM.iso2<-iso3[-c(HH),]
   lll<-sub(".csv","",csv[4])
    nom3<-paste0(lll,"-",HH,".csv")
    write.csv(LGM.iso2,file=paste0(nom3))
    
  }

liste.essai.loocv<-list.files(pattern=("LGM.iso_altCOR*"), recursive = F)

#computation of interpolation per file /// take time ...
for (i in 1:length(liste.essai.loocv)){
  
  iso <- Compute_Iso2(liste.essai.loocv[i], ElevEurope3, par2use = "MAT")
  gc()
  
}

# traitement SD 

liste.rds<-list.files(pattern="Europe_isoscape_VR_2024_", recursive = F)
df <- data.frame(filename=liste.rds, stringsAsFactors = FALSE)
for (i in 1:nrow(df)){
  test2<-readRDS(df[i,1])
  nom<-paste(df[i,1])
  nom<-gsub(".iso_altCOR","",nom)
  nom<-gsub(".rds","",nom)
  nom<-gsub(" ","",nom)
  nom<-gsub("Europe_isoscape_VR_2024_","LOOCV_",nom)
  nom<-gsub("-",".",nom)
  assign(nom, test2,envir = .GlobalEnv)
  df[i,2]<-nom
  
}

mylist<-df[,2]
LOOCv.bis<-NULL
temp<-get(mylist[1])
LOOCv.bis<-temp$isoscapes$mean
temp<-unwrap(temp$isoscapes)
LOOCv.bis<-temp$mean
for (i in 2:length(mylist)){
  temp<-get(mylist[i])
  temp<-unwrap(temp$isoscapes)
  LOOCv.bis<-c(LOOCv.bis,temp$mean)
}

#mylist.LGM<-mylist
sd.LGM<-terra::stdev(LOOCv.bis)
temp$isoscapes$mean<-sd.LGM$std
ii=5

title=paste0("SD-MAT-LGM")
png(file=paste(title,".png"),width =1139, height = 742)
plot_ISOSCAPE(temp, which = "mean",
              y_title = list(which = FALSE, title=title),
              sources = list(pch=17, cex=1.7, col=attr(temp,which='col'),lwd=2, draw=TRUE),
              borders = list(borders = NA, lwd = 0.5, col = "white"),
              mask    = list(list(mask = vect(liste.carte.bathi[[ii]]),lwd=0.5, col = "black", fill = "darkgrey"),
                             list(mask = vect(liste.carte.ice[[ii]]),lwd=0.5, col = "black",  fill = "paleturquoise1")),
              palette = list(range = c(0,2), step = 0.1,n_labels = 21, fn =grDevices::colorRampPalette(IsoriX::isopalette1, bias = 1), digits = 2))
dev.off()


