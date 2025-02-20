
################## Beyer2020 data -----
#Beyer, R. M., Krapp, M., & Manica, A. (2020). High-resolution terrestrial climate, bioclimate and vegetation for the last 120,000 years. Scientific data, 7(1), 236.

library(tidyr)
require("ncdf4")
require("lattice")
require("ggplot2")
library(raster)
library(rasterVis)

file <- "LateQuaternary_Environment.nc" ## Data file to find in the paper Beyer et al., 2020 #https://figshare.com/articles/dataset/LateQuaternary_Environment_nc/12293345/4

env_nc      <- ncdf4::nc_open(file)
longitude   <- ncdf4::ncvar_get(env_nc, "longitude")
latitude    <- ncdf4::ncvar_get(env_nc, "latitude")
years       <- ncdf4::ncvar_get(env_nc, "time")
months      <- ncdf4::ncvar_get(env_nc, "month")
temperature <- ncdf4::ncvar_get(env_nc, "temperature")
biome       <- ncdf4::ncvar_get(env_nc, "biome")
ncdf4::nc_close(env_nc) 

##extract raster de Beyer2020 at -21000
my_year      <- -21000

lat2 <- seq(min(latitude),max(latitude),0.5)
e <- extent(min(longitude)-0.25, max(longitude)+0.25, min(latitude)-0.25, max(latitude)+0.25)
r <- raster(nrow=length(lat2), ncol=length(longitude), ext=e)
m <- matrix(NA, nrow=length(lat2), ncol=length(longitude))

#extraction per month and transformation
v1<-temperature[,,months == 1,years == my_year]
v2<-temperature[,,months == 2,years == my_year]
v3<-temperature[,,months == 3,years == my_year]
v4<-temperature[,,months == 4,years == my_year]
v5<-temperature[,,months == 5,years == my_year]
v6<-temperature[,,months == 6,years == my_year]
v7<-temperature[,,months == 7,years == my_year]
v8<-temperature[,,months == 8,years == my_year]
v9<-temperature[,,months == 9,years == my_year]
v10<-temperature[,,months == 10,years == my_year]
v11<-temperature[,,months == 11,years == my_year]
v12<-temperature[,,months == 12,years == my_year]


v1 <- t(v1[,ncol(v1):1])
v2 <- t(v2[,ncol(v2):1])
v3 <- t(v3[,ncol(v3):1])
v4 <- t(v4[,ncol(v4):1])
v5 <- t(v5[,ncol(v5):1])
v6 <- t(v6[,ncol(v6):1])
v7 <- t(v7[,ncol(v7):1])
v8 <- t(v8[,ncol(v8):1])
v9 <- t(v9[,ncol(v9):1])
v10 <- t(v10[,ncol(v10):1])
v11 <- t(v11[,ncol(v11):1])
v12 <- t(v12[,ncol(v12):1])
r1<-r
r2<-r
r3<-r
r4<-r
r5<-r
r6<-r
r7<-r
r8<-r
r9<-r
r10<-r
r11<-r
r12<-r
r1$layer<-v1
r2$layer<-v2
r3$layer<-v3
r4$layer<-v4
r5$layer<-v5
r6$layer<-v6
r7$layer<-v7
r8$layer<-v8
r9$layer<-v9
r10$layer<-v10
r11$layer<-v11
r12$layer<-v12

r.beyer.max<-max(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12)
r.beyer.min<-min(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12)
r.beyer.mat<-mean(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12)
# save(r.beyer.max,file="beyer_raster_12000_max_temperature.Rdata")
# save(r.beyer.min,file="beyer_raster_12000_min_temperature.Rdata")
