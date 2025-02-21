################# Supplementary material for the work Late Pleistocene temperature patterns 
##in the Western Palearctic: insights from rodent associations compared with General Circulation Models
## by Aurélien Royer, Julien Crétat, Rémi Laffont, Sara Gamboa, Belén Luna, Iris Menéndez, 
##Benjamin Pohl, Sophie Montuire, Manuel Hernandez Fernandez 
## version 20/02/2025

#### function used 

plot_ISOSCAPE <- function (x, which = "mean", y_title = list(which = TRUE, title = bquote(delta^2 * H)),
                           sources = list(draw = TRUE, cex = 0.5, pch = 2, lwd = 1, col = "red"),
                           borders = list(borders = NA, lwd = 0.5, col = "black"),
                           mask = list(mask = NA, lwd = 0, col = "black", fill = "black"),
                           palette = list(step = NA, range = c(NA, NA), n_labels = 11, digits = 2, fn = NA), 
                           plot = TRUE, sphere = list(build = FALSE, keep_image = FALSE), ...) 
{
  if (!(any(class(x) %in% "ISOSCAPE"))) {
    stop("This function must be called on an object of class ISOSCAPE.")
  }
  simu <- "ISOSIM" %in% class(x)
  #IsoriX:::.complete_args(IsoriX:::plot.ISOSCAPE)
  if (!is.null(palette$fn) && !is.function(palette$fn) && is.na(palette$fn)) {
    isopalette1 <- NULL
    utils::data("isopalette1", envir = environment(), 
                package = "IsoriX")
    palette$fn <- grDevices::colorRampPalette(isopalette1, 
                                              bias = 1)
  }
  if (!is.null(borders$borders) && is.na(borders$borders)) {
    CountryBorders <- NULL
    borders$borders <- terra::readRDS(system.file("extdata/CountryBorders.rds", # modif 2024/03
                                                  package = "IsoriX"))
    #borders$borders <- CountryBorders
  }
  if (!is.null(mask$mask) && !inherits(mask$mask, "SpatVector") && #class(mask$mask) != "SpatialPolygons"
      identical(NA, mask$mask)) { #is.na(mask$mask)  # modif 2024/03 : conditions dans if
    OceanMask <- NULL
        # mask$mask <- terra::readRDS(system.file("extdata/OceanMask.rds", 
        #     package = "IsoriX"))
    mask$mask <- terra::readRDS(system.file("extdata/OceanMask.rds",  # modif 2024/03
                                            package = "IsoriX"))    
  }
  if (simu) {
    if (sources$draw) {
      sources$draw <- FALSE
      message("You have asked to plot sources, but it does not make sense for simulations as each raster cell is a source. The argument 'plot.sources' was thus considered to be FALSE.")
    }
    if (!(which %in% c("mean", "disp"))) {
      stop("For simulated data, the argument 'which' must be 'mean' or 'disp'.")
    }
  }
  else {
    if (!(which %in% c("mean", "mean_predVar", 
                       "mean_residVar", "mean_respVar", "disp", 
                       "disp_predVar", "disp_residVar", "disp_respVar",
                       "mean_corALT"                      #modif 2024 - altitude
                       ))) {
      stop("argument 'which' unknown")
    }
  }
  
  #print(x$isoscape[[which]])
  #print(str(palette))
  
  colours <- IsoriX:::.cut_and_color(var = x$isoscape[[which]], step = palette$step, 
                                     range = palette$range, palette = palette$fn, n_labels = palette$n_labels, 
                                     digits = palette$digits)
  
  
  Title <- ""
  if (simu) 
    Title <- "simulated"
  if (y_title$which) 
    Title <- paste(Title, sub("_", " ", which, 
                              fixed = TRUE))
  if (!is.null(y_title$title)) 
    Title <- bquote(.(Title) ~ .(y_title$title))
  map <- rasterVis::levelplot(x$isoscape[[which]], maxpixels = prod(dim(x$isoscape[[which]])[1:2]), 
                              margin = FALSE, col.regions = colours$all_cols, at = colours$at, 
                              colorkey = list(labels = list(at = colours$at_keys, labels = colours$at_labels)), 
                              main = Title)
  
  #print(x$isoscape[[which]]@file@name)
  #print(str(x$isoscape[[which]]@data))
  #print(str(map$panel.args.common))
  
  ##############
  if (!is.null(mask$mask)){
    decor <- IsoriX:::.build_additional_layers(x = x, sources = sources, calibs = NULL, borders = borders, mask = mask)
    complete_map <- map + decor$borders_layer + decor$mask_layer +decor$sources_layer    
  }else{
    
    Masks <- mask
    
    n_masks <- length(Masks)
    for (i in 1:n_masks){
      
      mask <- Masks[[i]]
      
      if (!is.null(mask$mask) && !inherits(mask$mask, "SpatVector") && #class(mask$mask) != "SpatialPolygons"
          identical(NA, mask$mask)) { #is.na(mask$mask)  # modif 2024/03 : conditions dans if
        OceanMask <- NULL
        # utils::data("OceanMask", envir = environment(), 
        #             package = "IsoriX")
        mask$mask <- terra::readRDS(system.file("extdata/OceanMask.rds",   # modif 2024/03 
                                                package = "IsoriX"))           
        mask$mask <- OceanMask
      }
      
      
      if (i==1){
        
        decor <- IsoriX:::.build_additional_layers(x = x, sources = sources, calibs = NULL, borders = borders, mask = mask)
        complete_map <- map + decor$borders_layer + decor$mask_layer +decor$sources_layer
      }else{
        decor <- IsoriX:::.build_additional_layers(x = x, sources = sources, calibs = NULL, borders = borders, mask = mask)
        complete_map <- complete_map + decor$borders_layer + decor$mask_layer +decor$sources_layer
      }
      
      
      
    }
  }
  #############
  
  
  
  
  if (plot & !sphere$build) {
    
    if (IsoriX:::.data_IsoriX$IsoriX_options$dont_ask) {
      
      options(example_ask = "FALSE")
    }
    
    print(complete_map)
  }
  rm(Title)
  
  if (sphere$build) {
    IsoriX:::.build_sphere(x$isoscape[[which]], colours = colours, 
                           decor = decor)
    if (!sphere$keep_image) {
      file.remove("IsoriX_world_image.png")
    }
  }
  return(invisible(complete_map))
}


MaskFromRaster <- function(Rast, ext_lon, ext_lat, step_lon, step_lat, thr, testSup = TRUE){
  
  # mask_topo <- MaskFromRaster(topo, c(-15, 65), c(30, 70), 2, 2, -120, TRUE)
  # mask_ice <- MaskFromRaster(ice, c(-15, 65), c(30, 70), 2, 2, 50, FALSE)
  
  if (testSup){
    mask <- Rast > thr
  }else{
    mask <- Rast < thr
  }
  
  
  
  cpt <- 0
  pb <- txtProgressBar(min = cpt+1, max = (diff(ext_lat)*diff(ext_lon))
                       / (step_lon*step_lat), style=3)
  for (lon in seq(ext_lon[1], ext_lon[2] - step_lon, by=step_lon)){
    for (lat in seq(ext_lat[1], ext_lat[2] - step_lat, by=step_lat)){

      cpt <- cpt+1
      setTxtProgressBar(pb, cpt)
      
      
      
      e <- extent(c(lon, lon+step_lon, lat, lat+step_lat))
      
      mask_pol <-crop(mask, e)
      mask_pol<-rasterToPolygons(mask_pol, fun=function(x) {x< 0.5})
      
      if (!is.null(mask_pol)){
        mask_pol <- aggregate(mask_pol)
        
        if (!exists("final_Mask")){
          final_Mask <- mask_pol
          
        }else{
          #https://gis.stackexchange.com/questions/180682/merge-a-list-of-spatial-polygon-objects-in-r
          
          mask_pol@polygons[[1]]@plotOrder <- mask_pol@polygons[[1]]@plotOrder + length(final_Mask@polygons[[1]]@plotOrder)
          spl <- list(final_Mask, mask_pol)
          final_Mask <- do.call(bind, spl)
          final_Mask<-aggregate(final_Mask)
          
        }
      }
    }
  }
  
  print("\n")
  
  return(final_Mask)
  
} 

make_source_cols <- function(x, vn = 3:9, vc = c("Yellow", "Red", "orchid3",
                              "green3", "bisque1", "cyan3", "deepskyblue1")){
  vc[match(x, vn)]
}

Compute_Iso2 <- function(csv_filnames, ElevMap, pattern_mean = ".moy",
                         par2use = "all"){
  
  n_dat <- length(csv_filnames) 
  
  LIso <- list()
  for (ii in 1:n_dat){
    
    csv_filname <- csv_filnames[ii]
    Data <- read.csv(csv_filname)
    var_names <- colnames(Data)
    idx_par <- grep(pattern_mean, var_names)
    par_names <- gsub(pattern_mean, "", var_names[idx_par])
    if (!identical(par2use, "all")){
      idx2do <- match(par2use, par_names)
      idx2do <- idx2do[!is.na(idx2do)]
      if (length(idx2do)<1){
        stop(cat("par2use should take values within: ", par_names, collapse=""))
      }
      idx_par <- idx_par[idx2do]
      par_names <- par_names[idx2do]
    }
    
    source_ID <- Data$couche
    
    cpt <- 0
    lIso <- list()
    for (jj in 1:length(idx_par)){
      
      cpt <- cpt + 1
      
      mean_source_value <- Data[, idx_par[jj]]
      se_fit <- Data[, match(paste0(par_names[jj], "_se.fit"), var_names)]
      n_source_value <- Data[, match(paste0(par_names[jj], "_df"), var_names)]
      rs <- Data[, match(paste0(par_names[jj], "_residual.scale"), var_names)]
      lat <- Data$lat
      long <- Data$long
      elev <- Data$alt
      
      var_source_value <- se_fit ^ 2 + rs ^ 2
      
      
      
      data_isorix <- data.frame(source_ID, mean_source_value, var_source_value,
                                n_source_value, lat, long, elev)
      
      Fit <- isofit(data = data_isorix, mean_model_fix = list(elev = FALSE, lat_abs = TRUE))
      
      #modif 2024
      if (class(ElevMap)== "SpatRaster") {
        Iso <- isoscape(raster = ElevMap, isofit = Fit)
      } else {
        
        Iso <- isoscape(raster = as(ElevMap, "SpatRaster") , isofit = Fit)
      }
      vcolors <- make_source_cols(Data$First_predite)
      attr(Iso, "col") <- vcolors
      lIso[[cpt]] <- Iso
      period1<-gsub(".csv","",paste0(csv_filnames[ii]))
      save(Fit, file = paste("EuropeFit_VR_2024_",par_names[[jj]],"_",period1,".rda"), compress = "xz")
      saveRDS_IsoriX(Iso, file = paste("Europe_isoscape_VR_2024_",par_names[[jj]],"_",period1,".rds"), compress = "xz")
    }
    names(lIso) <- par_names

    
    LIso[[ii]] <- lIso
    
  }
  
  names(LIso) <- csv_filnames
  
  return(LIso)
  
}

Save_IsoMap <- function(LIso, Lmasks = NULL, maskCol = c("darkgrey", "paleturquoise1"),
                        which = "mean", wd = getwd(), useCSVNames = TRUE,
                        PeriodNames = NULL, absColBar = TRUE, ncolors = 100, 
                        nlabs = ncolors / 5, width = 1139, height = 742,
                        y_title = list(which = FALSE), 
                        sources = list(pch = 17, cex = 1, lwd = 2, draw = TRUE),
                        borders = list(borders = NA, lwd = 0.5, col = "white"),
                        mask = list(lwd = 0, col = "black"),
                        palette = list(nsteps = 100, n_labels = 20, fn = NULL, digits = 2),
                        plot = TRUE,
                        sphere = list(build = FALSE, keep_image = FALSE)){
                        
                        
  
  n_dat <- length(LIso)
  
  Lrange <- list()
  if (absColBar){
    for (ii in 1:n_dat){
      lt <- LIso[[ii]]
      nv <- length(lt)
      for (jj in 1:nv){
        iso <- lt[[jj]]
        namevar <- names(lt)[[jj]]
        data <- extract(iso$isoscapes[[which]], 1:prod(dim(iso$isoscapes[[which]])))
        Lrange[[namevar]] <- range(data)
      }
    }
  }
  
  y_title_init <- y_title
  sources_init <- sources
  maskCol_init <- maskCol
  palette_init <- palette
  
  for (ii in 1:n_dat){
    lt <- LIso[[ii]]
    nv <- length(lt)
    source_colors <- attr(lt, "col")
    
    if (is.logical(useCSVNames) && useCSVNames){
      period_name <- names(LIso)[[ii]]
      period_name <- gsub(".csv", "", period_name)
    }else{
      period_name <- PeriodNames[ii]
    }
    
    for (jj in 1:nv){
      
      iso <- lt[[jj]]
      namevar <- names(lt)[[jj]]
      
      title <- paste("Isoscape - ", which, " annual", namevar, " - ", period_name)
      
      png(file=paste(title,".png"), width = width, height = height)
      
      if (!y_title$which && is.null(y_title$title)){
        y_title$title <- title
      }
      if (is.null(sources$col)){
        sources$col <- source_colors
      }
      if (!is.null(Lmasks)){

        lmasks <- Lmasks[[ii]]
        nmasks <- length(lmasks)
        maskArg <- list()
        if (length(maskCol) < length(lmasks)){
          maskCol <- 1:length(lmasks)
        }        
        for (kk in 1:nmasks){
          maskArg[[kk]] <- list(mask = lmasks[[kk]], fill = maskCol[[kk]],
                                lwd = mask$lwd, col = mask$col)
        }
      }
      if (is.null(palette$range)){
        if (absColBar){
          palette$range <- Lrange[[jj]]
        }else{
          palette$range <- c(iso$isoscape[[which]]@data@min,
                             iso$isoscape[[which]]@data@max)
        }
        palette$step <- diff(palette$range) / palette$nsteps
      }
      
      
      plot_ISOSCAPE(iso,
                    y_title = y_title,
                    sources = sources,
                    borders = borders,
                    mask    = maskArg,
                    palette = palette,
                    plot = plot,
                    sphere = sphere)
      
      
      dev.off()
      
      y_title <- y_title_init
      sources <- sources_init
      maskCol <- maskCol_init
      palette <- palette_init
      
    }
  }

  
  
}
