#!/usr/bin/env Rscript
##########################################
#### SPATIO-TEMPORAL MODELLING SCRIPT ####
##########################################

###Script to model plant phenology and invertebrate biomass across space and time

setwd("~/Documents/UniversityWork/PhD/ImageAnalysis/Code/")
args = commandArgs(trailingOnly=TRUE)

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Data Manipulation

##Loading data
suppressMessages(library(data.table))
#Pixel data
phenoweight <- fread("../Data/LongTermStudy/PixelData/PhenoWeight.csv", header=T)
#Insect biomass
insectbiomass <- fread("../../GroundSampling/InsectSampling/Data/InsectSampling.csv", header = T)
#Match data
match <- fread("../../GroundSampling/GroundPhenology/Data/VTreeMatching.csv")
#Location data
locations <- fread("../../GroundSampling/PointLocations/locationdistances.csv", header = T)

##Manipulation
#Insect biomass data with additional time variables
insectbiomass$Date <- as.Date(insectbiomass$Date, "%d/%m/%Y")
insectbiomass$DOY <- as.POSIXlt(insectbiomass$Date, format = "%Y/%m/%d")$yday+1
insectbiomass$Year <- as.POSIXlt(insectbiomass$Date, format = "%Y-%m-%d")$year
insectbiomass$Year <- paste0("20", substring(as.character(insectbiomass$Year),2,3))
insectbiomass$conseqDOY <- ifelse(insectbiomass$Year==2020, insectbiomass$DOY + 365, insectbiomass$DOY)
insectbiomass$Week <- strftime(as.POSIXlt(insectbiomass$Date, format = "%Y-%m-%d"),format="%W")
insectbiomass$Week <- as.numeric(as.character(insectbiomass$Week))
insectbiomass$conseqWeek <- ifelse(insectbiomass$Year==2020, insectbiomass$Week + 53, insectbiomass$Week)

##Get geometries and spatial information
suppressMessages(library(sf))
suppressMessages(library(pbapply))
suppressMessages(library(dplyr))

#Invertebrate biomass
insectgrid <- read_sf("../../GroundSampling/InsectSampling/Data/InsectSamplingGrid.shp")
insectgrid <- subset(insectgrid, select = c(name, geometry))
insectlonglats <- do.call(rbind, st_geometry(insectgrid)) %>% as_tibble() %>% setNames(c("lon","lat"))
insectlonglats$name <- insectgrid$name
insectgrid <- right_join(insectgrid, insectlonglats, by ="name")
add_geom_trap <- function(trap){
  singletrap <- subset(insectbiomass, Trap==trap)
  singletrap$geometry <- insectgrid$geometry[insectgrid$name==trap]
  singletrap$lon <- insectgrid$lon[insectgrid$name==trap]
  singletrap$lat <- insectgrid$lat[insectgrid$name==trap]
  return(singletrap)
}
insectgeo <- do.call("rbind", pblapply(unique(insectbiomass$Trap), FUN=add_geom_trap))

#Phenology
treeloc<- read_sf("../Data/LongTermStudy/ITCSegmentation/Manual/ManualDelineation.shp")
colnames(treeloc)[1:2] <- c("VTree", "Tree_Crown_ID")
treeloc <- treeloc[!is.na(treeloc$Tree_Crown_ID) ,]
treeloc$centroid <- st_centroid(treeloc$geometry)
treeloc$geometry <- st_transform(treeloc$geometry, "+proj=longlat +datum=WGS84")
treeloc$centroid <- st_transform(treeloc$centroid, "+proj=longlat +datum=WGS84")
treelonglats <- do.call(rbind, st_geometry(treeloc$centroid)) %>% as_tibble() %>% setNames(c("lon","lat"))
treelonglats$Tree_Crown_ID <- treeloc$Tree_Crown_ID
treeloc <- right_join(treeloc, treelonglats, by="Tree_Crown_ID")
#Add geometry to phenology data
phenogeo <- merge(phenoweight, treeloc[,c(2,5:9)], by = "Tree_Crown_ID")
colnames(phenogeo)[76] <- "TCArea"
#Calculate Percentage Area of Phenology 
phenogeo$PercAreaPheno <- phenogeo$TCArea*(as.numeric(phenogeo$PhenologyPerc)/100)

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Assessing correlation structures

##Temporal Autocorrelation - Rho
# suppressMessages(library(tseries)) #for removing NAs
# #Invertebrate Biomass
# acf(insectbiomass$Dry_Weight, na.action=na.remove, plot=F)$acf[2] #Not very temporally autocorrelated
# acf(insectbiomass[insectbiomass$Trap=="IFT14",]$Dry_Weight, na.action=na.remove, plot=F)$acf[2] #Within trap no temporal autocorrelation either
# #Phenology
# acf(phenogeo$PercAreaPheno, na.action=na.remove, plot=F)$acf[2] #Serially autocorrelated
# acf(phenogeo[phenogeo$Tree_Crown_ID==45,]$PercAreaPheno, lag.max = 380, na.action=na.remove, plot=F)$acf[2] #Within tree crown there is temporal autocorrelation

##Spatial Autocorrelation
# suppressMessages(library(spdep))
# #Invertebrate Biomass 
# k8 <- knn2nb(knearneigh(insectgrid$geometry, k=8, longlat = T), row.names = insectgrid$name)
# dlist <- unlist(nbdists(k8, insectgrid$geometry, longlat = T))
# kd <- dnearneigh(insectgrid$geometry, d1 = 0, d2 = max(dlist)*150, row.names = insectgrid$name)
# W.list <- nb2listw(kd, style = "W")
# W <- nb2mat(kd, style = "B")
# spatmodel <- lm(formula = log(Dry_Weight+1) ~ Trap, data = insectbiomass)
# moran.mc(x = residuals(spatmodel)[1:56], listw = W.list, nsim = 10000) #For one day
# #For all days
# a = 1
# b = 56
# day = 1
# values <- data.table()
# while(b<NROW(insectbiomass)){
#   model <- moran.mc(x = residuals(spatmodel)[a:b], listw = W.list, nsim = 10000)
#   coef <- as.data.table(cbind(model$p.value, model$statistic))
#   values <- rbindlist(list(values, coef))
#   # print(day, model$p.value)
#   day <- day + 1
#   a <- a + 56
#   b <- b + 56
# }
# mean(values$V1)
# mean(values$V2)
# #Bonferoni correction
# padjusted <- p.adjust(values$V1, "bonferroni")
# mean(padjusted)
# iadjusted <- p.adjust(values$V2, "bonferroni")
# mean(iadjusted)
# #No Spatial Autocorrelation

#Phenology
# k8 <- knn2nb(knearneigh(treeloc$centroid, k=8, longlat = T), row.names = treeloc$Tree_Crown_ID)
# dlist <- unlist(nbdists(k8, treeloc$centroid, longlat = T))
# kd <- dnearneigh(treeloc$centroid, d1 = 0, d2 = max(dlist)*150, row.names = treeloc$Tree_Crown_ID)
# W.list <- nb2listw(kd, style = "W")
# W <- nb2mat(kd, style = "B")
# phenoorder <- phenogeo[order(phenogeo$conseqDOY),]
# spatmodel1 <- lm(formula = PercAreaPheno ~ Tree_Crown_ID, data = phenoorder)
# moran.mc(x = residuals(spatmodel1)[1:1447], listw = W.list, nsim = 10000) #For one day 
# #For all days
# a = 1
# b = 1447
# day = 1
# values <- data.frame()
# while(b<NROW(phenoorder)){
#   model <- moran.mc(x = residuals(spatmodel1)[a:b], listw = W.list, nsim = 10000)
#   coef <- as.data.frame(cbind(model$p.value, model$statistic))
#   values <- rbindlist(list(values, coef))
#   # print(day, model$p.value)
#   day <- day + 1
#   a <- a + 1447
#   b <- b + 1447
# }
# mean(values$V1)
# mean(values$V2)
# #Bonferoni correction
# padjusted <- p.adjust(values$V1, "bonferroni")
# mean(padjusted)
# iadjusted <- p.adjust(values$V2, "bonferroni")
# mean(iadjusted)
#No Spatial Autocorrelation

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Area Based Neighbourhood

##Bounding box
#Square
suppressMessages(library(sf))
box <- st_make_grid(treeloc, n=1)
square<- as.data.frame(box[[1]][[1]])
colnames(square) <- c("x","y")
square$group <- 1

#Polygon
pixpoly <- read_sf("../Data/pix4d/Polygon.shp")
pixpoly$geometry <- st_transform(pixpoly$geometry, "+proj=longlat +datum=WGS84")
pixpolyxy <- as.data.frame(st_coordinates(pixpoly))[1:2]

# #Non GGplot
# suppressMessages(library(deldir))
# voro <- deldir(insectgrid$lon, insectgrid$lat)
# plot(insectgrid$lon, insectgrid$lat, type="n", asp=1)
# points(insectgrid$lon, insectgrid$lat, pch=20, col="red", cex=0.5)
# points(treeloc$lon,treeloc$lat,pch=2)
# plot(voro, wlines="tess", wpoints="none", number=F, add=T, lty=1)
# plot(box, add=T)

##Counting point in polygons
suppressMessages(library(spatstat))
suppressMessages(library(maptools))
window = owin(xrange=c(min(square$x), max(square$x)), yrange=c(min(square$y), max(square$y)))
vtess <- dirichlet(ppp(insectgrid$lon, insectgrid$lat, window=window))
treepoints <- ppp(treeloc$lon,treeloc$lat, window=window)
quadratcount(treepoints, tess=vtess)
#In planar
# window <- as(as_Spatial(st_transform(pixpoly$geometry, crs=32650)),"owin")
# insectgridplanar <- as.data.frame(st_coordinates(st_transform(insectgrid$geometry, crs=32650)))
# treelocplanar <- as.data.frame(st_coordinates(st_transform(treeloc$geometry, crs=32650)))
# vtess <- dirichlet(ppp(insectgridplanar$X, insectgridplanar$Y, window=window))
# treepoints <- ppp(treelocplanar$X,treelocplanar$Y, window=window)
# quadratcount(treepoints, tess=vtess)

##Convert objects to sp classes
#Segment to line to polygon
# suppressMessages(library(rgeos))
# boundary <- c(min(square$x), max(square$x), min(square$y), max(square$y))
# voro <- deldir(insectgeo$lon, insectgeo$lat, rw=boundary)
# #Convert data.frame of segment coordinates to a list of SpatialLines objects
# ll <- apply(voro$dirsgs, 1, FUN=function(X) {readWKT(sprintf("LINESTRING(%s %s, %s %s)", X[1], X[2], X[3], X[4]))})
# #Convert SpatialLines list to SpatialPolygons object
# vpolys <- gPolygonize(ll)

#Tesselation to polygon
suppressMessages(library(spatstat))
suppressMessages(library(maptools)) #For the spatial conversion
# window <- owin(xrange=c(min(square$x), max(square$x)), yrange=c(min(square$y), max(square$y)))
window <- as(as_Spatial(st_transform(pixpoly$geometry, crs=32650)),"owin")
insectgridplanar <- as.data.frame(st_coordinates(st_transform(insectgrid$geometry, crs=32650)))
vtess <- dirichlet(ppp(insectgridplanar$X, insectgridplanar$Y, window=window))
tess2SPdf <- function(x) { 
  stopifnot(is.tess(x)) 
  require(sp) 
  require(spatstat) 
  require(spatstat.utils)
  a <- tiles(x) 
  tess.labels <- names(a) 
  c <- list() 
  for(i in seq(a)){ 
    
    b <- as.polygonal(a[[i]]) 
    closering <- function(df) { 
      df[c(seq(nrow(df)), 1), ] 
    } 
    pieces <- lapply(b$bdry, 
                     function(p) { 
                       Polygon(coords=closering(cbind(p$x,p$y)), 
                               hole=is.hole.xypolygon(p)) 
                     }) 
    c[[i]] <- Polygons(pieces, tess.labels[i]) 
    
  } 
  
  d <- data.frame() 
  d <- d[seq(x$n),] 
  row.names(d) <- names(x$tiles) 
  return(SpatialPolygonsDataFrame(SpatialPolygons(c), data=d)) 
} 
vpolys <- tess2SPdf(vtess)
#Add planar CRS and convert to lon/lat
proj4string(vpolys) <- CRS("+init=epsg:32650")
vpolys <- spTransform(vpolys, CRS("+proj=longlat +datum=WGS84"))

##Which points are within polygons
#Create sf object of voronoi polygons
sf_vpolys <- st_set_crs(st_as_sf(vpolys, crs = st_crs(treeloc)), st_crs(treeloc))
#Add row names
sf_vpolys$names <- row.names(sf_vpolys)
#Get tree locations
treeregions <- treeloc[,c(2,8:9)]
#Transform from lat/lon to planar
sf_vpolys_planar <- st_transform(sf_vpolys, 32650)
treeregions_planar <- st_transform(treeregions, 32650)
#Find intersection and extract polygons name - sparse creates logical matrix
treeregions$region <- apply(st_intersects(sf_vpolys_planar, treeregions_planar$geometry, sparse = FALSE), 2, 
                        function(col) { sf_vpolys_planar[which(col),]$names})

##Calculate Area of Tree Crowns in each polygon
TotalTCArea <- data.frame(VPoly = 1:56,
                          TotalTCArea = pbsapply(1:56, function(x){sum(st_area(st_intersection(sf_vpolys_planar[x,], treeregions_planar$geometry)))}))

##Calculate Area of Polygons
PolyArea <- as.data.frame(do.call('rbind', lapply(sf_vpolys_planar[[1]], st_area)))
colnames(PolyArea) <- c("PolyArea")
PolyArea$VPoly <- row.names(PolyArea)

##Add region/neighbourhood to phenodata
phenogeonb <- as.data.frame(phenogeo)
add_neighbourhood <- function(tree){
  singletree <- subset(phenogeo, Tree_Crown_ID==tree)
  if(length(treeregions$region[treeregions$Tree_Crown_ID==tree][[1]])>1){
    singletree$VPoly <- list(treeregions$region[treeregions$Tree_Crown_ID==tree][[1]])
    singletree$MultiPoly <- length(treeregions$region[treeregions$Tree_Crown_ID==tree][[1]])
  } else {
    singletree$VPoly <- treeregions$region[treeregions$Tree_Crown_ID==tree][[1]]
    singletree$MultiPoly <- 1
  }
  return(singletree)
}
phenogeonb <- as.data.frame(rbindlist(pblapply(unique(phenogeonb$Tree_Crown_ID), FUN=add_neighbourhood)))

##For MultiPoly trees assign to polygon where most of tree resides
#Transform from lat/lon to planar
sf_vpolys_planar <- st_transform(sf_vpolys, 32650)
treeregions_planar <- st_transform(treeregions, 32650)
treeregions_planar$treearea <- st_area(treeregions_planar)
#Find which polygon each tree crown is in
assign_multipolys <- function(tree){
  
  #Subset data frames
  subcross <- subset(treeregions_planar, Tree_Crown_ID==tree)
  singletree <- subset(phenogeonb, Tree_Crown_ID==tree)
  
  #Find trees that cross regions
  intpolys <- sf_vpolys_planar$names[row.names(sf_vpolys_planar) %in% which(st_intersects(sf_vpolys_planar$geometry, subcross$geometry, sparse = F))]
  
  #Find intersection area
  intarea <- st_area(st_intersection(sf_vpolys_planar$geometry, subcross$geometry))
  tcarea <- subcross$treearea
  percarea <- (as.numeric(intarea/tcarea))*100
 
  #Select polygon with most area covered
  singletree$MainPoly <- intpolys[which(percarea==max(percarea))]
  singletree$PercPoly <- max(percarea)
  
  #Recombine
  return(singletree)
}
phenogeomp <- as.data.frame(rbindlist(pblapply(unique(phenogeonb$Tree_Crown_ID), FUN=assign_multipolys)))
# phenogeomplite <- dplyr::select(phenogeomp,
#                                 conseqWeek, conseqDOY, DOY, Tree_Crown_ID, 
#                                 Phenology, CPWeighting, BinPheno, PercAreaPheno)
# write.csv(phenogeomplite, "../Data/LongTermStudy/PixelData/PhenoGeoMPLite.csv")

#Find which tree crowns are in poly 1
# poly = 1
# inttc <- crosstreeregions$Tree_Crown_ID[row.names(crosstreeregions) %in% which(st_intersects(crosstreeregions$geometry, sf_vpolys$geometry[poly], sparse = F))]
# intarea <- st_area(st_intersection(crosstreeregions$geometry, sf_vpolys$geometry[1]))
# intdf <- as.data.frame(cbind(poly, inttc, intarea))
# intdf$tcarea <- crosstreeregions$treearea[crosstreeregions$Tree_Crown_ID %in% intdf$inttc]
# intdf$percarea <- (intdf$intarea/as.vector(intdf$tcarea))*100
# intdf$newpoly <- ifelse(intdf$percarea >= 50.0, intdf$poly, NA)

##Summarise Neighbourhood - maybe replace NA with None
# onepoly <- subset(phenogeonb, MultiPoly==1)
# onepoly$VPoly <- as.factor(as.numeric(onepoly$VPoly))
suppressMessages(library(dplyr)); suppressMessages(library(tidyr))
phenompsum <- phenogeomp %>%
  dplyr::select(conseqWeek, conseqDOY, DOY, Tree_Crown_ID, 
                Phenology, CPWeighting, BinPheno,
                TreeFamily, TreeGenus, TreeSpecies, TreeTaxon,
                geometry, centroid, lon, lat,
                TCArea, PercAreaPheno, VPoly, MultiPoly, MainPoly, PercPoly) %>%
  tidyr::replace_na(list(PercAreaPheno = 0)) %>%
  dplyr::group_by(conseqWeek, MainPoly, Phenology) %>%
  dplyr::summarise(PercAreaPheno = sum(PercAreaPheno)) %>%
  mutate(TotalAreaPheno = sum(PercAreaPheno))

##Add Neighbourhood PhenoStatus
add_pheno_status <- function(trap){
  
  #Convert Trap to number
  singletrap <- subset(insectgeo, Trap==trap)
  IFT <- as.numeric(strsplit(as.character(trap), "IFT")[[1]][2])
  singletrap$IFT <- IFT
  
  #For each week extract pheno characteristics
  each_week <- function(week){
    
    #Extract phenology area
    totalphenoarea <- unique(phenompsum$TotalAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week])
    lfphenoarea <- phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf flush"][which(!is.na(phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf flush"]))]
    flphenoarea <- phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Flowering"][which(!is.na(phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Flowering"]))]
    frphenoarea <- phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Fruiting"][which(!is.na(phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Fruiting"]))]
    lsphenoarea <- phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf senescence"][which(!is.na(phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf senescence"]))]
    llphenoarea <- phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf loss"][which(!is.na(phenompsum$PercAreaPheno[phenompsum$MainPoly==IFT & phenompsum$conseqWeek==week & phenompsum$Phenology=="Leaf loss"]))]
    
    #Adding missing values
    totalphenoarea <- ifelse(length(totalphenoarea)==0, 0, totalphenoarea)
    lfphenoarea <- ifelse(length(lfphenoarea)==0, 0, lfphenoarea)
    flphenoarea <- ifelse(length(flphenoarea)==0, 0, flphenoarea)
    frphenoarea <- ifelse(length(frphenoarea)==0, 0, frphenoarea)
    lsphenoarea <- ifelse(length(lsphenoarea)==0, 0, lsphenoarea)
    llphenoarea <- ifelse(length(llphenoarea)==0, 0, llphenoarea)
    
    #Recombine
    comb <- do.call("cbind", list(totalphenoarea, lfphenoarea, flphenoarea, 
                                  frphenoarea, lsphenoarea, llphenoarea))
    
  }
  PhenoArea <- do.call("rbind", lapply(unique(singletrap$conseqWeek), FUN=each_week))
  colnames(PhenoArea) <- c("TotalPhenoArea", "FlushArea", "FlowerArea", 
                           "FruitArea", "SenesArea", "LossArea")
  singletrap <- do.call("cbind", list(singletrap, PhenoArea))
  
  #Add TotalTCArea & PolyArea
  singletrap$TotalTCArea <- TotalTCArea$TotalTCArea[TotalTCArea$VPoly==IFT]
  singletrap$PolyArea <- PolyArea$PolyArea[PolyArea$VPoly==IFT]
  
  #Calculate PhenoStatus per Polygon
  singletrap$TotalPhenoStatus <- (singletrap$TotalPhenoArea/singletrap$PolyArea)*100
  singletrap$FlushStatus <- (singletrap$FlushArea/singletrap$PolyArea)*100
  singletrap$FlowerStatus <- (singletrap$FlowerArea/singletrap$PolyArea)*100
  singletrap$FruitStatus <- (singletrap$FruitArea/singletrap$PolyArea)*100
  singletrap$SenesStatus <- (singletrap$SenesArea/singletrap$PolyArea)*100
  singletrap$LossStatus <- (singletrap$LossArea/singletrap$PolyArea)*100
  
  #Calculate PhenoStatus across tree crowns
  singletrap$TCPhenoStatus <- (singletrap$TotalPhenoArea/singletrap$TotalTCArea)*100
  singletrap$TCFlushStatus <- (singletrap$FlushArea/singletrap$TotalTCArea)*100
  singletrap$TCFlowerStatus <- (singletrap$FlowerArea/singletrap$TotalTCArea)*100
  singletrap$TCFruitStatus <- (singletrap$FruitArea/singletrap$TotalTCArea)*100
  singletrap$TCSenesStatus <- (singletrap$SenesArea/singletrap$TotalTCArea)*100
  singletrap$TCLossStatus <- (singletrap$LossArea/singletrap$TotalTCArea)*100
  
  return(singletrap)
}
insectphenoarea <- as.data.frame(rbindlist(pblapply(unique(insectgeo$Trap), FUN=add_pheno_status)))

###--------------------------------------------------------------------------------###

###Mixed Model Analysis

##All phenology
suppressMessages(library(lme4)); suppressMessages(library(lmerTest))
insectphenoarea$scaleconseqWeek <- scale(insectphenoarea$conseqWeek)
summary(lmmmodel1 <- lmer(log(Dry_Weight+1) ~ TotalPhenoStatus * scaleconseqWeek +
                            #Effects of Week on DryWeight is different for different Traps and Tree crowns
                            (scaleconseqWeek | Trap), 
                          data = insectphenoarea, na.action = na.omit, REML=T))
summary(lmmmodel2 <- lmer(log(Dry_Weight+1) ~ TotalPhenoStatus * scaleconseqWeek +
                            #Effects of Week on DryWeight is different for different Traps and Tree crowns
                            (1| Trap), 
                          data = insectphenoarea, na.action = na.omit, REML=T))
anova(lmmmodel1, lmmmodel2)
anova(lmmmodel2)
acf(residuals(lmmmodel2))$acf[2]

###--------------------------------------------------------------------------------###

##Plot points
suppressMessages(library(ggplot2))
suppressMessages(library(ggvoronoi))
sfphenogeo <- st_as_sf(phenogeo)
p1 <- ggplot(insectgrid) +
  geom_sf(data = subset(sfphenogeo, conseqWeek==21), aes(fill=PercAreaPheno), lwd = 0) +
  geom_point(col = "blue", pch = 2, aes(x = lon, y = lat)) +
  geom_sf(data = sf_vpolys, alpha = 0.01) +
  # geom_path(stat = "voronoi", alpha=0.5, size=0.25, aes(x = lon, y = lat)) +
  # geom_path(data = square, aes(x, y, group=group)) +
  # geom_path(data = pixpolyxy,aes(x=X, y=Y)) +
  scale_fill_gradient(low="lightgreen", high="darkgreen", name="Area of Phenology") + 
  # scale_x_continuous(expand = c(0,0), limits=c(0,8), breaks = seq(min(0), max(8), by = 1), labels=seq(0,480,60)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
# png("../Results/SpatioTemporal/VoronoiPhenoArea.png", width=1800, height=1800, res=200)
# p1
# dev.off()

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Distance Based Neighbourhood Phenology

# #Matrix structure: Col, Row
# distmatrix <- st_distance(insectgrid, treeloc)
# max(distmatrix)

##Create Distance Weighting
segment=5
distrange <- seq(segment, ceiling(520/segment)*segment, segment)
#Support choice of bandwith and kernel based on ecology
#Inverse distance squared weighting
WeightDistSq <- 1/(distrange^2)
#Epanechnikov Kernel
WeightEpan <- 0.75*(1-(distrange^2))
#Weibull distribution
WeightWeib <- pweibull(distrange, shape=1)
#Kernel Bandwidth
#Weight of observations is zero beyond a certain bandwidth: #Bandwidth = radius
h = 5
#Uniform Kernel
WeightUni <- ifelse(distrange < h, 1, 0)
#BiSquare Kernel
WeightBiSq <- ifelse(distrange < h, ( 1 - (distrange/h)^2 )^2, 0)
#Gaussian Kernel
WeightGauss <- exp(-0.5*((distrange/h)^2))

##Data manipulation
#Transform to planar
insectgrid_planar <- st_transform(st_as_sf(insectgrid), crs=32650)
#Create smaller dataset 
phenogeomplite <- dplyr::select(phenogeomp,
                                conseqWeek, conseqDOY, DOY, Tree_Crown_ID, 
                                Phenology, CPWeighting, BinPheno, PercAreaPheno
                                # TreeFamily, TreeGenus, TreeSpecies, TreeTaxon,
                                # geometry, centroid, lon, lat,
                                # TCArea, VPoly, MultiPoly, MainPoly, PercPoly
                                )
#Convert to data table for speed
suppressMessages(library(data.table))
phenogeomplite <- as.data.table(phenogeomplite)

##Add distance based phenology area
suppressMessages(library(plyr))
dist_pheno_area <- function(week){
  
  ##For a specific doy
  singleweek <- subset(phenogeomplite, conseqWeek==week)
  
  ##For a specific trap
  each_trap <- function(trap){
    
    #Subset trap
    singletrap <- insectgrid_planar[insectgrid_planar$name==trap,]
    
    ##Find which trees are within X distance of trap - 5 m intervals
    find_trees <- function(dist){

      #Trees within
      treeswithin <- treeregions_planar$Tree_Crown_ID[which(st_is_within_distance(singletrap, treeregions_planar, dist, sparse=F))]
      numtrees <- length(treeswithin)
      percarea <- singleweek$PercAreaPheno[singleweek$Tree_Crown_ID %in% treeswithin]
      phen <- singleweek$Phenology[singleweek$Tree_Crown_ID %in% treeswithin]
      
      #Calculate number of phenological events
      phenonum <- length(which(!is.na(percarea)==T))
      phenos <- plyr::count(phen)
      flushnum <- phenos$freq[phenos$x == "Leaf flush"][which(!is.na(phenos$freq[phenos$x == "Leaf flush"]))]
      flownum <- phenos$freq[phenos$x == "Flowering"][which(!is.na(phenos$freq[phenos$x == "Flowering"]))]
      fruitnum <- phenos$freq[phenos$x == "Fruiting"][which(!is.na(phenos$freq[phenos$x == "Fruiting"]))]
      senesnum <- phenos$freq[phenos$x == "Leaf senescence"][which(!is.na(phenos$freq[phenos$x == "Leaf senescence"]))]
      lossnum <- phenos$freq[phenos$x == "Leaf loss"][which(!is.na(phenos$freq[phenos$x == "Leaf loss"]))]
      # deadnum <- phenos$freq[phenos$x == "Dead"]
      
      #Change empty to NA
      flushnum <- ifelse(length(flushnum)==0, 0, flushnum)
      flownum <- ifelse(length(flownum)==0, 0, flownum)
      fruitnum <- ifelse(length(fruitnum)==0, 0, fruitnum)
      senesnum <- ifelse(length(senesnum)==0, 0, senesnum)
      lossnum <- ifelse(length(lossnum)==0, 0, lossnum)
      
      #Sum Area Phenology and split by phenology
      sumphenoarea <- sum(percarea, na.rm=T)
      flushsum <- sum(percarea[which(phen=="Leaf flush")])
      flowsum <- sum(percarea[which(phen=="Flowering")])
      fruitsum <- sum(percarea[which(phen=="Fruiting")]) 
      senessum <- sum(percarea[which(phen=="Leaf senescence")])
      losssum <- sum(percarea[which(phen=="Leaf loss")])
      
      #Match distance to weights - 5m slots
      distw <- as.numeric(sapply(as.numeric(5*ceiling(st_distance(singletrap, treeregions_planar[treeswithin,])/5)),
                                 function(x){WeightDistSq[which(distrange==x)]}))
      
      #Distance Weighted Sum
      wsumphenoarea <- sum(percarea * distw, na.rm=T)
      flushwsum <- sum(percarea[which(phen=="Leaf flush")] * distw[which(phen=="Leaf flush")])
      flowwsum <- sum(percarea[which(phen=="Flowering")] * distw[which(phen=="Flowering")])
      fruitwsum <- sum(percarea[which(phen=="Fruiting")] * distw[which(phen=="Fruiting")])
      seneswsum <- sum(percarea[which(phen=="Leaf senescence")] * distw[which(phen=="Leaf senescence")])
      losswsum <- sum(percarea[which(phen=="Leaf loss")] * distw[which(phen=="Leaf loss")])
      
      #Calculate percentage of pheno area within distance
      # bufferarea <- as.numeric(st_area(st_buffer(singletrap, dist)))
      # percsumphenoarea <- (as.numeric(sumphenoarea/bufferarea))*100
      # percwsumphenoarea <- (as.numeric(wsumphenoarea/bufferarea))*100
      # 
      # percsumflusharea <- (as.numeric(flushsum/bufferarea))*100
      # percwsumflusharea <- (as.numeric(flushwsum/bufferarea))*100
      # percsumflowarea <- (as.numeric(flowsum/bufferarea))*100
      # percwsumflowarea <- (as.numeric(flowwsum/bufferarea))*100
      # percsumfruitarea <- (as.numeric(fruitsum/bufferarea))*100
      # percwsumfruitarea <- (as.numeric(fruitwsum/bufferarea))*100
      # percsumsenesarea <- (as.numeric(senessum/bufferarea))*100
      # percwsumsenesarea <- (as.numeric(seneswsum/bufferarea))*100
      # percsumlossarea <- (as.numeric(losssum/bufferarea))*100
      # percwsumlossarea <- (as.numeric(losswsum/bufferarea))*100
      
      #Combine
      comb <- as.data.table(do.call("cbind", list(week, trap, dist, numtrees,
                                                  phenonum, flushnum, flownum, fruitnum, senesnum, lossnum,
                                                  sumphenoarea, flushsum, flowsum, fruitsum, senessum, losssum,
                                                  wsumphenoarea, flushwsum, flowwsum, fruitwsum, seneswsum, losswsum
                                                  # ,bufferarea, percsumphenoarea, percwsumphenoarea,
                                                  # percsumflusharea, percwsumflusharea,
                                                  # percsumflowarea, percwsumflowarea,
                                                  # percsumfruitarea, percwsumfruitarea,
                                                  # percsumsenesarea, percwsumsenesarea,
                                                  # percsumlossarea, percwsumlossarea
                                                  )))
      return(comb)
    }
    eachdist <- as.data.table(rbindlist(lapply(distrange, FUN=find_trees)))
    return(eachdist)
  }
  alltraps <- as.data.table(rbindlist(lapply(unique(insectgrid_planar$name), FUN=each_trap)))
  return(alltraps)
}
suppressMessages(library(pbapply))
distphenoarea <- as.data.table(rbindlist(pblapply(unique(phenogeonb$conseqWeek), FUN=dist_pheno_area)))
colnames(distphenoarea) <- c("conseqWeek", "Trap", "Distance", "NumTrees",
                             "PhenoNum", "FlushNum", "FlowNum", "FruitNum", "SenesNum", "LossNum",
                             "SumPhenoArea", "FlushSum", "FlowSum", "FruitSum", "SenesSum", "LossSum",
                             "WeightSumPhenoArea", "FlushWSum", "FlowWSum", "FruitWSum", "SenesWSum", "LossWSum"
                             # ,"BufferArea", "PercSumPhenoArea", "PercWeightSumPhenoArea",
                             # "PercSumFlushArea", "PercWSumFlushArea",
                             # "PercSumFlowArea", "PercWSumFlowArea",
                             # "PercSumFruitArea", "PercWSumFruitArea",
                             # "PercSumSenesArea", "PercWSumSenesArea",
                             # "PercSumLossArea", "PercWSumLossArea"
                             )
# fwrite(distphenoarea, paste0("../Data/LongTermStudy/PixelData/DistPhenoArea-",segment,".csv"), row.names=F)
# DPA <- distphenoarea
suppressMessages(library(data.table))
DPA <- fread("../Data/LongTermStudy/PixelData/DistPhenoArea.csv") #Faster load

##Add Dry Weight and IFT column
colnames(DPA)[1] <- c("conseqWeek")
insectdw <- subset(insectbiomass, select=c(conseqWeek, Trap, Dry_Weight))
DPAdw <- merge(DPA, insectdw, by=c("conseqWeek", "Trap"), all.x=T)
DPAdw$IFT <- sapply(DPA$Trap, function(x){as.numeric(strsplit(x, "IFT")[[1]][2])})

##Calculate updated area of buffer geom
suppressMessages(library(sf))
allbuffergeoms <- st_set_crs(st_sfc(sapply(distrange,function(x){st_buffer(st_transform(insectgrid$geometry, crs=32650),x)})),32650)
allintbuffgeoms <- st_sfc(st_intersection(allbuffergeoms, st_transform(pixpoly, crs=32650)))
newbufferarea <- as.numeric(st_area(allintbuffgeoms))
updatebufferarea <- data.frame(NewBufferArea = newbufferarea,
                               Distance = unlist(lapply(distrange, function(x){rep(x,56)})),
                               IFT = rep(1:56,length(distrange)))
DPAdw <- merge(DPAdw, updatebufferarea, by=c("Distance", "IFT"))
#Update percentages
# DPAdw$BufferArea <- NULL
DPAdw$PercSumPhenoArea <- (as.numeric(DPAdw$SumPhenoArea/DPAdw$NewBufferArea))*100
DPAdw$PercWeightSumPhenoArea <- (as.numeric(DPAdw$WeightSumPhenoArea/DPAdw$NewBufferArea))*100
DPAdw$PercSumFlushArea <- (as.numeric(DPAdw$FlushSum/DPAdw$NewBufferArea))*100
DPAdw$PercWSumFlushArea <- (as.numeric(DPAdw$FlushWSum/DPAdw$NewBufferArea))*100
DPAdw$PercSumFlowArea <- (as.numeric(DPAdw$FlowSum/DPAdw$NewBufferArea))*100
DPAdw$PercWSumFlowArea <- (as.numeric(DPAdw$FlowWSum/DPAdw$NewBufferArea))*100
DPAdw$PercSumFruitArea <- (as.numeric(DPAdw$FruitSum/DPAdw$NewBufferArea))*100
DPAdw$PercWSumFruitArea <- (as.numeric(DPAdw$FruitWSum/DPAdw$NewBufferArea))*100
DPAdw$PercSumSenesArea <- (as.numeric(DPAdw$SenesSum/DPAdw$NewBufferArea))*100
DPAdw$PercWSumSenesArea <- (as.numeric(DPAdw$SenesWSum/DPAdw$NewBufferArea))*100
DPAdw$PercSumLossArea <- (as.numeric(DPAdw$LossSum/DPAdw$NewBufferArea))*100
DPAdw$PercWSumLossArea <- (as.numeric(DPAdw$LossWSum/DPAdw$NewBufferArea))*100

##Calculate TC area in each buffer geom
TotalBufferTCArea <- data.table(Distance = unlist(pblapply(distrange, function(x){rep(x,56)})),
                                IFT = rep(1:56,length(distrange)),
                                TotalBufferTCArea = pbsapply(1:length(distrange)*56, 
                                                             function(x){sum(st_area(st_intersection(allbuffergeoms[x], treeregions_planar$geometry)))}))
DPAdw <- merge(DPAdw, TotalBufferTCArea, by=c("Distance", "IFT"))
#Caluclate Percentage Areas for TC
DPAdw$TCPercSumPhenoArea <- (as.numeric(DPAdw$SumPhenoArea/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWeightSumPhenoArea <- (as.numeric(DPAdw$WeightSumPhenoArea/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercSumFlushArea <- (as.numeric(DPAdw$FlushSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWSumFlushArea <- (as.numeric(DPAdw$FlushWSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercSumFlowArea <- (as.numeric(DPAdw$FlowSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWSumFlowArea <- (as.numeric(DPAdw$FlowWSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercSumFruitArea <- (as.numeric(DPAdw$FruitSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWSumFruitArea <- (as.numeric(DPAdw$FruitWSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercSumSenesArea <- (as.numeric(DPAdw$SenesSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWSumSenesArea <- (as.numeric(DPAdw$SenesWSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercSumLossArea <- (as.numeric(DPAdw$LossSum/DPAdw$TotalBufferTCArea))*100
DPAdw$TCPercWSumLossArea <- (as.numeric(DPAdw$LossWSum/DPAdw$TotalBufferTCArea))*100

###--------------------------------------------------------------------------------###

####Model Analysis

###Exploring the data
complete <- DPAdw[complete.cases(DPAdw$Dry_Weight),]
summary(explmodel1 <- lm(log(Dry_Weight+1) ~ WeightSumPhenoArea * conseqWeek, data = complete))
par(mfrow=c(3,2))
plot(explmodel1, which = c(1), col = 1, add.smooth = FALSE, caption = "")
#Clear violation of heterogeneity
hist(resid(explmodel1), xlab = "Residuals", main = "")
qqnorm(resid(explmodel1)); qqline(resid(explmodel1))
#Not normally distributed
plot(as.factor(complete$conseqWeek), resid(explmodel1), xlab = "WOY", ylab = "Residuals")
plot(complete$WeightSumPhenoArea, residuals(explmodel1), xlab = "WeightSumPhenoArea", ylab = "Residuals")
#Smaller the WeightSumPhenoArea the larger the variance

###Generalised Least Squares Model
suppressMessages(library(nlme))
#Deals with heterogeneity (heteroscedacity)
#No nested structure
#Not dealt with non normal data
summary(glsmodel1 <- gls(log(Dry_Weight+1) ~ WeightSumPhenoArea * conseqWeek,
                         data = complete, na.action = na.omit))
#GLS with variance weighting
#Variance covariate = WeightSumPhenoArea
#Variance covariate is not nominal and includes 0 values so use varExp or varConstPower
vf <- varExp(form = ~ WeightSumPhenoArea) #Models the heteroscedacity in WSPA
# vf <- varComb(varExp(form = ~ WeightSumPhenoArea), varPower(form = ~ conseqWeek))
summary(glsmodel2 <- gls(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         data = complete, na.action = na.omit,
                         weights = vf, method = "REML"))
anova(glsmodel1, glsmodel2) #Model 2 has lower AIC

###Linear Mixed Effects Models
#Include nested random effects
#Includes compensation of heteroscedacity
##Compare inclusion of random effects
summary(glsmodel3 <- gls(log(Dry_Weight+1) ~ WeightSumPhenoArea * conseqWeek,
                         data = complete, na.action=na.omit, method = "REML"
                         ))
summary(lmmmodel1 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea * conseqWeek,
                         random = ~ 1 | Trap,
                         data = complete, na.action=na.omit, method = "REML"
                         ))
summary(lmmmodel2 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea * conseqWeek,
                         random = ~ conseqWeek | Trap,
                         data = complete, na.action=na.omit, method = "REML"
                         ))
anova(glsmodel3, lmmmodel1, lmmmodel2) #Model 2 lower AIC

##Drop values that are non significant
summary(lmmmodel3 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         random = ~ conseqWeek | Trap,
                         data = complete, na.action=na.omit, method = "REML"
                         ))
##Add variance structure for heteroscedacity
summary(lmmmodel4 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         random = ~ conseqWeek | Trap,
                         weights = varExp(form = ~ WeightSumPhenoArea),
                         data = complete, na.action=na.omit, method = "REML"
                         ))
##Centre the conseqWeek value
complete$scaleconseqWeek <- scale(complete$conseqWeek, center = T, scale = F)
summary(lmmmodel4 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         random = ~ scaleconseqWeek | Trap,
                         weights = varExp(form = ~ WeightSumPhenoArea),
                         data = complete, na.action=na.omit, method = "REML"
                         ))

##Checking for linearity
# library(lattice)
# p1 <- xyplot(residuals(lmmmodel4) ~ WeightSumPhenoArea | scaleconseqWeek,
#        data = complete, ylab = "Residuals", xlab = "WOY", panel = function(x,y){
#          panel.grid(h = -1, v = 2)
#          panel.points(x, y, col = 1)
#          panel.loess(x, y, span = 0.5, col = 1,lwd=2)})
# png("../Results/Sandbox/lattice.png", width=2400, height=2400, res=140)
# p1
# dev.off()
#They are linear so a additive model is not necessary

###Correlation Structures
##Reference structure
summary(glsmodel4 <- gls(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         data = complete, na.action = na.omit, method="REML"))
#Test for autocorrelation
acf(residuals(glsmodel4)) #Is autocorrelation

##Correlations
#Need to add grouping factor
complete$ID <- row.names(complete)
#corCompSymm - does not make ecological sense
#corAR1 - further away the residuals are in time the less their correlation is
#corARMA - try different p and q between 0 & 3 - p=2,q=1:crashes R
#Esimating parameters
suppressMessages(library(tseries))
coeff <- arma(sqrt(complete$Dry_Weight), order = c(1,0))$coef[[1]]
#Apply parameters
summary(glsmodel5 <- gls(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         correlation = corARMA(c(coeff), form = ~ scaleconseqWeek | Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
acf(residuals(glsmodel5))$acf[2]

##Include random effects
summary(lmmmodel5 <- lme(log(Dry_Weight+1) ~ WeightSumPhenoArea,
                         random = ~scaleconseqWeek|Trap,
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
acf(residuals(lmmmodel5))$acf[2]

##Inlcude variance structure
summary(lmmmodel6 <- lme(log(Dry_Weight+1) ~ log(WeightSumPhenoArea+1),
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ WeightSumPhenoArea),
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
acf(resid(lmmmodel6, type = "normalized"))
qqnorm(resid(lmmmodel6, type = "normalized"))
qqline(resid(lmmmodel6, type = "normalized"))

##Allow for non-normal distribution and overdispersion
coeff <- arma(complete$scaleconseqWeek, order = c(1,0))$coef[[1]]
library(moments)
skewness(log(complete$Dry_Weight+1), na.rm = TRUE)
skewness(sqrt(complete$Dry_Weight), na.rm = TRUE)
skewness(sign(complete$Dry_Weight) * abs(complete$Dry_Weight)^(1/3), na.rm = TRUE)
complete$cubDW <- sign(complete$Dry_Weight) * abs(complete$Dry_Weight)^(1/3)
complete$cubWSPA <- sign(complete$WeightSumPhenoArea) * abs(complete$WeightSumPhenoArea)^(1/3)
boxcox.
summary(lmmmodel7 <- lme(cubDW ~ cubWSPA,
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ cubWSPA),
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
qqnorm(resid(lmmmodel7, type = "normalized"))
qqline(resid(lmmmodel7, type = "normalized"))

###Aggregation of Data
test <- complete %>%
  dplyr::group_by(conseqWeek, Trap) %>%
  dplyr::summarise(MeanPhenoArea=mean(WeightSumPhenoArea))
test2 <- merge(test, insectdw, by=c("conseqWeek", "Trap"), all.x=T)
test2$scaleconseqWeek <- scale(test2$conseqWeek, center = T, scale = F)
test2$ID <- row.names(test2)
summary(lmmmodel8 <- lme(log(Dry_Weight+1) ~ log(MeanPhenoArea+1),
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ log(MeanPhenoArea+1)),
                         # correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = test2, na.action = na.omit, method="REML", control = lmeControl(opt='optim')
                         ))
hist(resid(lmmmodel8, type = "normalized"))
plot(test2$MeanPhenoArea, resid(lmmmodel8, type = "normalized"))

##Individual phenology
suppressMessages(library(broom.mixed))
summary(lmmflush <- lme(log(Dry_Weight+1) ~ log(FlushWSum+1),
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ log(FlushWSum+1)),
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
summary(lmmflow <- lme(log(Dry_Weight+1) ~ log(FlowWSum+1),
                        random = ~scaleconseqWeek|Trap,
                        weights = varExp(form = ~ log(FlowWSum+1)),
                        correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                        data = complete, na.action = na.omit, method="REML"
                        ))
summary(lmmfruit <- lme(log(Dry_Weight+1) ~ log(FruitWSum+1),
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ log(FruitWSum+1)),
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
summary(lmmsenes <- lme(log(Dry_Weight+1) ~ log(SenesWSum+1),
                         random = ~scaleconseqWeek|Trap,
                         weights = varExp(form = ~ log(SenesWSum+1)),
                         correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                         data = complete, na.action = na.omit, method="REML"
                         ))
summary(lmmloss <- lme(log(Dry_Weight+1) ~ log(LossWSum+1),
                       random = ~scaleconseqWeek|Trap,
                       weights = varExp(form = ~ log(LossWSum+1)),
                       correlation = corARMA(value=c(coeff), form = ~ scaleconseqWeek|Trap/ID, p = 1, q = 0),
                       data = complete, na.action = na.omit, method="REML"
))

suppressMessages(library(broom.mixed))
mmphenotable <- rbindlist(list(tidy(lmmflush), tidy(lmmflow), tidy(lmmfruit), tidy(lmmsenes), tidy(lmmloss)))
colnames(mmphenotable) <- c("Effect", "Group", "Variable", "Estimate", "SE", "t", "DF", "P")
# write.csv(mmphenotable, "../Results/SpatioTemporal/PhenoMEMTable.csv")
# summary(lmmmodel4 <- lmer(log(Dry_Weight+1) ~ FlushWSum * FlowWSum * FruitWSum * SenesWSum * LossWSum * scaleconseqWeek +
#                           #Effects of Week on DryWeight is different for different Traps
#                           (scaleconseqWeek | Trap), 
#                         data = DPAdw, na.action = na.omit, REML=T))
# anova(lmmmodel4)

###--------------------------------------------------------------------------------###

###Plotting

##Create individual buffer polygons for trap 1
suppressMessages(library(sf))

allbuffergeoms <- st_set_crs(st_sfc(sapply(distrange,function(x){st_buffer(st_transform(insectgrid$geometry, crs=32650),x)})),32650)
newbuffergeoms <- st_sfc(st_intersection(allbuffergeoms, st_transform(pixpoly, crs=32650)))
# updatebuffer <- data.frame(Distance = unlist(lapply(distrange, function(x){rep(x,56)})),
#                            IFT = rep(1:56,length(distrange)))
distrange <- seq(5,515,5)
trapsubdist <- data.frame(Distance = distrange)
each_trap <- function(trap){
  
  ##Extract buffer geometries for single trap
  trapbuffergeoms <- st_set_crs(st_sfc(sapply(distrange,function(x){st_buffer(st_transform(insectgrid$geometry[trap], crs=32650),x)})),32650)
  
  each_week <- function(week){
    trapsub <- subset(DPAdw, IFT==trap & conseqWeek == week)
    weekcomb <- merge(trapsub, trapsubdist, by="Distance")
    x=2
    y=1
    segmentgeoms <- trapbuffergeoms[1]
    while(x<=NROW(trapsubdist) & y<=(NROW(trapsubdist)-1)){
      seggeom <- st_difference(trapbuffergeoms[x], trapbuffergeoms[y])
      segmentgeoms <- do.call("rbind", list(segmentgeoms, seggeom))
      x <- x+1
      y <- y+1
    }
    st_geometry(weekcomb) <- st_sfc(segmentgeoms)
    return(weekcomb)
  }
  weekcombined <- as.data.table(rbindlist(lapply(unique(DPAdw$conseqWeek), FUN = each_week)))
  return(weekcombined)
}
totalcombined <- as.data.table(rbindlist(pblapply(1:56, FUN = each_trap)))
st_geometry(totalcombined) <- st_sfc(totalcombined$geometry)

##Plot points - animate
suppressMessages(library(ggplot2))
suppressMessages(library(ggvoronoi))
suppressMessages(library(gganimate))
suppressMessages(library(gifski))
suppressMessages(library(png))
sfphenogeomp <- st_as_sf(phenogeomp)
p2 <- ggplot(data = insectgrid[1,]) +
  # geom_sf(data = subset(sfphenogeomp, conseqWeek==21), aes(fill=PercAreaPheno), lwd = 0) +
  geom_sf(data = totalcombined, aes(fill=WeightSumPhenoArea)) + 
  geom_point(aes(x = lon, y = lat), col = "blue", pch = 2) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 6)) +
  scale_fill_gradient(low="lightgreen", high="darkgreen", name="Weighted Area\n of Phenology") + 
  xlab("Longitude") + ylab("Latitude") +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) + 
  transition_manual(frames=conseqWeek) +
  labs(title = "Week: {current_frame}", fill = "Weighted Area of Phenology") +
  ease_aes("linear")
animp2 <- animate(plot = p2, fps = 2, renderer = gifski_renderer(), width = 800, height = 700)
anim_save(filename = "DistPhenoArea-Trap1.gif", animation = animp2, path="../Results/SpatioTemporal/")

# png(paste0("../Results/SpatioTemporal/DistPhenoArea-",segment,".png"), width=1800, height=1800, res=200)
p2
# dev.off()


###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###