#!/usr/bin/env Rscript
#####################################
#### CREATING PHENOWEIGHT SCRIPT ####
#####################################

### Merging phenological data with invertebrate biomass data

##Setting working directory
setwd("~/Documents/UniversityWork/PhD/ImageAnalysis/Code")
#For terminal arguments
args = commandArgs(trailingOnly=TRUE)

##Loading data
#Pixel data
suppressMessages(library(data.table))
# rawdata <- fread("../Data/LongTermStudy/PixelData/ManualITC.csv", header=T)
rawdata <- fread(args[1], header=T)

#Invertebrate data
insectbiomass <- read.csv("../../GroundSampling/InsectSampling/Data/InsectSamplingFinal.csv", header=T)

#Match data
match <- read.csv("../../GroundSampling/GroundPhenology/Data/VTreeMatching.csv")
#Location data
locations <- read.csv("../../GroundSampling/PointLocations/locationdistances.csv", header = T)

##Combining outlier data
tsdata <- read.csv("../Results/PhenologyValidation/ARIMA/ValidatedPhenologyARIMA.csv", header=T)
tsdata$Type <- "TS"
cpdata <- read.csv("../Results/PhenologyValidation/CP/ValidatedPhenologyCP.csv", header=T)
cpdata$Type <- "CP"
gmmdata <- read.csv("../Results/PhenologyValidation/GMM/ValidatedPhenologyGMM.csv", header=T)
gmmdata$ColourCoordinate <- "ALL"
gmmdata$Type <- "GMM"
validdata <- rbind(tsdata, cpdata, gmmdata)
#Adding frequency
validdata$Freq <- 1
#Adding DOYdiff
validdata$DOYdiff <- abs(validdata$DOY - validdata$ClosestDOY)

##Dead List
deadlist <- na.omit(match$Tree_Crown_ID[match$Dead=="Yes"])

###--------------------------------------------------------------------------------###

##Combining function
combining_datasets <- function(itctree){
  
  #Subset tree
  singletree <- subset(rawdata, Tree_Crown_ID==itctree)
  singlevalid <- subset(validdata, Tree_Crown_ID==itctree)
  
  ##Subset outlier type
  singlets <- subset(singlevalid, Type=="TS")
  singlecp <- subset(singlevalid, Type=="CP")
  singlegmm <- subset(singlevalid, Type=="GMM")
  
  #Apply function to add Phenology
  singletree$Phenology <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlevalid$DOY), as.character(singlevalid$ClosestPhenology[singlevalid$DOY==x]), NA)})
  singletree$CPWeighting <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlevalid$DOY), as.character(singlevalid$ClosestDOYWeighting[singlevalid$DOY==x]), NA)})
  singletree$PhenologyPerc <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlevalid$DOY), as.character(singlevalid$ClosestPhenologyPerc[singlevalid$DOY==x]), NA)})
  singletree$PhenologyTS <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlets$DOY), as.character(singlets$ClosestPhenology[singlets$DOY==x]), NA)})
  singletree$PhenologyCP <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlecp$DOY), as.character(singlecp$ClosestPhenology[singlecp$DOY==x]), NA)})
  singletree$PhenologyGMM <- sapply(singletree$conseqDOY, function(x){ifelse(any(x==singlegmm$DOY), as.character(singlegmm$ClosestPhenology[singlegmm$DOY==x]), NA)})
  
  #Creating Binary Phenology column
  singletree$BinPheno <- ifelse(!is.na(singletree$Phenology), 1, 0)
  
  ##Add Taxonomy
  if(length(match$Genus[match$Tree_Crown_ID %in% itctree])>0){
    if(is.na(match$Genus[match$Tree_Crown_ID %in% itctree])){
      if(is.na(match$Family[match$Tree_Crown_ID %in% itctree])){
        singletree$TreeFamily <- NA
        singletree$TreeGenus <- NA
        singletree$TreeSpecies <- NA
        singletree$TreeTaxon <- NA
      } else {
        singletree$TreeFamily <- match$Family[match$Tree_Crown_ID%in%itctree]
        singletree$TreeGenus <- NA
        singletree$TreeSpecies <- NA
        singletree$TreeTaxon <- NA
      }
    } else {
      singletree$TreeFamily <- match$Family[match$Tree_Crown_ID%in%itctree]
      singletree$TreeGenus <- match$Genus[match$Tree_Crown_ID%in%itctree]
      singletree$TreeSpecies <- match$Species[match$Tree_Crown_ID%in%itctree]
      #Add full taxonomic name
      if(unique(singletree$TreeSpecies)=="sp."){
        singletree$TreeTaxon <- paste0(unique(singletree$TreeGenus), " ", unique(singletree$TreeSpecies))
      } else {
        singletree$TreeTaxon <- paste0(substring(unique(singletree$TreeGenus), 1, 1), ". ", unique(singletree$TreeSpecies))
      }
    }
  } else {
    singletree$TreeFamily <- NA
    singletree$TreeGenus <- NA
    singletree$TreeSpecies <- NA
    singletree$TreeTaxon <- NA
  }

  ##Find Closest Trap and Weight
  #Matching Tree_Crown_ID (itctree) with VTree ID
  if(length(match$Tree_Crown_ID[match$Tree_Crown_ID%in%itctree])<1){
    vtree <- NA
    ClosestTrap <- rep(NA, NROW(singletree))
    DryWeight <- rep(NA, NROW(singletree))
    WetWeight <- rep(NA, NROW(singletree))
    singletree <- cbind(singletree, ClosestTrap, DryWeight, WetWeight)
    # print("Tree Unseen By Drone or No Matching Validation Tree")
  }else{
    if(match$Dead[match$Tree_Crown_ID %in% itctree] == "Yes"){
      # print(paste0("Tree Died on: ", unique(match$Date[match$Tree_Crown_ID==itctree])[1]))
      deadtree <- unique(match$Date[match$Tree_Crown_ID %in% itctree])
      deaddoy <- as.POSIXlt(deadtree, format = "%d/%m/%Y")$yday+1
      vtree <- match$VTree[match$Tree_Crown_ID %in% itctree]
      NewVTree <- ifelse(nchar(as.character(vtree)) == 1,  paste0("VTree00", vtree), 
                         ifelse(nchar(as.character(vtree)) == 2, paste0("VTree0", vtree), paste0("VTree",vtree)))
      subtree <- subset(locations, VTree == NewVTree)
      ClosestTrap <- subtree$IFT[subtree$HaversineDistance == min(subtree$HaversineDistance)][1]
      IFT <- subset(insectbiomass, Trap == ClosestTrap)
      IFT$timepoint <- 1:NROW(IFT)
      find_closest_weights <- function(doy){
        if(any(doy == singletree$conseqDOY[singletree$BinPheno==1])){
          tp <- which(abs(IFT$conseqDOY-doy) == min(abs(IFT$conseqDOY - doy)))
          DryWeight <- IFT$Dry_Weight[IFT$timepoint %in% tp]
          WetWeight <- IFT$Wet_Weight[IFT$timepoint %in% tp]
          return(cbind(doy, ClosestTrap, DryWeight, WetWeight))
        }else{
          DryWeight <- NA
          WetWeight <- NA
          return(cbind(doy, ClosestTrap, DryWeight, WetWeight))
        }
      }
      Weights <- as.data.frame(do.call("rbind", lapply(singletree$conseqDOY, FUN=find_closest_weights)))
      colnames(Weights)[1] <- "conseqDOY"
      Weights[,c(1,3,4)] <- sapply(Weights[,c(1,3,4)],as.numeric)
      singletree <- merge(singletree, Weights, by="conseqDOY")
      
    }else{
      deadtree <- NA
      deaddoy <- NA
      #Match VTree
      vtree <- match$VTree[match$Tree_Crown_ID %in% itctree]
      
      #Converting VTree format
      NewVTree <- ifelse(nchar(as.character(vtree)) == 1,  paste0("VTree00", vtree), 
                         ifelse(nchar(as.character(vtree)) == 2, paste0("VTree0", vtree), paste0("VTree",vtree)))
      
      #Subset the vtree in locations
      subtree <- subset(locations, VTree == NewVTree)
      
      #Find the closest trap
      ClosestTrap <- subtree$IFT[subtree$HaversineDistance %in% min(subtree$HaversineDistance)]
      
      #Subset the closest trap
      IFT <- subset(insectbiomass, Trap == ClosestTrap)
      
      #Adding timepoint
      IFT$timepoint <- 1:NROW(IFT)
      
      #Find Closest DOY and Weight of closest trap
      find_closest_weights <- function(doy){
        if(any(doy == singletree$conseqDOY[singletree$BinPheno==1])){
          tp <- which(abs(IFT$conseqDOY-doy) == min(abs(IFT$conseqDOY - doy)))
          DryWeight <- IFT$Dry_Weight[IFT$timepoint %in% tp]
          WetWeight <- IFT$Wet_Weight[IFT$timepoint %in% tp]
          return(cbind(doy, ClosestTrap, DryWeight, WetWeight))
        }else{
          DryWeight <- NA
          WetWeight <- NA
          return(cbind(doy, ClosestTrap, DryWeight, WetWeight))
        }
      }
      Weights <- as.data.frame(do.call("rbind", lapply(singletree$conseqDOY, FUN=find_closest_weights)))
      colnames(Weights)[1] <- "conseqDOY"
      Weights[,c(1,3,4)] <- sapply(Weights[,c(1,3,4)],as.numeric)
      singletree <- merge(singletree, Weights, by="conseqDOY")
    }
  }
  
  ##Changing Phenology for dead trees
  if(any(itctree == deadlist)){
    deaddate <- match$Date[match$Tree_Crown_ID==itctree][which(!is.na(match$Date[match$Tree_Crown_ID==itctree]))]
    deaddoy <- as.POSIXlt(deaddate, format = "%d/%m/%Y")$yday+1
    Year <- paste0("20", substring(as.character(as.POSIXlt(deaddate, format = "%d/%m/%Y")$year),2,3))
    conseqdeaddoy <- ifelse(Year==2020, deaddoy + 365, deaddoy)
    beforedead <- subset(singletree, conseqDOY < conseqdeaddoy)
    afterdead <- subset(singletree, conseqDOY >= conseqdeaddoy)
    afterdead$Phenology <- "Dead"
    afterdead$PhenologyTS <- "Dead"
    afterdead$PhenologyCP <- "Dead"
    afterdead$PhenologyGMM <- "Dead"
    singletree <- rbindlist(list(beforedead, afterdead))
  }
  
  #Recombining
  return(singletree)
}

##Adding phenology, phenology for each outlier, binary phenology and closest trap weight
suppressMessages(library(pbapply))
pbo <- pboptions(type="txt")
# phenoweight <- do.call("rbind", pblapply(unique(rawdata$Tree_Crown_ID), FUN=combining_datasets))
phenoweight <- as.data.frame(rbindlist(pblapply(unique(rawdata$Tree_Crown_ID), FUN=combining_datasets), use.names=T))

##Saving data
fwrite(phenoweight, "../Data/LongTermStudy/PixelData/PhenoWeight.csv", row.names = F, na="NA")
