#!/usr/bin/env Rscript
################################################
#### IDENTIFYING TREE CROWN OUTLIERS SCRIPT ####
################################################

### Using multiple methods to identify outliers in spectral time series for individual tree crowns

##VERSION 5

##Setting working directory
setwd("~/Documents/UniversityWork/PhD/ImageAnalysis/Code")
#For terminal arguments
args = commandArgs(trailingOnly=TRUE)
suppressMessages(library(parallel))

##Loading the SAFE data
print("Loading data.")
suppressMessages(library(data.table))
# rawdata <- fread("../Data/LongTermStudy/PixelData/Manual_Pixel_Data.csv", header = T)
rawdata <- fread(args[1], header=T)
# match <- read.csv("../../GroundSampling/GroundPhenology/Data/VTreeMatching.csv")

###--------------------------------------------------------------------------------###

###Data Wrangling
print("Data Wrangling.")

##Removing NAs
data <- na.omit(rawdata)

##Changing all dates to same format
suppressMessages(library(pbapply))
pbo <- pboptions(type="txt")
datevec <- character()
formdate <- function(i){
# for(i in as.character(data$Date)){
  if(grepl("/", i)==T){
    slashdate <- as.Date(i, "%d/%m/%Y")
    datevec <- c(datevec, as.character(slashdate))
  } else {
    dashdate <- as.Date(i, "%d-%m-%y")
    datevec <- c(datevec, as.character(dashdate))
  }
}
data$Date <- pbsapply(as.character(data$Date), FUN=formdate)

##Coverting to day of year
data$DOY <- (as.POSIXlt(data$Date, format = "%Y-%m-%d")$yday)+1

##Adding year
data$Year <- as.POSIXlt(data$Date, format = "%Y-%m-%d")$year
data$Year <- paste0("20", substring(as.character(data$Year),2,3))
#Subsetting 2019
data <- subset(data, !Year == "2018")

##Making DOY consequetive
data$conseqDOY <- ifelse(data$Year==2020, data$DOY + 365, data$DOY)

##Adding Week and consequtive week
data$Week <- strftime(as.POSIXlt(data$Date, format = "%Y-%m-%d"),format="%W")
data$Week <- as.numeric(as.character(data$Week))
data$conseqWeek <- ifelse(data$Year==2020, data$Week + 53, data$Week)

##Converting factors to numeric
if(!is.numeric(data$GCC)){data$GCC <- as.numeric(levels(data$GCC))[data$GCC]}
if(!is.numeric(data$RCC)){data$RCC <- as.numeric(levels(data$RCC))[data$RCC]}
if(!is.numeric(data$BCC)){data$BCC <- as.numeric(levels(data$BCC))[data$BCC]}
if(!is.numeric(data$ExG)){data$ExG <- as.numeric(levels(data$ExG))[data$ExG]}

##Renaming
allITC <- data

##Creating column of previous values
print("Calculating previous values for CCs.")
prev_vals <- function(tree){
  
  #Subset
  singletree <- subset(allITC, Tree_Crown_ID==tree)
  
  #Sort by consequtive DOY
  singletree <- singletree[order(singletree$conseqDOY),]
  
  #Previous values
  singletree$rccpreval <- c(NA, singletree$RCC[1:NROW(singletree)-1])
  singletree$gccpreval <- c(NA, singletree$GCC[1:NROW(singletree)-1])
  singletree$bccpreval <- c(NA, singletree$BCC[1:NROW(singletree)-1])
  
  return(singletree)
}
suppressMessages(library(pbapply))
pbo <- pboptions(type="txt")
temp <- do.call("rbind", pblapply(unique(allITC$Tree_Crown_ID), FUN=prev_vals))
allITC <- temp

##Calculate change for each colour coordinate
allITC$RCCchange <- allITC$RCC-allITC$rccpreval
allITC$GCCchange <- allITC$GCC-allITC$gccpreval
allITC$BCCchange <- allITC$BCC-allITC$bccpreval

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Find outlier using ARIMA models

print("Starting initial ARIMA model time series analysis.")
##Identify the outliers by... 
#Uses ARIMA models

##Time Series Outliers Package
suppressMessages(library("tsoutliers"))

##Function to identify ts outlier days for each tree
ts_outlier <- function(tree){
  
  #Subsetting tree
  singletree <- subset(allITC, Tree_Crown_ID==tree)
  
  #Order by conseqDOY
  singletree <- singletree[order(singletree$conseqDOY),]
  
  #Add time point
  singletree$timepoint <- 1:NROW(singletree)
  
  #Create featurelist
  featurelist <- as.list(singletree[,c("R_Mean", "G_Mean", "B_Mean",
                                       "R_StDev", "G_StDev", "B_StDev",
                                       "RCC", "GCC", "BCC",
                                       "RCCchange", "GCCchange", "BCCchange")])
  
  #Run ts outlier for each feature
  each_feature <- function(feature){
    
    #Create time series for each CC
    ts <- ts(unlist(featurelist[feature], use.names = F), frequency=1)
    
    #Find outliers in time series
    ts.try <- try(tso(ts), silent = T)
    
    #If error just add 0
    if(is(ts.try,"try-error")){
      TSOutlier <- rep(0, length(singletree$timepoint))
    }else{
      #Extracting times
      tsoutliers <- ts.try$times
      
      #Return out times if not 0
      addouts <- function(i){
        ifelse(is.null(i), return(0), return(i))
      }
      outs <- do.call("rbind", lapply(tsoutliers, FUN=addouts))
      
      #Finding outliers in each time series and adding one
      binscore <- function(i){
        ifelse(any(i==outs), return(1), return(0))
      }
      TSOutlier <- do.call("rbind", lapply(singletree$timepoint, FUN=binscore))
    }
    return(TSOutlier)
  }
  TSOutliers <- lapply(1:length(featurelist), FUN=each_feature)
  
  #Adding feature names
  names(TSOutliers) <- paste0(names(featurelist), "TSOutlier")
  
  #Combining with singletree
  singletree <- cbind(singletree, as.data.frame(TSOutliers))
  
  #Recombining
  return(singletree)
}

##Apply function to all trees and combine dataframe
# temp1 <- do.call("rbind", pblapply(unique(allITC$Tree_Crown_ID), FUN=ts_outlier))
# temp1 <- as.data.table(rbindlist(mclapply(unique(allITC$Tree_Crown_ID), FUN=ts_outlier, mc.cores = detectCores())))
temp1 <- as.data.table(rbindlist(pblapply(unique(allITC$Tree_Crown_ID), FUN=ts_outlier)))
allITC <- temp1

###--------------------------------------------------------------------------------###

##Identifying outlier Colour Coordinates
allITC$OutCC <- ifelse(allITC$RCCTSOutlier==1 & allITC$GCCTSOutlier==1 & allITC$BCCTSOutlier==1, "ALL",
                       ifelse(allITC$RCCTSOutlier==1 & allITC$GCCTSOutlier==1, "RCC-GCC",
                              ifelse(allITC$RCCTSOutlier==1 & allITC$BCCTSOutlier==1, "RCC-BCC",
                                     ifelse(allITC$GCCTSOutlier==1 & allITC$BCCTSOutlier==1, "GCC-BCC",
                                            ifelse(allITC$RCCTSOutlier==1, "RCC",
                                                   ifelse(allITC$GCCTSOutlier==1, "GCC",
                                                          ifelse(allITC$BCCTSOutlier==1, "BCC",
                                                                 ifelse(allITC$GCCTSOutlier==1, "ExG", NA))))))))

##Giving value of outlier to each CC
allITC$rccoutvalue <- ifelse(allITC$OutCC=="RCC", allITC$RCC,
                             ifelse(allITC$OutCC=="RCC-GCC", allITC$RCC,
                                    ifelse(allITC$OutCC=="RCC-BCC", allITC$RCC,
                                           ifelse(allITC$OutCC=="ALL", allITC$RCC, NA))))
allITC$gccoutvalue <- ifelse(allITC$OutCC=="GCC", allITC$GCC,
                             ifelse(allITC$OutCC=="RCC-GCC", allITC$GCC,
                                    ifelse(allITC$OutCC=="GCC-BCC", allITC$GCC,
                                           ifelse(allITC$OutCC=="ALL", allITC$GCC, NA))))
allITC$bccoutvalue <- ifelse(allITC$OutCC=="BCC", allITC$BCC,
                             ifelse(allITC$OutCC=="RCC-BCC", allITC$BCC,
                                    ifelse(allITC$OutCC=="GCC-BCC", allITC$BCC,
                                           ifelse(allITC$OutCC=="ALL", allITC$BCC,NA))))

##Giving change in value of outlier
allITC$RCCTSChange <- ifelse(allITC$OutCC=="RCC" |
                               allITC$OutCC=="RCC-GCC" |
                               allITC$OutCC=="RCC-BCC" |
                               allITC$OutCC=="ALL",
                             allITC$RCC-allITC$rccpreval, NA)
allITC$GCCTSChange <- ifelse(allITC$OutCC=="GCC" |
                               allITC$OutCC=="RCC-GCC" |
                               allITC$OutCC=="GCC-BCC" |
                               allITC$OutCC=="ALL",
                             allITC$GCC-allITC$gccpreval, NA)
allITC$BCCTSChange <- ifelse(allITC$OutCC=="BCC" |
                               allITC$OutCC=="RCC-BCC" |
                               allITC$OutCC=="GCC-BCC" |
                               allITC$OutCC=="ALL",
                             allITC$BCC-allITC$bccpreval, NA)

print("Finished initial ARIMA model time series analysis.")

###--------------------------------------------------------------------------------###

###Iteratively Fill ARIMA Model Outlier Gaps

print("Starting gap filling ARIMA model time series analysis.")

##Find first outliers
tsoutliers <- subset(allITC, 
                     R_MeanTSOutlier==1 | G_MeanTSOutlier==1 | B_MeanTSOutlier==1 | 
                       R_StDevTSOutlier==1 | G_StDevTSOutlier==1 | B_StDevTSOutlier==1 |
                       RCCTSOutlier==1 | GCCTSOutlier==1 | BCCTSOutlier==1| 
                       RCCchangeTSOutlier==1 | GCCchangeTSOutlier==1 | BCCchangeTSOutlier==1
)

##Remove previous outliers
final <- subset(allITC,
                R_MeanTSOutlier==1 | G_MeanTSOutlier==1 | B_MeanTSOutlier==1 | 
                  R_StDevTSOutlier==1 | G_StDevTSOutlier==1 | B_StDevTSOutlier==1 |
                  RCCTSOutlier==1 | GCCTSOutlier==1 | BCCTSOutlier==1| 
                  RCCchangeTSOutlier==1 | GCCchangeTSOutlier==1 | BCCchangeTSOutlier==1
)
remove <- allITC[!(allITC$R_MeanTSOutlier==1 | allITC$G_MeanTSOutlier==1 | allITC$B_MeanTSOutlier==1 | 
                     allITC$R_StDevTSOutlier==1 | allITC$G_StDevTSOutlier==1 | allITC$B_StDevTSOutlier==1 |
                     allITC$RCCTSOutlier==1 | allITC$GCCTSOutlier==1 | allITC$BCCTSOutlier==1| 
                     allITC$RCCchangeTSOutlier==1 | allITC$GCCchangeTSOutlier==1 | allITC$BCCchangeTSOutlier==1
),]
removed <- remove

###Iterate until finished.
repeat{
  ##Remove columns for next iteration
  iterallITC <- subset(removed, select = -c(R_MeanTSOutlier, G_MeanTSOutlier, B_MeanTSOutlier, 
                                            R_StDevTSOutlier, G_StDevTSOutlier, B_StDevTSOutlier,
                                            RCCTSOutlier, GCCTSOutlier, BCCTSOutlier, 
                                            RCCchangeTSOutlier, GCCchangeTSOutlier, BCCchangeTSOutlier,
                                            OutCC, rccoutvalue, gccoutvalue, bccoutvalue,        
                                            RCCTSChange, GCCTSChange, BCCTSChange))
  
  ##Run TSO again
  suppressMessages(library("tsoutliers"))
  ts_outlier <- function(tree){
    
    #Subsetting tree
    singletree <- subset(iterallITC, Tree_Crown_ID==tree)
    
    #Order by conseqDOY
    singletree <- singletree[order(singletree$conseqDOY),]
    
    #Add time point
    singletree$timepoint <- 1:NROW(singletree)
    
    #Create featurelist
    featurelist <- as.list(singletree[,c("R_Mean", "G_Mean", "B_Mean",
                                         "R_StDev", "G_StDev", "B_StDev",
                                         "RCC", "GCC", "BCC",
                                         "RCCchange", "GCCchange", "BCCchange")])
    
    #Run ts outlier for each feature
    each_feature <- function(feature){
      
      #Create time series for each CC
      ts <- ts(unlist(featurelist[feature], use.names = F), frequency=1)
      
      #Find outliers in time series
      ts.try <- try(tso(ts), silent = T)
      
      #If error just add 0
      if(is(ts.try,"try-error")){
        TSOutlier <- rep(0, length(singletree$timepoint))
      }else{
        #Extracting times
        tsoutliers <- ts.try$times
        
        #Return out times if not 0
        addouts <- function(i){
          ifelse(is.null(i), return(0), return(i))
        }
        outs <- do.call("rbind", lapply(tsoutliers, FUN=addouts))
        
        #Finding outliers in each time series and adding one
        binscore <- function(i){
          ifelse(any(i==outs), return(1), return(0))
        }
        TSOutlier <- do.call("rbind", lapply(singletree$timepoint, FUN=binscore))
      }
      return(TSOutlier)
    }
    TSOutliers <- lapply(1:length(featurelist), FUN=each_feature)
    
    #Adding feature names
    names(TSOutliers) <- paste0(names(featurelist), "TSOutlier")
    
    #Combining with singletree
    singletree <- cbind(singletree, as.data.frame(TSOutliers))
    
    #Recombining
    return(singletree)
  }
  print("Run ARIMA Model")
  suppressMessages(library("parallel"))
  iterallITC <- as.data.table(rbindlist(mclapply(unique(iterallITC$Tree_Crown_ID), FUN=ts_outlier, mc.cores = detectCores())))
  
  ##Identifying outlier Colour Coordinates
  iterallITC$OutCC <- ifelse(iterallITC$RCCTSOutlier==1 & iterallITC$GCCTSOutlier==1 & iterallITC$BCCTSOutlier==1, "ALL",
                             ifelse(iterallITC$RCCTSOutlier==1 & iterallITC$GCCTSOutlier==1, "RCC-GCC",
                                    ifelse(iterallITC$RCCTSOutlier==1 & iterallITC$BCCTSOutlier==1, "RCC-BCC",
                                           ifelse(iterallITC$GCCTSOutlier==1 & iterallITC$BCCTSOutlier==1, "GCC-BCC",
                                                  ifelse(iterallITC$RCCTSOutlier==1, "RCC",
                                                         ifelse(iterallITC$GCCTSOutlier==1, "GCC",
                                                                ifelse(iterallITC$BCCTSOutlier==1, "BCC",
                                                                       ifelse(iterallITC$GCCTSOutlier==1, "ExG", NA))))))))
  
  ##Giving value of outlier to each CC
  iterallITC$rccoutvalue <- ifelse(iterallITC$OutCC=="RCC"|
                                     iterallITC$OutCC=="RCC-GCC"|
                                     iterallITC$OutCC=="RCC-BCC"|
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$RCC, NA)
  iterallITC$gccoutvalue <- ifelse(iterallITC$OutCC=="GCC" |
                                     iterallITC$OutCC=="RCC-GCC"|
                                     iterallITC$OutCC=="GCC-BCC"|
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$GCC, NA)
  iterallITC$bccoutvalue <- ifelse(iterallITC$OutCC=="BCC" |
                                     iterallITC$OutCC=="RCC-BCC"|
                                     iterallITC$OutCC=="GCC-BCC"|
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$BCC,NA)
  ##Giving change in value of outlier
  iterallITC$RCCTSChange <- ifelse(iterallITC$OutCC=="RCC" | 
                                     iterallITC$OutCC=="RCC-GCC" | 
                                     iterallITC$OutCC=="RCC-BCC" |
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$RCC-iterallITC$rccpreval, NA)
  iterallITC$GCCTSChange <- ifelse(iterallITC$OutCC=="GCC" | 
                                     iterallITC$OutCC=="RCC-GCC" | 
                                     iterallITC$OutCC=="GCC-BCC" |
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$GCC-iterallITC$gccpreval, NA)
  iterallITC$BCCTSChange <- ifelse(iterallITC$OutCC=="BCC" | 
                                     iterallITC$OutCC=="RCC-BCC" | 
                                     iterallITC$OutCC=="GCC-BCC" |
                                     iterallITC$OutCC=="ALL", 
                                   iterallITC$BCC-iterallITC$bccpreval, NA)
  
  #Subset new outliers
  newtsoutliers <- subset(iterallITC, 
                          R_MeanTSOutlier==1 | G_MeanTSOutlier==1 | B_MeanTSOutlier==1 | 
                            R_StDevTSOutlier==1 | G_StDevTSOutlier==1 | B_StDevTSOutlier==1 |
                            RCCTSOutlier==1 | GCCTSOutlier==1 | BCCTSOutlier==1| 
                            RCCchangeTSOutlier==1 | GCCchangeTSOutlier==1 | BCCchangeTSOutlier==1
  )
  
  ##Break if zero
  if(NROW(newtsoutliers)==0){
    break
  }
  
  ##Remove previous outliers
  print("Subset Outliers")
  addfinal <- subset(iterallITC,  
                     R_MeanTSOutlier==1 | G_MeanTSOutlier==1 | B_MeanTSOutlier==1 | 
                       R_StDevTSOutlier==1 | G_StDevTSOutlier==1 | B_StDevTSOutlier==1 |
                       RCCTSOutlier==1 | GCCTSOutlier==1 | BCCTSOutlier==1| 
                       RCCchangeTSOutlier==1 | GCCchangeTSOutlier==1 | BCCchangeTSOutlier==1
  )
  final <- rbind(final, addfinal)
  removed <- iterallITC[!(iterallITC$R_MeanTSOutlier==1 | iterallITC$G_MeanTSOutlier==1 | iterallITC$B_MeanTSOutlier==1 | 
                            iterallITC$R_StDevTSOutlier==1 | iterallITC$G_StDevTSOutlier==1 | iterallITC$B_StDevTSOutlier==1 |
                            iterallITC$RCCTSOutlier==1 | iterallITC$GCCTSOutlier==1 | iterallITC$BCCTSOutlier==1| 
                            iterallITC$RCCchangeTSOutlier==1 | iterallITC$GCCchangeTSOutlier==1 | iterallITC$BCCchangeTSOutlier==1
  ),]
}

temp2 <- rbind(removed, final)
temp2 <- temp2[!duplicated(temp2),]
allITC <- temp2

###--------------------------------------------------------------------------------###

###Finding Changepoints

print("Starting Changepoint Analysis")
##Loading the changepoint analysis package
suppressMessages(library(changepoint))

#Remove timepoint 
allITC$timepoint <- NULL

##Function to identify changepoint days for each tree
cpt_outlier <- function(tree){
  
  #Subset
  singletree <- subset(allITC, Tree_Crown_ID==tree)
  
  #Order by conseqDOY
  singletree <- singletree[order(singletree$conseqDOY),]
  
  #Add time point
  singletree$timepoint <- 1:NROW(singletree)
  
  #Create featurelist
  featurelist <- as.list(singletree[,c("R_Mean", "G_Mean", "B_Mean",
                                       "R_StDev", "G_StDev", "B_StDev",
                                       "RCC", "GCC", "BCC",
                                       "RCCchange", "GCCchange", "BCCchange")])
  
  #Run changepoint analysis for each feature
  each_feature <- function(feature){
    
    #Create time series for each CC
    ts <- ts(unlist(featurelist[feature], use.names = F), frequency=1)
    
    #Find changepoints in time series using mean and variance
    if(is.na(ts[1])){
      ts[1] <- 0 #to remove NAs
      mvvalue <- cpt.meanvar(ts, method="PELT")
    } else {
      mvvalue <- cpt.meanvar(ts, method="PELT")
    }
    
    #Lag changepoint by one to find outlier
    cpts <- cpts(mvvalue)+1
    
    #Matching changepoints
    binscore <- function(i){
      ifelse(any(i==cpts), return(1), return(0))
    }
    cpts <- do.call("rbind", lapply(singletree$timepoint, FUN=binscore))
    
    #Recombining
    return(cpts)
  }
  CPOutliers <- lapply(1:length(featurelist), FUN=each_feature)
  
  #Adding feature names
  names(CPOutliers) <- paste0(names(featurelist), "CPOutlier")
  
  #Recombining
  singletree <- cbind(singletree, as.data.frame(CPOutliers))
  
  return(singletree)
}

##Apply function to all trees and combine dataframe
# temp4 <- as.data.table(rbindlist(mclapply(unique(allITC$Tree_Crown_ID), FUN=cpt_outlier, mc.cores=detectCores())))
temp4 <- as.data.table(rbindlist(pblapply(unique(allITC$Tree_Crown_ID), FUN=cpt_outlier)))
allITC <- temp4

###--------------------------------------------------------------------------------###

##Identifying outliers CC
allITC$cptsCC <- ifelse(allITC$RCCCPOutlier==1 & allITC$GCCCPOutlier==1 & allITC$BCCCPOutlier==1, "ALL",
                        ifelse(allITC$RCCCPOutlier==1 & allITC$GCCCPOutlier==1, "RCC-GCC",
                               ifelse(allITC$RCCCPOutlier==1 & allITC$BCCCPOutlier==1, "RCC-BCC",
                                      ifelse(allITC$RCCCPOutlier==1 & allITC$BCCCPOutlier==1, "GCC-BCC",
                                             ifelse(allITC$RCCCPOutlier==1, "RCC",
                                                    ifelse(allITC$GCCCPOutlier==1, "GCC",
                                                           ifelse(allITC$BCCCPOutlier==1, "BCC", NA)))))))

##Giving value of outlier
allITC$rcccpvalue <- ifelse(allITC$cptsCC=="RCC" | allITC$cptsCC=="RCC-GCC" |
                              allITC$cptsCC=="RCC-BCC" | allITC$cptsCC=="ALL", 
                            allITC$RCC, NA)
allITC$gcccpvalue <- ifelse(allITC$cptsCC=="GCC" | allITC$cptsCC=="RCC-GCC" | 
                              allITC$cptsCC=="GCC-BCC" | allITC$cptsCC=="ALL", 
                            allITC$GCC, NA)
allITC$bcccpvalue <- ifelse(allITC$cptsCC=="BCC" | allITC$cptsCC=="RCC-BCC" |
                              allITC$cptsCC=="GCC-BCC" | allITC$cptsCC=="ALL", 
                            allITC$BCC,NA)

##Giving change in value of outlier
allITC$RCCCPChange <- ifelse(allITC$cptsCC=="RCC" | 
                               allITC$cptsCC=="RCC-GCC" | 
                               allITC$cptsCC=="RCC-BCC" |
                               allITC$cptsCC=="ALL", 
                             allITC$RCC-allITC$rccpreval, NA)
allITC$GCCCPChange <- ifelse(allITC$cptsCC=="GCC" | 
                               allITC$cptsCC=="RCC-GCC" | 
                               allITC$cptsCC=="GCC-BCC" |
                               allITC$cptsCC=="ALL", 
                             allITC$GCC-allITC$gccpreval, NA)
allITC$BCCCPChange <- ifelse(allITC$cptsCC=="BCC" | 
                               allITC$cptsCC=="RCC-BCC" | 
                               allITC$cptsCC=="GCC-BCC" |
                               allITC$cptsCC=="ALL", 
                             allITC$BCC-allITC$bccpreval, NA)
print("Finished Changepoint Analysis")

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

###Multivariate analysis using Gaussian Mixture Modelling

print("Starting GMM Analysis.")
##Load packages and python script
suppressMessages(library(reticulate))
source_python("GMM.py")

##Replace NA with zero for change values
suppressMessages(library(data.table))
allITC$RCCchange <- nafill(allITC$RCCchange, fill=0)
allITC$GCCchange <- nafill(allITC$GCCchange, fill=0)
allITC$BCCchange <- nafill(allITC$BCCchange, fill=0)

##Finding anomolies using Gaussian Mixture Models
temp5 <- Run_GMM_Pheno(allITC)

##Match GMM Results by Tree_Crown_ID and DOY
allITC <- merge(allITC, temp5, by=c("Tree_Crown_ID", "DOY"))

##Run ARIMA on GMM Log Likelihoods
gmm_ts_outlier <- function(tree){
  
  #Subsetting tree
  singletree <- subset(allITC, Tree_Crown_ID==tree)
  
  #Order by conseqDOY
  singletree <- singletree[order(singletree$conseqDOY),]
  
  #Add timepoint
  singletree$timepoint <- 1:NROW(singletree)
  
  #Create featurelist
  featurelist <- as.list(singletree[,c("GMMprobscore_CC", "GMMprobscore_Mean", "GMMprobscore_StDev",
                                       "GMMprobscore_Change", "GMMprobscore_All")])
  
  #Run ts outlier for each feature
  each_feature <- function(feature){
    
    #Create time series for each CC
    ts <- ts(unlist(featurelist[feature], use.names = F), frequency=1)
    
    #Find outliers in time series
    ts.try <- try(tso(ts), silent = T)
    
    #If error just add 0
    if(is(ts.try,"try-error")){
      GMMTSOutlier <- rep(0, length(singletree$timepoint))
    }else{
      #Extracting times
      tsoutliers <- ts.try$times
      
      #Return out times if not 0
      addouts <- function(i){
        ifelse(is.null(i), return(0), return(i))
      }
      outs <- do.call("rbind", lapply(tsoutliers, FUN=addouts))
      
      #Finding outliers in each time series and adding one
      binscore <- function(i){
        ifelse(any(i==outs), return(1), return(0))
      }
      GMMTSOutlier <- do.call("rbind", lapply(singletree$timepoint, FUN=binscore))
    }
    return(GMMTSOutlier)
  }
  GMMTSOutliers <- lapply(1:length(featurelist), FUN=each_feature)
  
  #Adding feature names
  names(GMMTSOutliers) <- paste0("GMMTSOutlier", "_", unlist(strsplit(names(featurelist), "_"))[c(2,4,6,8,10)])
  
  #Combining with singletree
  singletree <- cbind(singletree, as.data.frame(GMMTSOutliers))
  
  #Recombining
  return(singletree)
}
temp6 <- as.data.table(rbindlist(pblapply(unique(allITC$Tree_Crown_ID), FUN=gmm_ts_outlier)))
allITC <- temp6

##Giving change in value of outlier
allITC$rccgmmvalue <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$RCC, NA)
allITC$gccgmmvalue <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$GCC, NA)
allITC$bccgmmvalue <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$BCC, NA)

##Giving change in value of outlier
allITC$RCCGMMChange <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$RCC-allITC$rccpreval, NA)
allITC$GCCGMMChange <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$GCC-allITC$gccpreval, NA)
allITC$BCCGMMChange <- ifelse(allITC$GMMTSOutlier_CC==1, allITC$BCC-allITC$bccpreval, NA)

print("Finished GMM Analysis.")

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###

##Saving file
# write.csv(allITC, file = "../Data/LongTermStudy/PixelData/allITC.csv", row.names = F)
# write.csv(allITC, file = "../Data/LongTermStudy/PixelData/VTreeITC.csv", row.names = F)
# write.csv(allITC, file = "../Data/LongTermStudy/PixelData/ManualITC.csv", row.names = F)
# write.csv(allITC, file = args[2], row.names = F)

###--------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------###
