##########################################################################################################################################
################# Combination D Multi Surface Redo Optimization 1 #######################################################################################
##########################################################################################################################################

### STEP 1: Load packages

options(bitmapType='cairo')

library(sp)
library(doParallel)
library(spdep)
library(reshape2)
library(raster)
library(devtools)
library(ResistanceGA)
library(parallel)

### STEP 2: Set working directory

setwd(".")

### STEP 3: Load the sample location data and create a spatial points object

Vv.localities <- read.csv("Data/Vv_Coord_Unrelated.csv")
Vv.localities.spatial <- SpatialPoints(Vv.localities[,c(3,4)]) # sample.locales


### STEP 4: Create vector of genetic distances

Ar.matrix <- read.csv("Data/Ar_Distances_Unrelated.csv", row.names = 1)
Ar.matrix.subset <- subset(melt(Ar.matrix), value!=0)
Ar.matrix.vector <- Ar.matrix.subset[,2]


### STEP 5: Create source folder for copying ASCII files

source.folder <- paste(getwd(), "/Data/ASCII_Files", sep="")

### STEP 6: Create new directory for results

dir.create(file.path(paste(getwd(), "/Results/Multi_Unrelated/", sep=""), "MS_Unrelated_Combo_D_Results_Opt1"))

### STEP 7: Paste .asc files into the newly created sub-directory

combo.d <- c("TPI_final.asc", "TRI_final.asc")
ss.results.opt1 <- paste(getwd(), "/Results/Multi_Unrelated/MS_Unrelated_Combo_D_Results_Opt1/", sep="")

file.copy(file.path(source.folder, combo.d), ss.results.opt1)

### STEP 8: Set directory to write results of first optimization

write.dir.opt1 <- paste(getwd(), "/Results/Multi_Unrelated/MS_Unrelated_Combo_D_Results_Opt1/", sep="")


### STEP 9: Prepare data for optimization No. 1
  
gdist.inputs <- gdist.prep(n.Pops = length(Vv.localities.spatial),
                           samples = Vv.localities.spatial,
                           response = Ar.matrix.vector,
                           method = 'commuteDistance')

GA.inputs <- GA.prep(method = "LL",
                     ASCII.dir = write.dir.opt1,
                     Results.dir = write.dir.opt1,
                     max.cat = 10000,
                     max.cont = 10000,
                     select.trans = list('A', 'A'), #TPI, TRI
                     parallel = 32) 


### STEP 10: Run all single surface optimization

MS.results.opt1 <- MS_optim(gdist.inputs = gdist.inputs,
                             GA.inputs = GA.inputs)
