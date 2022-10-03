#Evaluating ecological niche of the genus Tympanoctomys and their species 
#Line 8-118: Tympanoctomys
#Line 120-200: T. barrerae
#Line 230-330:T. kirchnerorum (current records) 
#Line 334-430:T. kirchnerorum (fossil records, CCSM4) 
#Line 430-520:T. kirchnerorum (fossil records, MIROC-ESM) 

#Genus Tympanoctomys
rm(list=ls()) 
library(sp)
library(raster)
library(dismo)# Needs maxent.jar in Java carpet
library(rJava)
library(rgdal)
library(maptools)
library(rgeos)
library(spam)
library(grid)
library(maps)
library(fields)
library(parallel)
library(ENMeval)
options(digits=5)
setwd("C:/Users/USUARIO/Google Drive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data S3/Tympanoctomys")

#Configurating
options(java.parameters = "-Xmx1g")

#Upload current ocurrences of the genus
Tympanoctomys <- read.csv("Tympanoctomys.csv",header=TRUE, sep=',', stringsAsFactors=F)
Tympanoctomys
locTympanoctomys <-cbind(Tympanoctomys$x,Tympanoctomys$y)
plot(locTympanoctomys)

loc_Tympanoctomys <- list.files(pattern = "*.csv")

for (i in 1:length(loc_Tympanoctomys)) assign(loc_Tympanoctomys[i], read.csv(loc_Tympanoctomys[i])) 

#Apply spatial thinning if exists GBIF data

#Import bioclimatic data (predictor variables)
env_layers <- list.files(pattern='asc', full.names=TRUE)
env_layers
env_predictors <- stack(env_layers) # making a stack
projection(env_predictors) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') 
env_predictors # resume
names(env_predictors) # names of vars
plot(env_predictors,1) # see that stack
points(cbind(locTympanoctomys, pch=24,col='black',cex=0.7))

#Create a MCP
PolMC <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2]) #chull is from grDevices
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

PMC_Tympanoctomys <- PolMC(locTympanoctomys) # Tympanoctomys barrerae
plot(PMC_Tympanoctomys) # Visualizar el PMC

#Creating a buffer
Buff_PMC_Tympanoctomys <- gBuffer(PMC_Tympanoctomys, width = .9) # 0.9 degrees
plot (Buff_PMC_Tympanoctomys)

#All visualizations
plot(env_predictors,1) 
plot(PMC_Tympanoctomys,add=T) # visualizar el PMC
plot(Buff_PMC_Tympanoctomys, add=T) # visualizar el PMC y su buffer
points(locTympanoctomys,pch=24,col='black',cex=0.5) # visualizar los puntos de presencia de una especie

#Cutting environments from the buffer
pred_Tympanoctomys_box <- crop(env_predictors, extent(Buff_PMC_Tympanoctomys)) 
pred_Tympanoctomys <- mask(pred_Tympanoctomys_box, Buff_PMC_Tympanoctomys) 
plot(pred_Tympanoctomys) # visualize

#Saving
writeRaster(pred_Tympanoctomys, "C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data SD3/Tympanoctomys/Tympanoctomys_vars", bylayer=TRUE, format="ascii")

#Put maxent.jar file (from the dowloaded software) in Rversion/dismo/java

#Creating range of settings, current environments
Tympanoctomys_test <- ENMevaluate(locTympanoctomys, pred_Tympanoctomys, bg.coords = NULL, 
                                  RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "LQH" ), method = "block", overlap=F) # FC: (L)lineal, (Q)quadratic and (H)hinge

##*** Running ENMevaluate using maxnet v.0.1.2 ***
##  Doing evaluations using spatial blocks...
##|======================================================================================| 100%
##ENMeval completed in 19 minutes 49.3 seconds.

metricsvalues_Tympanoctomys <- Tympanoctomys_test@results
metricsvalues_Tympanoctomys

#Write a table
write.table(metricsvalues_Tympanoctomys, file="metrics_Tympanoctomys.csv", row.names=FALSE, col.names=TRUE)
metrics_Tympanoctomys <- read.csv("metrics_Tympanoctomys.csv",header=TRUE, sep=',', stringsAsFactors=F)
View(metrics_Tympanoctomys)

#Checking results
Tympanoctomys_test@predictions #prediction propierties
Tympanoctomys_test@occ.pts # presence points for the models
Tympanoctomys_test@bg.pts # background cells

# Which settings gave delta.AICc < 2?
aicmods <- which(Tympanoctomys_test@results$delta.AICc < 2)
Tympanoctomys_test@results[aicmods,]
##settings features  rm train.AUC avg.test.AUC var.test.AUC avg.diff.AUC var.diff.AUC
##3  LQH_0.5      LQH 0.5 0.9475394     0.626225    0.1072642    0.3206112    0.1280434
##avg.test.orMTP var.test.orMTP avg.test.or10pct var.test.or10pct     AICc delta.AICc w.AIC
##3      0.4423077      0.1405325            0.625       0.09307199 1796.984          0     1
##parameters
##3         53

# View predictions in geographic space for these models
plot(Tympanoctomys_test@predictions[[aicmods]])

# Plot delta.AICc Vs. RM's
eval.plot(Tympanoctomys_test@results)

#========================================
#T. barrerae
rm(list=ls())
options(digits=5)
setwd("C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data S3/Tbarrerae")

library(sp)
library(raster)
library(dismo)# Needs maxent.jar in Java carpet
library(rJava)
library(rgdal)
library(maptools)
library(rgeos)
library(spam)
library(grid)
library(maps)
library(fields)
library(parallel)
library(ENMeval)

#Configurating
options(java.parameters = "-Xmx1g")

#Upload all ocurrences of the species: T.barrerae
Tbarrerae <- read.csv("Tbarrerae.csv",header=TRUE, sep=',', stringsAsFactors=F)
Tbarrerae
locbarrerae <-cbind(Tbarrerae$x,Tbarrerae$y)
plot(locbarrerae)

loc_barrerae <- list.files(pattern = "*.csv")

for (i in 1:length(loc_barrerae)) assign(loc_barrerae[i], read.csv(loc_barrerae[i])) 

#Apply spatial thinning if exists GBIF data

#Import bioclimatic data (predictor variables)
env_layers <- list.files(pattern='asc', full.names=TRUE)
env_layers
env_predictors <- stack(env_layers) # making a stack
projection(env_predictors) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') 
env_predictors # resume
names(env_predictors) # names of vars
plot(env_predictors,1) # see that stack
points(cbind(locbarrerae, pch=24,col='black',cex=0.7))

#Create a MCP
PolMC <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2]) #chull is from grDevices
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

PMC_barrerae <- PolMC(locbarrerae) # Tympanoctomys barrerae
plot(PMC_barrerae) # Visualizar el PMC

#Creating a buffer
Buff_PMC_barrerae <- gBuffer(PMC_barrerae, width = .9) # 0.9 degrees
plot (Buff_PMC_barrerae)

#All visualizations
plot(env_predictors,1) 
plot(PMC_barrerae,add=T) # visualizar el PMC
plot(Buff_PMC_barrerae, add=T) # visualizar el PMC y su buffer
points(locbarrerae,pch=24,col='black',cex=0.5) # visualizar los puntos de presencia de una especie

#Cutting environments from the buffer
pred_barrerae_box <- crop(env_predictors, extent(Buff_PMC_barrerae)) 
pred_barrerae <- mask(pred_barrerae_box, Buff_PMC_barrerae) 
plot(pred_barrerae) # visualize

#Saving
writeRaster(pred_barrerae, "C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data SD3/Tbarrerae/barrerae_vars", bylayer=TRUE, format="ascii")

#Put maxent.jar file (from the dowloaded software) in Rversion/dismo/java

#Creating range of settings 
barrerae_test <- ENMevaluate(locbarrerae, pred_barrerae, bg.coords = NULL, 
                             RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH"), method = "block", overlap=F) # FC: lineal, quadratic and hinge

##** Running ENMevaluate using maxnet v.0.1.2 ***
##  Doing evaluations using spatial blocks...
##|==================                                                                    |  100%
#ENMeval completed in 23 minutes 1.5 seconds.

metricsvalues_barrerae <- barrerae_test@results
metricsvalues_barrerae

#Write a table
write.table(metricsvalues_barrerae, file="metrics_barrerae.csv", row.names=FALSE, col.names=TRUE)
metrics_barrerae <- read.csv("metrics_barrerae.csv",header=TRUE, sep=',', stringsAsFactors=F)
View(metrics_barrerae)

#Checking results
barrerae_test@predictions #prediciton propierties
barrerae_test@occ.pts # presence points for the models
barrerae_test@bg.pts # background cells

# Which settings gave delta.AICc < 2?
aicmods <- which(barrerae_test@results$delta.AICc < 2)
barrerae_test@results[aicmods,]
##settings features  rm train.AUC avg.test.AUC var.test.AUC avg.diff.AUC var.diff.AUC
##4  LQH_0.5      LQH 0.5    0.9662       0.5378     0.030967      0.42826     0.030938
##avg.test.orMTP var.test.orMTP avg.test.or10pct var.test.or10pct   AICc delta.AICc w.AIC
##4        0.53084       0.044381          0.80465          0.07307 1322.5          0     1
##parameters
##4         45

# View predictions in geographic space for these models
plot(barrerae_test@predictions[[aicmods]])

# Plot delta.AICc Vs. RM's
eval.plot(barrerae_test@results)

#================================================================================
#Species T. kirchnerorum with current records
rm(list=ls())
setwd("C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data S3/Tkirchnerorum")
library(sp)
library(raster)
library(dismo)# Needs maxent.jar in Java carpet
library(rJava)
library(rgdal)
library(maptools)
library(rgeos)
library(spam)
library(grid)
library(maps)
library(fields)
library(parallel)
library(ENMeval)
#Configurating
options(java.parameters = "-Xmx1g")

#Upload all ocurrences of the species: T.kirchnerorum
Tkirchnerorum <- read.csv("Tkirchnerorum.csv",header=TRUE, sep=',', stringsAsFactors=F)
Tkirchnerorum
lockirchnerorum <-cbind(Tkirchnerorum$x,Tkirchnerorum$y)
plot(lockirchnerorum)

loc_kirchnerorum <- list.files(pattern = "*.csv")

for (i in 1:length(loc_kirchnerorum)) assign(loc_kirchnerorum[i], read.csv(loc_kirchnerorum[i])) 

#Apply spatial thinning if exists GBIF data

#Import bioclimatic data (predictor variables)
env_layers <- list.files(pattern='asc', full.names=TRUE)
env_layers
env_predictors <- stack(env_layers) # making a stack
projection(env_predictors) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # Caracteristicas del stack
env_predictors # resume
names(env_predictors) # names of vars
plot(env_predictors,1) # see that stack
points(cbind(lockirchnerorum, pch=24,col='black',cex=0.7))

#Create a MCP
PolMC <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2]) #chull is from grDevices
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

PMC_kirchnerorum <- PolMC(lockirchnerorum) # Tympanoctomys kirchnerorum
plot(PMC_kirchnerorum) # Visualizar el PMC

#Creating a buffer
Buff_PMC_kirchnerorum<- gBuffer(PMC_kirchnerorum, width = .9) # 0.9 degrees
plot (Buff_PMC_kirchnerorum)

#All visualizations
plot(env_predictors,1) 
plot(PMC_kirchnerorum,add=T) # visualizar el PMC
plot(Buff_PMC_kirchnerorum, add=T) # visualizar el PMC y su buffer
points(lockirchnerorum,pch=24,col='black',cex=0.5) # visualizar los puntos de presencia de una especie

#Cutting environments from the buffer
pred_kirchnerorum_box <- crop(env_predictors, extent(Buff_PMC_kirchnerorum)) 
pred_kirchnerorum <- mask(pred_kirchnerorum_box, Buff_PMC_kirchnerorum) 
plot(pred_kirchnerorum) # visualize

#Saving
writeRaster(pred_kirchnerorum, "C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data SD3/Tkirchnerorum/kirchnerorum_vars", bylayer=TRUE, format="ascii")
#Put maxent.jar file (from the dowloaded software) in Rversion/dismo/java

#Creating range of settings
kirchnerorum_test <- ENMevaluate(lockirchnerorum, pred_kirchnerorum, bg.coords = NULL, 
                                 RMvalues=seq(0.5, 4, 0.5), fc=c("L"), method = "block", overlap=F) # FC: lineal, quadratic and hinge
##*** Running ENMevaluate using maxnet v.0.1.2 ***
##  Doing evaluations using spatial blocks...
##|======================================================================================| 100%
##ENMeval completed in 2 minutes 45.6 seconds.

metricsvalues_kirchnerorum <- kirchnerorum_test@results
metricsvalues_kirchnerorum

#Write a table
write.table(metricsvalues_kirchnerorum, file="metrics_kirchnerorum_jackknife.csv", row.names=FALSE, col.names=TRUE)

#Checking results
kirchnerorum_test@predictions #prediction propierties
kirchnerorum_test@occ.pts # presence points for the models
kirchnerorum_test@bg.pts # background cells

# Which settings gave delta.AICc < 2?
aicmods <- which(kirchnerorum_test@results$delta.AICc < 2)
kirchnerorum_test@results[aicmods,]
##settings features rm train.AUC avg.test.AUC var.test.AUC avg.diff.AUC var.diff.AUC
##16     LQ_4       LQ  4 0.9688296    0.8861396   0.07161257   0.08234086   0.06102016
##avg.test.orMTP var.test.orMTP avg.test.or10pct var.test.or10pct     AICc delta.AICc
##16          0.125         0.0625            0.125           0.0625 127.6252          0
##w.AIC parameters
##16 0.8776257          5

# View predictions in geographic space for these models
plot(kirchnerorum_test@predictions[[aicmods]])

# Plot delta.AICc Vs. RM's
eval.plot(kirchnerorum_test@results)

#========================================
#Upload only fossil ocurrences of the species: T.kirchnerorum (CCSM4 model)
setwd("C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data S3/Tkirchnerorum/MH/ccsm")
Tkirchnerorumf <- read.csv("Tkirchnerorumf.csv",header=TRUE, sep=',', stringsAsFactors=F)
Tkirchnerorumf
lockirchnerorumf <-cbind(Tkirchnerorumf$x,Tkirchnerorumf$y)
plot(lockirchnerorumf)

loc_kirchnerorumf <- list.files(pattern = "*.csv")

for (i in 1:length(loc_kirchnerorumf)) assign(loc_kirchnerorumf[i], read.csv(loc_kirchnerorumf[i])) 

#Apply spatial thinning if exists GBIF data

#Import bioclimatic data (predictor variables)
env_layers <- list.files(pattern='asc', full.names=TRUE)
env_layers
env_predictors <- stack(env_layers) # making a stack
projection(env_predictors) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') 
env_predictors # resume
names(env_predictors) # names of vars
plot(env_predictors,1) # see that stack
points(cbind(lockirchnerorumf, pch=24,col='black',cex=0.7))

#Create a MCP
PolMC <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2]) #chull is from grDevices
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

PMC_lockirchnerorumf <- PolMC(lockirchnerorumf) # Tympanoctomys barrerae
plot(PMC_lockirchnerorumf) # Visualizar el PMC

#Creating a buffer
Buff_PMC_lockirchnerorumf <- gBuffer(PMC_lockirchnerorumf, width = .9) # 0.9 degrees
plot (Buff_PMC_lockirchnerorumf)

#All visualizations
plot(env_predictors,1) 
plot(PMC_lockirchnerorumf,add=T) # visualizar el PMC
plot(Buff_PMC_lockirchnerorumf, add=T) # visualizar el PMC y su buffer
points(lockirchnerorumf,pch=24,col='black',cex=0.5) # visualizar los puntos de presenica de una especie

#Cutting environments from the buffer
pred_lockirchnerorumf_box <- crop(env_predictors, extent(Buff_PMC_lockirchnerorumf)) 
pred_lockirchnerorumf <- mask(pred_lockirchnerorumf_box, Buff_PMC_lockirchnerorumf) 
plot(pred_lockirchnerorumf) # visualize

#Saving
writeRaster(pred_lockirchnerorumf, "C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data SD3/Tkirchnerorum/MH/ccsm/kirchnerorumf_vars", bylayer=TRUE, format="ascii")

#Put maxent.jar file (from the dowloaded software) in Rversion/dismo/java

#Creating range of settings 
kirchnerorumf_test <- ENMevaluate(lockirchnerorumf, pred_lockirchnerorumf, bg.coords = NULL, 
                             RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH"), method = "jackknife", overlap=F) # FC: lineal, quadratic and hinge
##*** Running ENMevaluate using maxnet v.0.1.2 ***
##Doing evaluations using k-1 jackknife...
##|==========================================================================| 100%
##ENMeval completed in 147 minutes 4.5 seconds.

metricsvalues_kirchnerorumf <- kirchnerorumf_test@results
metricsvalues_kirchnerorumf

#Write a table
write.table(metricsvalues_kirchnerorumf, file="metrics_kirchnerorumf_jackknife.csv", row.names=FALSE, col.names=TRUE)

#Checking results
kirchnerorumf_test@predictions #prediciton propierties
kirchnerorumf_test@occ.pts # presence points for the models
kirchnerorumf_test@bg.pts # background cells

# Which settings gave delta.AICc < 2?
aicmods <- which(kirchnerorumf_test@results$delta.AICc < 2)
kirchnerorumf_test@results[aicmods,]
##  settings features rm train.AUC avg.test.AUC var.test.AUC avg.diff.AUC
##8    LQH_1      LQH  1 0.9787091    0.9691636    0.1362451   0.02536201
##var.diff.AUC avg.test.orMTP var.test.orMTP avg.test.or10pct var.test.or10pct
##8    0.1264567     0.04545455     0.04545455       0.09090909       0.08658009
##AICc delta.AICc     w.AIC parameters
##8 339.2062          0 0.9981182         14

# View predictions in geographic space for these models
plot(kirchnerorumf_test@predictions[[aicmods]])

# Plot delta.AICc Vs. RM's
eval.plot(kirchnerorumf_test@results)

#========================================
#Upload only fossil ocurrences of the species: T.kirchnerorum (MIROC-ESM model)
setwd("C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data S3/Tkirchnerorum/MH/miroc")
Tkirchnerorumf <- read.csv("Tkirchnerorumf.csv",header=TRUE, sep=',', stringsAsFactors=F)
Tkirchnerorumf
lockirchnerorumf <-cbind(Tkirchnerorumf$x,Tkirchnerorumf$y)
plot(lockirchnerorumf)

loc_kirchnerorumf <- list.files(pattern = "*.csv")

for (i in 1:length(loc_kirchnerorumf)) assign(loc_kirchnerorumf[i], read.csv(loc_kirchnerorumf[i])) 
#Apply spatial thinning if exists GBIF data

#Import bioclimatic data (predictor variables)
env_layers <- list.files(pattern='asc', full.names=TRUE)
env_layers
env_predictors <- stack(env_layers) # making a stack
projection(env_predictors) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') 
env_predictors # resume
names(env_predictors) # names of vars
plot(env_predictors,1) # see that stack
points(cbind(lockirchnerorumf, pch=24,col='black',cex=0.7))

#Create a MCP
PolMC <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2]) #chull is from grDevices
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

PMC_lockirchnerorumf <- PolMC(lockirchnerorumf) # Tympanoctomys barrerae
plot(PMC_lockirchnerorumf) # Visualizar el PMC

#Creating a buffer
Buff_PMC_lockirchnerorumf <- gBuffer(PMC_lockirchnerorumf, width = .9) # 0.9 degrees
plot (Buff_PMC_lockirchnerorumf)

#All visualizations
plot(env_predictors,1) 
plot(PMC_lockirchnerorumf,add=T) # visualizar el PMC
plot(Buff_PMC_lockirchnerorumf, add=T) # visualizar el PMC y su buffer
points(lockirchnerorumf,pch=24,col='black',cex=0.5) # visualizar los puntos de presenica de una especie

#Cutting environments from the buffer
pred_lockirchnerorumf_box <- crop(env_predictors, extent(Buff_PMC_lockirchnerorumf)) 
pred_lockirchnerorumf <- mask(pred_lockirchnerorumf_box, Buff_PMC_lockirchnerorumf) 
plot(pred_lockirchnerorumf) # visualize

#Saving
writeRaster(pred_lockirchnerorumf, "C:/Users/USUARIO/OneDrive/Spatial_ecology_RVR/PAPERS/JM/REVISION/Supplementary_files/Supplementary Data SD3/Tkirchnerorum/MH/miroc/kirchnerorumf_vars", bylayer=TRUE, format="ascii")
#Put maxent.jar file (from the dowloaded software) in Rversion/dismo/java

#Creating range of settings 
kirchnerorumf_test <- ENMevaluate(lockirchnerorumf, pred_lockirchnerorumf, bg.coords = NULL, 
                                  RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH"), method = "jackknife", overlap=F) # FC: lineal, quadratic and hinge
##*** Running ENMevaluate using maxnet v.0.1.2 ***
##  Doing evaluations k-1 jackknife...
##There are 75 background records with NA for at least one predictor variable.
##Removing these records from analysis, resulting in 9925 records...
##|============================================================================|   100%
##77 grid cells found with at least one NA value: these cells were excluded from raster predictions.
##ENMeval completed in 147 minutes 4.6 seconds.

metricsvalues_kirchnerorumf <- kirchnerorumf_test@results
metricsvalues_kirchnerorumf

#Write a table
write.table(metricsvalues_kirchnerorumf, file="metrics_kirchnerorumf_jackknife.csv", row.names=FALSE, col.names=TRUE)

#Checking results
kirchnerorumf_test@predictions #prediciton propierties
kirchnerorumf_test@occ.pts # presence points for the models
kirchnerorumf_test@bg.pts # background cells

# Which settings gave delta.AICc < 2?
aicmods <- which(kirchnerorumf_test@results$delta.AICc < 2)
kirchnerorumf_test@results[aicmods,]
## settings features rm train.AUC avg.test.AUC var.test.AUC avg.diff.AUC var.diff.AUC
##8    LQH_1      LQH  1 0.9800092    0.9706435    0.1388969   0.02264713    0.1339742
##avg.test.orMTP var.test.orMTP avg.test.or10pct var.test.or10pct     AICc
##8     0.04545455     0.04545455       0.09090909       0.08658009 322.0852
##delta.AICc     w.AIC parameters
##8          0 0.9532234         12

# View predictions in geographic space for these models
plot(kirchnerorumf_test@predictions[[aicmods]])

# Plot delta.AICc Vs. RM's
eval.plot(kirchnerorumf_test@results)

#========================================