# ENMs
library(biomod2)
library(raster)
library(sp)
library(rgeos)


# Data sources
# https://chelsa-climate.org/bioclim/ # current
# https://chelsa-climate.org/future/ # future
# http://www.paleoclim.org # past (using late holocene 4.2-0.3 ka)


# set some path variables
setwd("/Users/chrisbrauer/GoogleDrive/git_repos/NER/EnvironmentalNicheModels/")
path <- "/Users/chrisbrauer/GoogleDrive/git_repos/NER/EnvironmentalNicheModels/"
# set path to MaxEnt
MaxEnt.path <- "/Users/chrisbrauer/Analysis/maxent/"
#setwd("/PATH/TO/EnvironmentalNicheModels/")


# define function to fix potential errors in rasters by masking out any cells that contain NA for any layer
# intersect_mask <- function(x){
#   values_x <- getValues(x)
#   inter_x <- values_x %*% rep(1,nlayers(x))
#   mask <- setValues(subset(x,1),values = (inter_x>0))
#   return(mask)
# }


# define study area to crop all rasters
poly <- readWKT("POLYGON((142.5 -10.6, 144.7 -10.6, 151 -20, 153 -26, 149.2 -26, 145.8 -23.7, 143.1 -17.7, 141.4 -10.6, 142.5 -10.6))", p4s = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


#format and align rasters (not run)
######################################################################
# # current climate data
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/current/envicloud/chelsa/chelsa_V1/climatologies/bio")
# CHELSA_current.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawCurrent <- stack(c(CHELSA_current.files))
# current.rasters <- crop(RawCurrent, poly)
# current.rasters <- mask(current.rasters,poly)
# # mask out any cells that contain NA for any layer
# current.rasters <- stack(mask(current.rasters, intersect_mask(current.rasters)))
# 
# 
# ######### Past climate data (three time periods spanning 300ya-12ka)
# ######### late holocene (4.2-0.3 ka)
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/paleo/late_holocene/LH_v1_2_5m")
# CHELSA_LH.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawLH <- stack(c(CHELSA_LH.files))
# # clip to study area polygon
# LH.rasters <- crop(RawLH, poly)
# LH.rasters <- mask(LH.rasters,poly)
# # resample to match 30 second resolution of current and future rasters
# RawLH.hires <- resample(LH.rasters, current.rasters)
# 
# # mask out any cells that contain NA for any layer
# LH.rasters = stack(mask(RawLH.hires, intersect_mask(RawLH.hires)))
# 
# # check extents match
# stopifnot(all(extent(LH.rasters) == extent(current.rasters)))
# 
# ######### mid holocene (8.326-4.2 ka)
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/paleo/mid_holocene/MH_v1_2_5m")
# CHELSA_MH.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawMH <- stack(c(CHELSA_MH.files))
# # clip to study area polygon
# MH.rasters <- crop(RawMH, poly)
# MH.rasters <- mask(MH.rasters,poly)
# # resample to match 30 second resolution of current and future rasters
# RawMH.hires <- resample(MH.rasters, current.rasters)
# 
# # mask out any cells that contain NA for any layer
# MH.rasters = stack(mask(RawMH.hires, intersect_mask(RawMH.hires)))
# 
# # check extents match
# stopifnot(all(extent(MH.rasters) == extent(current.rasters)))
# 
# 
# 
# ######### early holocene (11.7-8.326 ka)
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/paleo/early_holocene/EH_v1_2_5m")
# CHELSA_EH.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawEH <- stack(c(CHELSA_EH.files))
# # clip to study area polygon
# EH.rasters <- crop(RawEH, poly)
# EH.rasters <- mask(EH.rasters,poly)
# # resample to match 30 second resolution of current and future rasters
# RawEH.hires <- resample(EH.rasters, current.rasters)
# 
# # mask out any cells that contain NA for any layer
# EH.rasters = stack(mask(RawEH.hires, intersect_mask(RawEH.hires)))
# 
# # check extents match
# stopifnot(all(extent(EH.rasters) == extent(current.rasters)))
# 
# 
# ######## future scenarios based on ACCESS1 circulation model
# # 2070 rcp45
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/2070/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/rcp45")
# rcp45_2070.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawRcp45 <- stack(c(rcp45_2070.files))
# rcp45_2070.rasters <- crop(RawRcp45, poly)
# rcp45_2070.rasters <- mask(rcp45_2070.rasters,poly)
# # mask out any cells that contain NA for any layer
# rcp45_2070.rasters <- stack(mask(rcp45_2070.rasters, intersect_mask(rcp45_2070.rasters)))
# 
# # check extents match
# stopifnot(all(extent(rcp45_2070.rasters) == extent(current.rasters)))
# 
# 
# # 2070 rcp85
# setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/2070/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/rcp85")
# rcp85_2070.files <- list.files(all.files = TRUE, recursive=F, full.names=T, pattern = "\\.tif$")
# RawRcp85 <- stack(c(rcp85_2070.files))
# rcp85_2070.rasters <- crop(RawRcp85, poly)
# rcp85_2070.rasters <- mask(rcp85_2070.rasters,poly)
# # mask out any cells that contain NA for any layer
# rcp85_2070.rasters <- stack(mask(rcp85_2070.rasters, intersect_mask(rcp85_2070.rasters)))
# 
# # check extents match
# stopifnot(all(extent(rcp85_2070.rasters) == extent(current.rasters)))
# 
# 
# # make layer names the same for all data (two sets as original data not all ordered the same way (CHELSA current data has padded 0s, future and PaleoClim data does not))
# # i.e. bio01, bio02... vs. bio1, bio10...
# bioclim.CHELSA.names <- c("bio_01", "bio_02", "bio_03", "bio_04", "bio_05", "bio_06", "bio_07", "bio_08", "bio_09",
#                            "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18", "bio_19")
# 
# bioclim.PaleoClim.names <- c("bio_01",  "bio_10", "bio_11", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_18",
#                              "bio_19", "bio_02",  "bio_03",  "bio_04",  "bio_05",  "bio_06",  "bio_07",  "bio_08",  "bio_09")
# 
# 
# names(current.rasters) <- bioclim.CHELSA.names
# names(LH.rasters) <- bioclim.PaleoClim.names
# names(MH.rasters) <- bioclim.PaleoClim.names
# names(EH.rasters) <- bioclim.PaleoClim.names
# names(rcp45_2070.rasters) <- bioclim.PaleoClim.names
# names(rcp85_2070.rasters) <- bioclim.PaleoClim.names
# 
# # save cropped and aligned data
# for (i in c("EH", "MH", "LH", "current", "rcp45_2070", "rcp85_2070")) {
#   dir.create(paste0(i),  showWarnings = FALSE)
# }
# writeRaster(EH.rasters, "EH/EH.rasters.grd", bylayer=TRUE, suffix=bioclim.PaleoClim.names)
# writeRaster(MH.rasters, "MH/MH.rasters.grd", bylayer=TRUE, suffix=bioclim.PaleoClim.names)
# writeRaster(LH.rasters, "LH/LH.rasters.grd", bylayer=TRUE, suffix=bioclim.PaleoClim.names)
# writeRaster(current.rasters, "current/current.rasters.grd", bylayer=TRUE, suffix=bioclim.CHELSA.names)
# writeRaster(rcp45_2070.rasters, "rcp45_2070/rcp45_2070.rasters.grd", bylayer=TRUE, suffix=bioclim.PaleoClim.names)
# writeRaster(rcp85_2070.rasters, "rcp85_2070/rcp85_2070.rasters.grd", bylayer=TRUE, suffix=bioclim.PaleoClim.names)
#########################################

# load data
#Import occurrence data for all species
DataSpecies <- read.csv("occurrence_data.csv", header=TRUE)

# import rasters
for (i in c("current", "EH", "MH", "LH", "rcp45_2070", "rcp85_2070")) {
  dir.create(paste0(i),  showWarnings = FALSE)
  
raster.files <- list.files(path = paste0(path,"environmental_rasters/",i), all.files = TRUE, recursive=F, full.names=T, pattern = "\\.grd$")
rasters <- stack(c(raster.files))
assign(paste0(i,".rasters"),rasters)

}


# make dir tree
for (i in c("splendida", "eachamensis", "utcheensis", "Malanda")) {
  dir.create(paste0(i),  showWarnings = FALSE)


# aggregate occurence data in each dir
setwd(paste0(path,i))
occurrence <- DataSpecies[DataSpecies[i] == 1, ]
myRespName <- i
myResp <- as.numeric(occurrence[,i])
myRespXY <- occurrence[,c("X","Y")]



# current climate
coords <- data.frame(x=myRespXY$X,y=myRespXY$Y)
points <- SpatialPoints(coords, proj4string = current.rasters@crs)
env.values <- extract(current.rasters,points)

df <- cbind.data.frame(coordinates(points),env.values)
write.csv(df, paste0(i,"_occurrences_current_climate.csv"), row.names = FALSE)


# get retained variable names to subset temporal rasterstacks
retained.env <- c("bio_05", "bio_19")


current.env <- as.data.frame(cbind(df$bio_05, df$bio_19))
colnames(current.env) <- c("bio_05", "bio_19")

# data for ENMs

reduced.env <- subset(current.rasters, retained.env)


#####################################################
### run models

#Set up modelling, including adding pseudoabsences
myBiomodData <- BIOMOD_FormatingData(resp.var= myResp,
                                     expl.var = reduced.env,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName, 
                                     na.rm = TRUE, 
                                     PA.nb.rep = 3, 
                                     PA.nb.absences = 383, 
                                     PA.strategy = "random")
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT.Phillips = list(memory_allocated=2048, path_to_maxent.jar=MaxEnt.path))

#Run models and evaluations
# maxent + three most common algorithms from:
# Hao, T., Elith, J., Guillera‐Arroita, G., & Lahoz‐Monfort, J. J. (2019). A review of evidence about use and performance of
# species distribution modelling ensembles like BIOMOD. Diversity and Distributions, 25(5), 839-852.
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c("MAXENT.Phillips","GLM", "GBM", "RF"),
                                    models.options = myBiomodOption,
                                    NbRunEval = 3,
                                    DataSplit = 80,
                                    Prevalence = 0.5,
                                    VarImport = 3,
                                    models.eval.meth = c("ROC"),
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE,
                                    do.full.models = FALSE,
                                    modeling.id = myRespName)

#Evaluate fit and distribution of models
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
myBiomodModelEval["ROC","Testing.data",,,] #reports ROC for all 3 pseudoreplicates for all models
var_imp <- get_variables_importance(myBiomodModelOut) #reports importance of bioclim variables
write.csv(var_imp, "variable_importance.csv")



#Generate ensemble model of from all methods
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by = 'all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.8),
  prob.mean = TRUE,
  prob.cv = TRUE,
  prob.ci = TRUE,
  prob.ci.alpha = 0.05,
  prob.median = TRUE,
  committee.averaging = TRUE,
  prob.mean.weight = TRUE,
  prob.mean.weight.decay = 'proportional')
myBiomodEM
myBiomodEMEval <- get_evaluations(myBiomodEM, as.data.frame=T)

write.csv(myBiomodEMEval, "ensemble_eval.csv")

#Project map to current climatic conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.env,
  proj.name = "current",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = FALSE,
  output.format = ".grd")
myCurrentProj <- get_predictions(myBiomodProj)

#Generate ensemble maps and outputs
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj)
proj.stk <- get_predictions(myBiomodEF)
writeRaster(subset(proj.stk, 1), filename=paste0(i,"_CurrentMeanEM.asc"), overwrite = TRUE)



# plot weighted mean by ROC model
MyXFrm <- 140
MyXTo <- 155
MyXStp <- 5
MyYFrm <- -26
MyYTo <- -10
MyYStp <- 2



setEPS()
postscript(file = paste0(i,"_ENM_current.eps"), height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stk[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "(a)", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=xmin(proj.stk),las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_current_all_plots.eps"), height = 10, width = 20)
plot(proj.stk)
dev.off()


#####################################################
#####################################################
# late Holocene projection
LH.values <- extract(LH.rasters,points)
LH.df <- cbind.data.frame(coordinates(points),LH.values)
write.csv(LH.df, "occurrences_LH_climate.csv", row.names = FALSE)

LH.env <- as.data.frame(cbind(LH.df$bio_05, LH.df$bio_19))
colnames(LH.env) <- c("bio_05", "bio_19")

# data for ENMs
reduced.LH <- subset(LH.rasters, retained.env)


myBiomodProjLH <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.LH,
  proj.name = "LateHolocene",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = TRUE,
  output.format = ".grd")

#Run projections under past climatic conditions using ensemble model
myBiomodEFLH <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjLH)
proj.stkLH <- get_predictions(myBiomodEFLH)
writeRaster(subset(proj.stkLH, 1), filename=paste0(i,"_LateHoloMeanEM.asc"), overwrite=TRUE)


setEPS()
postscript(file = paste0(i,"_ENM_LH.eps"), height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stkLH[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "LH", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_LH_all_plots.eps"), height = 10, width = 20)
plot(proj.stkLH)
dev.off()



#####################################################
#####################################################
# mid Holocene projection
MH.values <- extract(MH.rasters,points)
MH.df <- cbind.data.frame(coordinates(points),MH.values)
write.csv(MH.df, "occurrences_MH_climate.csv", row.names = FALSE)

MH.env <- as.data.frame(cbind(MH.df$bio_05, MH.df$bio_19))
colnames(MH.env) <- c("bio_05", "bio_19")

# data for ENMs
reduced.MH <- subset(MH.rasters, retained.env)


myBiomodProjMH <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.MH,
  proj.name = "MidHolocene",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = TRUE,
  output.format = ".grd")

#Run projections under past climatic conditions using ensemble model
myBiomodEFMH <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjMH)
proj.stkMH <- get_predictions(myBiomodEFMH)
writeRaster(subset(proj.stkMH, 1), filename=paste0(i,"_MidHoloMeanEM.asc"), overwrite=TRUE)


setEPS()
postscript(file = "MsplendidaENM_MH.eps", height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stkMH[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "(a)", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_MH_all_plots.eps"), height = 10, width = 20)
plot(proj.stkMH)
dev.off()



#####################################################
#####################################################
# Early Holocene projection
EH.values <- extract(EH.rasters,points)
EH.df <- cbind.data.frame(coordinates(points),EH.values)
write.csv(EH.df, "occurrences_EH_climate.csv", row.names = FALSE)

EH.env <- as.data.frame(cbind(EH.df$bio_05, EH.df$bio_19))
colnames(EH.env) <- c("bio_05", "bio_19")

# data for ENMs
reduced.EH <- subset(EH.rasters, retained.env)


myBiomodProjEH <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.EH,
  proj.name = "EarlyHolocene",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = TRUE,
  output.format = ".grd")

#Run projections under past climatic conditions using ensemble model
myBiomodEFEH <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjEH)
proj.stkEH <- get_predictions(myBiomodEFEH)
writeRaster(subset(proj.stkEH, 1), filename=paste0(i,"_EarlyHoloMeanEM.asc"), overwrite=TRUE)


setEPS()
postscript(file = "MsplendidaENM_EH.eps", height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stkEH[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "(a)", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_EH_all_plots.eps"), height = 10, width = 20)
plot(proj.stkEH)
dev.off()




#####################################################
## future projections
#####################################################
# rcp45_2070
rcp45_2070.values <- extract(rcp45_2070.rasters,points)
rcp45_2070.df <- cbind.data.frame(coordinates(points),rcp45_2070.values)
write.csv(rcp45_2070.df, "occurrences_rcp45_2070_climate.csv", row.names = FALSE)

rcp45_2070.env <- as.data.frame(cbind(rcp45_2070.df$bio_05, rcp45_2070.df$bio_19))
colnames(rcp45_2070.env) <- c("bio_05", "bio_19")

# data for ENMs
reduced.rcp45_2070 <- subset(rcp45_2070.rasters, retained.env)


myBiomodProjrcp45_2070 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.rcp45_2070,
  proj.name = "rcp45_2070",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = TRUE,
  output.format = ".grd")

#Run projections under future climatic conditions using ensemble model
myBiomodEFrcp45_2070 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjrcp45_2070)
proj.stkrcp45_2070 <- get_predictions(myBiomodEFrcp45_2070)
writeRaster(subset(proj.stkrcp45_2070, 1), filename=paste0(i,"_rcp45_2070MeanEM.asc"), overwrite=TRUE)


setEPS()
postscript(file = paste0(i,"_ENM_rcp45_2070.eps"), height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stkrcp45_2070[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "rcp45_2070", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_rcp45_2070_all_plots.eps"), height = 10, width = 20)
plot(proj.stkrcp45_2070)
dev.off()




# rcp85_2070
rcp85_2070.values <- extract(rcp85_2070.rasters,points)
rcp85_2070.df <- cbind.data.frame(coordinates(points),rcp85_2070.values)
write.csv(rcp85_2070.df, "occurrences_rcp85_2070_climate.csv", row.names = FALSE)

rcp85_2070.env <- as.data.frame(cbind(rcp85_2070.df$bio_05, rcp85_2070.df$bio_19))
colnames(rcp85_2070.env) <- c("bio_05", "bio_19")

# data for ENMs
reduced.rcp85_2070 <- subset(rcp85_2070.rasters, retained.env)


myBiomodProjrcp85_2070 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = reduced.rcp85_2070,
  proj.name = "rcp85_2070",
  selected.models = "all",
  binary.meth = "ROC",
  compress = "xz",
  clamping.mask = TRUE,
  output.format = ".grd")

#Run projections under past climatic conditions using ensemble model
myBiomodEFrcp85_2070 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjrcp85_2070)
proj.stkrcp85_2070 <- get_predictions(myBiomodEFrcp85_2070)
writeRaster(subset(proj.stkrcp85_2070, 1), filename=paste0(i,"_rcp85_2070MeanEM.asc"), overwrite=TRUE)


setEPS()
postscript(file = paste0(i,"_ENM_rcp85_2070.eps"), height = 6, width = 6)
par(oma=c(6, 5, 5, 6) + 0.1, mar = c(0.5,0.5,2,0.5) + 0.1)#c(bottom, left, top, right)
plot(proj.stkrcp85_2070[[1]], axes=FALSE, box=FALSE, legend=TRUE, xlab="",ylab="")
title(main = "rcp85_2070", adj = 0, cex.main = 1.5)
# Plot the axes
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)

title(xlab = "Longitude",
      ylab = "Latitude",
      outer = TRUE, line = 3, cex.lab = 0.8)


dev.off()

setEPS()
postscript(file = paste0(i,"_rcp85_2070_all_plots.eps"), height = 10, width = 20)
plot(proj.stkrcp85_2070)
dev.off()



####################################################
### plot composite figure for all ENMs ###
####################################################

setEPS()
postscript(file = paste0(i,"_ENMsMS.eps"), height = 6, width = 9)
par(mfrow = c(2,3), oma=c(5, 5, 5, 6) + 0.1, mar = c(2,0.1,2,0.1) + 0.1)#it goes c(bottom, left, top, right)
z <- c(0,1)# set legend limits

# EH
plot(proj.stkEH[[1]]/1000, main = "EH", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)


# MH
plot(proj.stkMH[[1]]/1000, main = "MH", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)


# LH
plot(proj.stkLH[[1]]/1000, main = "LH", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)


# current
plot(proj.stk[[1]]/1000, main = "current", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)


# rcp45_2070
plot(proj.stkrcp45_2070[[1]]/1000, main = "rcp45_2070", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)


# rcp85_2070
plot(proj.stkrcp85_2070[[1]]/1000, main = "rcp85_2070", zlim=z, legend = F, axes = F, box = F)
axis(1,tick=TRUE,pos=ymin(proj.stk),las=1,at=seq(MyXFrm,MyXTo ,MyXStp ), cex=0.8) 
axis(2,tick=TRUE,pos=MyXFrm,las=1,at=seq(MyYFrm ,MyYTo ,MyYStp ), cex=0.8)



#reset par to single plot to place single legend and title
par(mfrow=c(1,1),new=FALSE, oma=c(0,0,0,0))
plot(proj.stkrcp85_2070[[1]]/1000,legend.only=TRUE ,legend.shrink=0.5, legend.width=0.3, zlim=z)
title(main = paste0(i," ENMs"))

dev.off()






# set presence threshold at 70% for binary representation
EH.binary_rep <- BinaryTransformation(proj.stkEH[[1]], 700)
MH.binary_rep <- BinaryTransformation(proj.stkMH[[1]], 700)
LH.binary_rep <- BinaryTransformation(proj.stkLH[[1]], 700)
current.binary_rep <- BinaryTransformation(proj.stk[[1]], 700)
rcp45_2070.binary_rep <- BinaryTransformation(proj.stkrcp45_2070[[1]], 700)
rcp85_2070.binary_rep <- BinaryTransformation(proj.stkrcp85_2070[[1]], 700)


# calculate range at each time point
EH.range <- as.matrix(table(as.numeric(as.vector(EH.binary_rep))))
EH.range

MH.range <- as.matrix(table(as.numeric(as.vector(MH.binary_rep))))
MH.range

LH.range <- as.matrix(table(as.numeric(as.vector(LH.binary_rep))))
LH.range

current.range <- as.matrix(table(as.numeric(as.vector(current.binary_rep))))
current.range

rcp45_2070.range <- as.matrix(table(as.numeric(as.vector(rcp45_2070.binary_rep))))
rcp45_2070.range

rcp85_2070.range <- as.matrix(table(as.numeric(as.vector(rcp85_2070.binary_rep))))
rcp85_2070.range


# table of range change over time
range <- cbind(EH.range, MH.range, LH.range, current.range, rcp45_2070.range, rcp85_2070.range)
colnames(range) <- c("EH", "MH", "LH", "Current", "rcp45_2070", "rcp85_2070")
write.csv(range, paste0(i,"_range.csv"))

}
