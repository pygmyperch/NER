# ENMs
library(biomod2)
library(raster)
library(sp)
library(rgeos)
library(scales)


# Raw data sources
# https://chelsa-climate.org/bioclim/ # current
# https://chelsa-climate.org/future/ # future
# http://www.paleoclim.org # past (using late holocene 4.2-0.3 ka)


# set some path variables
# path to working directory
setwd("/PATH/TO/WD/")
path <- "/PATH/TO/WD/"

# path to directory containing maxent.jar
MaxEnt.path <- "/PATH/TO/maxent_dir/"


#####################################################################################
############################ get formatted raster stacks ############################
# path to rasters
setwd("/PATH/TO/data/environmental_rasters/")

for (i in c("EH", "MH", "LH", "current", "rcp45_2070", "rcp85_2070")) {
  files <- list.files(path = i, all.files = TRUE, recursive=F, full.names=T, pattern = "\\.grd$")
  rasters <- stack(c(files))
  assign(paste(i,".rasters", sep=""),rasters)
}

#####################################################################################
############################ get occurrence data ############################
# path to data directory
setwd("/PATH/TO/data")
#Import occurrence data for all species
DataSpecies <- read.csv("occurrence_data.csv", header=TRUE)


#######################################################################################
################################### Run analyses ##########################
# path to working directory
setwd("/PATH/TO/WD")

# Run analyses for each species
for (i in c("splendida", "eachamensis", "utcheensis", "malanda")) {
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



# set plot extent
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


# calculate relative range at each time point
EH.range <- as.matrix(table(as.numeric(as.vector(EH.binary_rep))))
MH.range <- as.matrix(table(as.numeric(as.vector(MH.binary_rep))))
LH.range <- as.matrix(table(as.numeric(as.vector(LH.binary_rep))))
current.range <- as.matrix(table(as.numeric(as.vector(current.binary_rep))))
rcp45_2070.range <- as.matrix(table(as.numeric(as.vector(rcp45_2070.binary_rep))))
rcp85_2070.range <- as.matrix(table(as.numeric(as.vector(rcp85_2070.binary_rep))))



# table of range changes over time
if (!dir.exists(paste0(path,"range"))) {dir.create(paste0(path,"range"))}
setwd(paste0(path,"range"))
range <- cbind(EH.range, MH.range, LH.range, current.range, rcp45_2070.range, rcp85_2070.range)
colnames(range) <- c("EH", "MH", "LH", "Current", "rcp45_2070", "rcp85_2070")
write.csv(range, paste0(i,"_range.csv"))
setwd(path)

}


# Plot Figure 2: changes in suitable habitat over time for each species
setwd(paste0(path,"range"))

malanda.range <- read.csv("malanda_range.csv", row.names = 1)
eachamensis.range <- read.csv("eachamensis_range.csv", row.names = 1)
utcheensis.range <- read.csv("utcheensis_range.csv", row.names = 1)
splendida.range <- read.csv("splendida_range.csv", row.names = 1)

dat <- as.matrix(rbind((malanda.range[2,]), (eachamensis.range[2,]), (utcheensis.range[2,]), (splendida.range[2,])))
colnames(dat) <- c("Early-Holocene", "Mid-Holocene", "Late-Holocene", "Current", "2070 RCP4.5", "2070 RCP8.5")

NE_colors <- ggpubfigs::friendly_pal("ito_seven")[c(5,2,1,4,3)]
cols <- NE_colors[c(1,2,4,3)]

leg.labs <- c(substitute(paste("Malanda")),
              substitute(paste(italic("M. eachamensis"))),
              substitute(paste(italic("M. utcheensis"))),
              substitute(paste(italic("M. splendida"))))


x <- 10 ^ (2:5)
labs <- scientific_format(1)(x)

setwd(path)

#setEPS()
pdf(file = "habitat_suitabilityMS.pdf", height = 6, width = 8)
par(oma = c(4, 2, 1, 1), bg = "#fdf6e3", fg="#073642")

barplot(dat, log="y",
        col=cols, mgp=c(2,-0.5,.5),
        border="#fdf6e3", 
        cex.names = 0.6,
        cex.axis = 0.6, ylim=c(10,100000),
        beside=TRUE, axes = FALSE,
        legend=FALSE,
        cex.lab=0.8,
        col.lab="#073642", col.axis="#073642")
axis(2, cex.axis = 0.6, at = c(10,100,1000,10000,100000), labels=c(0,labs), las=1, col="#073642", col.lab="#073642", col.axis="#073642")
title(ylab=expression("log(habitat area)"), line=2.5, cex.lab=0.7, col="#073642", col.lab="#073642")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
inset.x <- 0.1
inset.y <- 0.23

# plot legend
for (i in 1:4) {
  # start inset y at 0.15 and place each successive legend element 0.03 to the right
  inset.x <- inset.x+0.15
  
  legend("bottomleft", xpd = TRUE, horiz = TRUE, inset = c(inset.x,inset.y), legend = leg.labs[[i]], col=cols[i], 
         bty = "n", pch = 15, cex = 0.6, pt.cex = 1)
}

dev.off()




