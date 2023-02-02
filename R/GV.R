# genomic vulnerability


library(adegenet)
library(raster)
library(AlleleShift)
library(dplyr)
library(stringi)
library(tidyverse)



# define modified AlleleShift plotting function to use custom colours etc.
population.shift_cb <- function (baseline.env.data, future.env.data, title, from, to, option = c("PCA", "RDA"), vector.multiply = 1) {
  require(grid)
  message("Checking data sets with BiodiversityR::check.datasets")
  BiodiversityR::check.datasets(baseline.env.data, future.env.data)
  pca.input <- rbind(baseline.env.data, future.env.data)
  np <- nrow(baseline.env.data)
  pca.env <- data.frame(climate = factor(c(rep("baseline", 
                                               np), rep("future", np))))
  if (option == "PCA") {
    message("Fitting PCA with vegan::rda and no explanatory variables")
    pca.result <- vegan::rda(pca.input, scale = TRUE)
    summary(pca.result)
    axis.long1 <- BiodiversityR::axis.long(pca.result, choices = c(1, 
                                                                   2))
    plot1 <- vegan::ordiplot(pca.result, xlim = c(-1, 3), ylim = c(-2.5, 1.5))
  }
  if (option == "RDA") {
    message("Fitting RDA with vegan::rda and time period (baseline / future) as explanatory variable")
    rda.result <- vegan::rda(pca.input ~ climate, data = pca.env, 
                             scale = TRUE)
    summary(rda.result)
    axis.long1 <- BiodiversityR::axis.long(rda.result, choices = c(1, 
                                                                   2))
    plot1 <- vegan::ordiplot(rda.result)
  }
  species.long1 <- BiodiversityR::species.long(plot1)
  pca.env <- cbind(pca.input, pca.env)
  sites.long1 <- BiodiversityR::sites.long(plot1, env.data = pca.env)
  segment.long1 <- sites.long1[sites.long1$climate == "baseline", 
                               c("labels", "axis1", "axis2")]
  segment.long2 <- sites.long1[sites.long1$climate == "future", 
                               c("labels", "axis1", "axis2")]
  segment.long3 <- data.frame(cbind(segment.long1, segment.long2))
  BioR.theme <- ggplot2::theme(plot.background = ggplot2::element_rect(fill="#fdf6e3"), panel.background = ggplot2::element_rect(fill="#fdf6e3"), 
                               panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), 
                               axis.line = ggplot2::element_line("#073642"), text = ggplot2::element_text(size = 12), 
                               axis.text = ggplot2::element_text(size = 10, colour = "#073642"), 
                               axis.title = ggplot2::element_text(size = 10, colour = "#073642"), 
                               plot.title = ggplot2::element_text(hjust = 0.5),
                               legend.title = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 10), 
                               legend.key = ggplot2::element_blank())
  plotggx <- ggplot2::ggplot() + ggplot2::ggtitle(title) + ggplot2::geom_vline(xintercept = c(0), 
                                                                               color = "#073642", linetype = 2) + ggplot2::geom_hline(yintercept = c(0), 
                                                                                                                                      color = "#073642", linetype = 2) + ggplot2::xlab(axis.long1[1, 
                                                                                                                                                                                                  "label"]) + ggplot2::ylab(axis.long1[2, "label"]) + ggplot2::scale_x_continuous(breaks = c(-1, 0.0, 1.0, 2.0), limits=c(-1, 2.5), sec.axis = ggplot2::dup_axis(breaks=NULL, labels = NULL, 
                                                                                                                                                                                                                                                                                                                                                                 name = NULL)) + ggplot2::scale_y_continuous(breaks = c(-1.6, -0.8, 0.0, 0.8, 1.6), limits=c(-1.6, 1.7), sec.axis = ggplot2::dup_axis(breaks=NULL, labels = NULL, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      name = NULL)) + ggforce::geom_mark_ellipse(data = sites.long1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ggplot2::aes(x = sites.long1$axis1, y = sites.long1$axis2, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              colour = sites.long1$climate), fill = ggplot2::alpha("#07364266", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0.4), expand = 0, size = 0.5, show.legend = FALSE) + 
    
    ggplot2::geom_segment(data = species.long1, ggplot2::aes(x = 0, 
                                                             y = 0, xend = species.long1$axis1 * vector.multiply, 
                                                             yend = species.long1$axis2 * vector.multiply), colour = "#073642", alpha = 0.5,
                          size = 0.5, arrow = ggplot2::arrow(angle = 20, length = unit(0.1, "inches"), type = "open")) + ggplot2::geom_segment(data = segment.long3, 
                                                                                                                                               ggplot2::aes(x = segment.long3$axis1, y = segment.long3$axis2, 
                                                                                                                                                            xend = segment.long3$axis1.1, yend = segment.long3$axis2.1), 
                                                                                                                                               colour = ind_cols, size = 0.3, arrow = ggplot2::arrow(angle = 20, length = unit(0.05, "inches"), type = "closed")) + 
    ggplot2::geom_text(data = species.long1, ggplot2::aes(x = species.long1$axis1 * 
                                                            vector.multiply, y = species.long1$axis2 * vector.multiply, 
                                                          label = species.long1$labels), colour = "#073642", size = 3, position = position_nudge(x = 0.2, y = c(0.1, 0.1)),
                       fontface = "plain", 
                       show.legend = FALSE) + ggplot2::geom_text(data = subset(sites.long1, 
                                                                               sites.long1$climate == "baseline"), ggplot2::aes(x = subset(sites.long1, 
                                                                                                                                           sites.long1$climate == "baseline")$axis1, y = subset(sites.long1, 
                                                                                                                                                                                                sites.long1$climate == "baseline")$axis2, label = subset(sites.long1, 
                                                                                                                                                                                                                                                         sites.long1$climate == "baseline")$labels, colour = subset(sites.long1, 
                                                                                                                                                                                                                                                                                                                    sites.long1$climate == "baseline")$climate), alpha = 0.7, 
                                                                 colour = ind_cols, size = 3, show.legend = FALSE) + BioR.theme + 
    
    
    # legend
    ggplot2::scale_colour_discrete(name  ="Climate",
                                   breaks=c("baseline", "future"),
                                   labels=c(from, to)) +
    ggplot2::scale_shape_discrete(name  ="Climate",
                                  breaks=c("baseline", "future"),
                                  labels=c(from, to)) +
    
    ggplot2::coord_fixed(ratio = 1)
  return(plotggx)
}

# define function to extract allele frequency deltas
get.deltaAF <- function (freq.future) {
  
  # function modified from AlleleShift::shift.dot.ggplot
  
  np <- length(unique(freq.future$Pop))
  
  freq.columns <- cbind(freq.future$Allele.freq, freq.future$Freq.e2, 
                        freq.future$LCL, freq.future$UCL)
  freq.means <- stats::aggregate(freq.columns ~ freq.future$Pop, 
                                 FUN = stats::median)
  names(freq.means)[1] <- "Pop"
  freq.means <- dplyr::left_join(freq.future[, c("Pop", 
                                                 "Allele.freq")], freq.means, by = "Pop")
  freq.future$Allele.freq <- freq.means$V1
  freq.future$Freq.e2 <- freq.means$V2
  freq.future$LCL <- freq.means$V3
  freq.future$UCL <- freq.means$V4
  freq.future <- freq.future[1:np, ]
  freq.future$increasing <- freq.future$Freq.e2 > freq.future$Allele.freq
  freq.future$Allele <- factor(rep("Alleles aggregated", 
                                   np))
  
  #freq.future$Pop <- droplevels(freq.future$Pop)
  
  #deltaAF <- as.data.frame(cbind(levels(freq.future$Pop), freq.future$Allele.freq-freq.future$Freq.e2))
  #colnames(deltaAF) <- c("Pop", "deltaAF")
  
  deltaAF <- as.data.frame(freq.future$Allele.freq-freq.future$Freq.e2)
  colnames(deltaAF) <- c("deltaAF")
  
  return(deltaAF)
}

setwd("/PATH/TO/data")
load("genomic_vulnerability.RData")


############# get admixture proportions and define hybrid and pure inds #############
#####
# get site averages from ADMIXTURE Qmatrix and reorder as per genetic data
admix <- read.table("NE342MS_LD300_bed.5.Q")
admix <- cbind(NEpops[,2],admix)
colnames(admix) <- c("site", "splendida", "Malanda", "eachamensis", "utcheensis", "Tully")


# calculate mean q values per site
mean_admix <- admix %>% group_by(site) %>% summarize(splendida = mean(splendida), Malanda = mean(Malanda), 
                                                     eachamensis = mean(eachamensis), utcheensis = mean(utcheensis), Tully = mean(Tully))

# check rows match (doesn't work with pop number AF matrix)
#stopifnot(all(rownames(NEaf) == mean_admix[,1]))
#stopifnot(all(rownames(NEaf) == as.character(pop.data[,1])))
stopifnot(all(mean_admix[,1] == as.character(pop.data[,1])))

# get hybrid and pure pops
eacham_splendida.hybrids <- mean_admix[mean_admix$splendida > 0.1 & mean_admix$eachamensis > 0.1 & mean_admix$Tully < 0.05 & mean_admix$utcheensis < 0.05 & mean_admix$Malanda < 0.05, ]
Malanda_splendida.hybrids <- mean_admix[mean_admix$splendida > 0.1 & mean_admix$eachamensis < 0.05 & mean_admix$Tully < 0.05 & mean_admix$utcheensis < 0.05 & mean_admix$Malanda > 0.1, ]
utchee_splendida.hybrids <- mean_admix[mean_admix$splendida > 0.1 & mean_admix$eachamensis < 0.05 & mean_admix$Tully < 0.05 & mean_admix$utcheensis > 0.1 & mean_admix$Malanda < 0.05, ]
Tully_splendida.hybrids <- mean_admix[mean_admix$splendida > 0.1 & mean_admix$eachamensis < 0.05 & mean_admix$Tully > 0.1 & mean_admix$utcheensis < 0.05 & mean_admix$Malanda < 0.05, ]

splendida.pops <- mean_admix[mean_admix$splendida > 0.95, ]
eacham.pops <- mean_admix[mean_admix$eachamensis > 0.95, ]
Malanda.pops <- mean_admix[mean_admix$Malanda > 0.95, ]
utchee.pops <- mean_admix[mean_admix$utcheensis > 0.95, ]
Tully.pops <- mean_admix[mean_admix$Tully > 0.95, ]

pure.pops <- as.data.frame(c(splendida.pops$site, eacham.pops$site, Malanda.pops$site, utchee.pops$site, Tully.pops$site))
hybrid.pops <- as.data.frame(c(eacham_splendida.hybrids$site, Malanda_splendida.hybrids$site, utchee_splendida.hybrids$site))


#####################################################################################
############################ get formatted raster stacks ############################
setwd("/PATH/TO/data/environmental_rasters/")

for (i in c("EH", "MH", "LH", "current", "rcp45_2070", "rcp85_2070")) {
  files <- list.files(path = i, all.files = TRUE, recursive=F, full.names=T, pattern = "\\.grd$")
  rasters <- stack(c(files))
  assign(paste(i,".rasters", sep=""),rasters)
}



# format environmental data for each time period
# current
current_var <- cbind(pop.data$bio_05, pop.data$bio_19)
colnames(current_var) <- c("bio_05", "bio_19")
rownames(current_var) <- pop.data$popnum
current_var <- as.data.frame(current_var)


# EH
env.site.EH <- raster::extract(EH.rasters,sites)
env.site.EH <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.EH)
write.csv(env.site.EH, "EH_site_env_data.csv", row.names = FALSE)

# Re-order env data to pop afs
EH.data <- cbind(unique(NEpops[,2]),env.site.EH)
EH.data <- EH.data[order(EH.data[,1]),]


EH_var <- cbind(EH.data$bio_05, EH.data$bio_19)
colnames(EH_var) <- c("bio_05", "bio_19")
rownames(EH_var) <- pop.data$popnum
EH_var <- as.data.frame(EH_var)


# MH
env.site.MH <- raster::extract(MH.rasters,sites)
env.site.MH <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.MH)
write.csv(env.site.MH, "MH_site_env_data.csv", row.names = FALSE)

# Re-order env data to pop afs
MH.data <- cbind(unique(NEpops[,2]),env.site.MH)
MH.data <- MH.data[order(MH.data[,1]),]


MH_var <- cbind(MH.data$bio_05, MH.data$bio_19)
colnames(MH_var) <- c("bio_05", "bio_19")
rownames(MH_var) <- pop.data$popnum
MH_var <- as.data.frame(MH_var)

# LH
env.site.LH <- raster::extract(LH.rasters,sites)
env.site.LH <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.LH)
write.csv(env.site.LH, "LH_site_env_data.csv", row.names = FALSE)

# Re-order env data to pop afs
LH.data <- cbind(unique(NEpops[,2]),env.site.LH)
LH.data <- LH.data[order(LH.data[,1]),]


LH_var <- cbind(LH.data$bio_05, LH.data$bio_19)
colnames(LH_var) <- c("bio_05", "bio_19")
rownames(LH_var) <- pop.data$popnum
LH_var <- as.data.frame(LH_var)



# rcp45_2070
env.site.rcp45_2070 <- raster::extract(rcp45_2070.rasters,sites)
env.site.rcp45_2070 <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.rcp45_2070)
write.csv(env.site.rcp45_2070, "rcp45_2070_site_env_data.csv", row.names = FALSE)

# Re-order env data to pop afs
rcp45_2070.data <- cbind(unique(NEpops[,2]),env.site.rcp45_2070)
rcp45_2070.data <- rcp45_2070.data[order(rcp45_2070.data[,1]),]


rcp45_2070_var <- cbind(rcp45_2070.data$bio_05, rcp45_2070.data$bio_19)
colnames(rcp45_2070_var) <- c("bio_05", "bio_19")
rownames(rcp45_2070_var) <- pop.data$popnum
rcp45_2070_var <- as.data.frame(rcp45_2070_var)




# rcp85_2070
env.site.rcp85_2070 <- raster::extract(rcp85_2070.rasters,sites)
env.site.rcp85_2070 <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.rcp85_2070)
write.csv(env.site.rcp85_2070, "rcp85_2070_site_env_data.csv", row.names = FALSE)

# Re-order env data to pop afs
rcp85_2070.data <- cbind(unique(NEpops[,2]),env.site.rcp85_2070)
rcp85_2070.data <- rcp85_2070.data[order(rcp85_2070.data[,1]),]


rcp85_2070_var <- cbind(rcp85_2070.data$bio_05, rcp85_2070.data$bio_19)
colnames(rcp85_2070_var) <- c("bio_05", "bio_19")
rownames(rcp85_2070_var) <- pop.data$popnum
rcp85_2070_var <- as.data.frame(rcp85_2070_var)



#####################################################################################
######## Run AlleleShift genomic vulnerability models for each time period ##########

# Calibrate the models using current environmental data
genp.count.model <- count.model(genp.candidates, env.data=current_var, ordistep=F, cca.model=F)

genp.pred.baseline <- count.pred(genp.count.model, env.data=current_var)

genp.freq.model <- freq.model(genp.pred.baseline)
genp.freq.baseline <- freq.pred(genp.freq.model, count.predicted=genp.pred.baseline)


# 4. Check how well the models predict baseline allele frequencies

# get AlleleShift model fit per population
np <- 38
lm.Rsq <- data.frame(Pop = genp.freq.baseline[1:np, "Pop"], Pop.label = as.character(c(1:np)), 
                     N = genp.freq.baseline[1:np, "N"]/2, GAM.rsq = numeric(np))
res <- for (p1 in 1:np) {
  Pop.focal <- genp.freq.baseline[p1, "Pop"]
  long.p <- genp.freq.baseline[genp.freq.baseline$Pop == Pop.focal, 
  ]
  m <- stats::lm(Freq.e2 ~ Allele.freq, data = long.p)
  lm.Rsq[p1, "GAM.rsq"] <- round(summary(m)$r.squared, 2)
}

model.fit <- lm.Rsq[order(lm.Rsq$GAM.rsq, decreasing = TRUE), ]
model.fit
write.csv(model.fit, "GV211_calibration_model_fit.csv")

# identify list of pops for which the model performs poorly (R2<0.5) to omit from the final analyses
retained.pops <- droplevels(as.data.frame(model.fit$Pop[model.fit$GAM.rsq>0.5]))
colnames(retained.pops) <- "pop"
retained.pops$popNumber <- sprintf("%02d", as.numeric(gsub("[[:alpha:]]", "", as.character(retained.pops$pop))))


# 5. Predict past/future allele frequencies

# late holocene
#predict allele counts
genp.pred.LH <- count.pred(genp.count.model, env.data=LH_var)
head(genp.pred.LH)
# transform predicted allele counts to allele frequencies and compare with base allele frequencies
genp.freq.LH <- freq.pred(genp.freq.model, count.predicted=genp.pred.LH)
head(genp.freq.LH)

# reduce to retained pops
reduced.genp.freq.LH <- genp.freq.LH[genp.freq.LH$Pop %in% retained.pops[,1],]


# rcp45_2070
genp.pred.rcp45_2070 <- count.pred(genp.count.model, env.data=rcp45_2070_var)
head(genp.pred.rcp45_2070)

genp.freq.rcp45_2070 <- freq.pred(genp.freq.model, count.predicted=genp.pred.rcp45_2070)
head(genp.freq.rcp45_2070)

# reduce to retained pops
reduced.genp.freq.rcp45_2070 <- genp.freq.rcp45_2070[genp.freq.rcp45_2070$Pop %in% retained.pops[,1],]


# rcp85_2070
genp.pred.rcp85_2070 <- count.pred(genp.count.model, env.data=rcp85_2070_var)
head(genp.pred.rcp85_2070)

genp.freq.rcp85_2070 <- freq.pred(genp.freq.model, count.predicted=genp.pred.rcp85_2070)
head(genp.freq.rcp85_2070)

# reduce to retained pops
reduced.genp.freq.rcp85_2070 <- genp.freq.rcp85_2070[genp.freq.rcp85_2070$Pop %in% retained.pops[,1],]
head(reduced.genp.freq.rcp85_2070)

# subset coords to retained pops
retained.coords <- pop.data[ ,1:3]
colnames(retained.coords) <- c("pop", "x", "y")
retained.coords <- retained.coords[retained.coords$pop %in% retained.pops[,1],]


# early to late holocene
#predict allele counts
genp.pred.LH <- count.pred(genp.count.model, env.data=LH_var)
head(genp.pred.LH)
# transform predicted allele counts to allele frequencies and compare with base allele frequencies
genp.freq.LH <- freq.pred(genp.freq.model, count.predicted=genp.pred.LH)
head(genp.freq.LH)

# reduce to retained pops
reduced.genp.freq.LH <- genp.freq.LH[genp.freq.LH$Pop %in% retained.pops[,1],]



genp.pred.MH <- count.pred(genp.count.model, env.data=MH_var)
head(genp.pred.MH)
# transform predicted allele counts to allele frequencies and compare with base allele frequencies
genp.freq.MH <- freq.pred(genp.freq.model, count.predicted=genp.pred.MH)
head(genp.freq.MH)

# reduce to retained pops
reduced.genp.freq.MH <- genp.freq.MH[genp.freq.MH$Pop %in% retained.pops[,1],]



genp.pred.EH <- count.pred(genp.count.model, env.data=EH_var)
head(genp.pred.EH)
# transform predicted allele counts to allele frequencies and compare with base allele frequencies
genp.freq.EH <- freq.pred(genp.freq.model, count.predicted=genp.pred.EH)
head(genp.freq.EH)

# reduce to retained pops
reduced.genp.freq.EH <- genp.freq.EH[genp.freq.EH$Pop %in% retained.pops[,1],]

# summarise number of maladapted alleles fixed per population for all pops
af.cand <- makefreq(genp.candidates, missing="0")

# get adaptive allele freqs
adapt.allele.freq <- af.cand[, colnames(af.cand) %in% unique(reduced.genp.freq.rcp85_2070$Allele)]
# count number of loci where adaptive allele is missing from pop
res3 <- as.data.frame(apply(adapt.allele.freq, 1, function(x) sum(x == 0)))
colnames(res3) <- "adapt.allele.missing"
barplot(res3$adapt.allele.missing)

rownames(res3) %in% hybrid.pops[,1]
pure.hybrid <- as.data.frame(ifelse(rownames(res3) %in% hybrid.pops[,1], "Hybrid", "Pure"))
pure.hybrid.species <- cbind(res3, species.pop, pure.hybrid)
colnames(pure.hybrid.species) <- c("adapt.allele.missing", "species", "hybrid")
write.csv(pure.hybrid.species, "pure_hybrid_species.csv")

# mean number of missing alleles
adapt.allele.missing <- aggregate(pure.hybrid.species$adapt.allele.missing ~ species + hybrid, data=pure.hybrid.species, mean)
colnames(adapt.allele.missing) <- c("species", "hybrid", "adapt.allele.missing")
adapt.allele.missing$percent <- adapt.allele.missing$adapt.allele.missing/211*100
adapt.allele.missing <- adapt.allele.missing[order(adapt.allele.missing$hybrid, decreasing = T), ]
adapt.allele.missing <- adapt.allele.missing[order(adapt.allele.missing$species),]
write.csv(adapt.allele.missing, "adapt_allele_missing.csv")


# plot candidate alleles missing hybrid vs pure for each species
NE_colors <- ggpubfigs::friendly_pal("ito_seven")[c(2,5,1,3,4)]
scales::show_col(NE_colors)
names(NE_colors) = levels(NE_species)



bar_colors.agg <- c(rep(NE_colors[1], sum(adapt.allele.missing$species == "eachamensis")),
                    rep(NE_colors[2], sum(adapt.allele.missing$species == "Malanda")),
                    rep(NE_colors[3], sum(adapt.allele.missing$species == "splendida")),
                    rep(NE_colors[4], sum(adapt.allele.missing$species == "Tully")),
                    rep(NE_colors[5], sum(adapt.allele.missing$species == "utcheensis")))



pdf(file = "adapt_alleles_absent_percent.pdf", height = 6, width = 6)
par(mar = c(7, 4, 2, 2) + 0.2, bg = "#fdf6e3", fg="#073642")
x <- barplot(adapt.allele.missing$percent, col=bar_colors.agg,
             ylim=c(0,30), width = 0.9,
             border="#fdf6e3", xaxt="n",
             cex.names = 0.8,
             cex.axis = 1,
             beside=TRUE, xpd=T,axes = F,
             legend=FALSE,
             cex.lab=0.8)
axis(1, at=x, labels = adapt.allele.missing$hybrid, cex.axis = 1, tick=F, line = -1, col="#073642", col.lab="#073642", col.axis="#073642")
axis(2, lwd = 2, cex.axis = 1.2, at = seq(0, 50, 10), las=1, col="#073642", col.lab="#073642", col.axis="#073642")
title(ylab=expression("Percentage of adaptive alleles absent"), line=2, cex.lab=1.2, col="#073642", col.lab="#073642", col.axis="#073642")

dev.off()



# format df
reduced.popNumber <-  as.character(gsub("[[:alpha:]]", "", as.character(levels(droplevels(as.factor(reduced.genp.freq.LH$Pop))))))
reduced.pop <-  as.character(levels(droplevels(as.factor(reduced.genp.freq.LH$Pop))))
reduced.species <- as.character(gsub("[[:digit:]]+", "", levels(droplevels(as.factor(reduced.pop)))))
reduced.pure.hybrid.species <- pure.hybrid.species[rownames(pure.hybrid.species) %in% reduced.pop, ]


# get median allele frequency changes between time periods per site 
deltaAF.EH <- get.deltaAF(reduced.genp.freq.EH)
deltaAF.MH <- get.deltaAF(reduced.genp.freq.MH)
deltaAF.LH <- get.deltaAF(reduced.genp.freq.LH)
deltaAF.rcp45_2070 <- get.deltaAF(reduced.genp.freq.rcp45_2070)
deltaAF.rcp85_2070 <- get.deltaAF(reduced.genp.freq.rcp85_2070)
temporal.deltaAF <- as.data.frame(cbind(deltaAF.EH$deltaAF, deltaAF.MH$deltaAF, deltaAF.LH$deltaAF,deltaAF.rcp45_2070$deltaAF,deltaAF.rcp85_2070$deltaAF))
colnames(temporal.deltaAF) <- c("EarlyHolocene_Current", "MidHolocene_Current", "LateHolocene_Current", "Current_rcp45_2070", "Current_rcp85_2070")
temporal.deltaAF$pop <- reduced.pop
temporal.deltaAF$popNumber <-  reduced.popNumber
temporal.deltaAF$species <-  reduced.species

# subset mean_admix df to same pops as reduced model
reducded.mean_admix <- mean_admix[mean_admix$site %in% reduced.pop,]
temporal.deltaAF$splendida_ancestry <- reducded.mean_admix$splendida


#set colours and labels for plotting
species.pop <- as.data.frame(pop.data[,4])
colnames(species.pop) <- "species"
unique_species <- as.data.frame(sort(unique(as.factor(species.pop$species))))
NE_species <- factor(unique_species[,1])

NE_colors <- ggpubfigs::friendly_pal("ito_seven")[c(2,5,1,3,4)]
scales::show_col(NE_colors)
names(NE_colors) = levels(NE_species)

# map colours to points
ind_cols <- mapvalues(pop.data$`NEXY$species`, names(NE_colors), unique(NE_colors))
names(ind_cols) = pop.data[,1]
reduced.ind_cols <- ind_cols[names(ind_cols) %in% reduced.pop]
temporal.deltaAF$colour <- reduced.ind_cols

NE_colors <- temporal.deltaAF$colour

# absolute allele frequency change
abs.af <- abs(temporal.deltaAF$Current_rcp85_2070)
#abs.af <- temporal.deltaAF$Current_rcp85_2070

# get lm stats
mod <- lm(abs.af~temporal.deltaAF$splendida_ancestry)
modsum <- summary(mod)
r2 = modsum$adj.r.squared

leg.labs <- c(substitute(paste(italic("M. eachamensis"))),
              substitute(paste("Malanda rainbowfish")),
              substitute(paste(italic("M. splendida"))),
              substitute(paste(italic("M. utcheensis"))))


temporal.deltaAF$hybrid <- reduced.pure.hybrid.species$hybrid
# pop18 is eacham/malanda hybrid, so add manuall
temporal.deltaAF$hybrid[18] <- "Hybrid"

temporal.deltaAF <- temporal.deltaAF %>% mutate(pch = case_when(stri_detect_fixed(hybrid, "Hybrid") ~ 4,
                                                                stri_detect_fixed(hybrid, "Pure") ~ 1))

write.table(temporal.deltaAF, "temporal_deltaAF.csv", sep = ",", row.names = FALSE)


pdf("deltaAFvsMsplendida_ancestry.pdf", width = 6, height = 6)
par(bg = "#fdf6e3", fg="#073642")
plot(temporal.deltaAF$splendida_ancestry, abs.af, type = "n", xlab="", ylab="", cex.axis=1.2, col="#073642", col.lab="#073642", col.axis="#073642")
box(lwd=2)
points(temporal.deltaAF$splendida_ancestry, abs.af, col = NE_colors, pch=temporal.deltaAF$pch, cex=2, lwd = 3)
title(xlab="splendida ancestry", line=2, cex.lab=1.2, col="#073642", col.lab="#073642", col.axis="#073642")
title(ylab=expression(paste(Delta," allele frequency")), line=2, cex.lab=1.2, col="#073642", col.lab="#073642", col.axis="#073642")
abline(mod, col="#073642", lwd=2)
dev.off()

