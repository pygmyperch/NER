# NER GEA and GV MS
library(gridExtra)
library(vegan)
library(cowplot)
library(gridGraphics)
library(plyr)
library(dplyr)
library(adegenet)
library(raster)
library(packfor)
library(psych)
library(AlleleShift)
library(ggplot2)
library(ggpubr)
library(memgene)
library(tidyverse)


# get genotypes
setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/Input_file/MS/MS_files")
NE_LD300 <- read.structure("NE342LD300.stru", n.ind=342, n.loc=13734, onerowperind=FALSE, col.lab=1, col.pop = 1, col.others=FALSE, row.marknames=1, NA.char="-9", ask=TRUE, quiet=FALSE)
NEpops <- read.table("NE342_popmap.txt", header = FALSE)
NE_LD300$pop <- as.factor(NEpops[,2])
NE_LD300


# impute missing data
geno <- NE_LD300
# get allele counts
alleles <- geno@tab
alleles[1:10,1:10]


# get genotypes (counts of reference allele) and clean up locus names
snps <- alleles[,seq(1,ncol(alleles),2)]
colnames(snps) <- locNames(geno)
snps[1:20,1:20]

# check total % missing data
(sum(is.na(snps)))/(dim(snps)[1]*dim(snps)[2])*100

# impute missing data with most common genotype
snps <- apply(snps, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

# check % missing data again
(sum(is.na(snps)))/(dim(snps)[1]*dim(snps)[2])*100
snps[1:10,1:10]




#####################################################################################
############################ get formatted raster stacks ############################
setwd("/Users/chrisbrauer/Projects/narrow_endemics/MS/niche_models/biomod/CHELSA_data/NER_formatted_data/")

for (i in c("EH", "MH", "LH", "current", "rcp45_2070", "rcp85_2070")) {
  files <- list.files(path = i, all.files = TRUE, recursive=F, full.names=T, pattern = "\\.grd$")
  rasters <- stack(c(files))
  assign(paste(i,".rasters", sep=""),rasters)
}

#####################################################################################
#####################################################################################

setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/GEA/GEA_GV_MS")
save.image(file='pop_GEA_MS.RData')
load(file='pop_GEA_MS.RData')


# get coords and extract env data
NEXY <- read.table("NE_GEA_XY.txt", header = TRUE)
coords <- data.frame(x=NEXY$x,y=NEXY$y)
sites <- SpatialPoints(coords, proj4string = current.rasters@crs)
env.site.current <- raster::extract(current.rasters,sites)
env.site.current <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.current)
write.csv(env.site.current, "site_current_env_data.csv", row.names = FALSE)



################################### impute missing genotypes ##########################
alleles <- NE_LD300@tab

# get genotypes (counts of reference allele) and clean up locus names
snps <- alleles[,seq(1,ncol(alleles),2)]
colnames(snps) <- locNames(NE_LD300)
# check total % missing data
(sum(is.na(snps)))/(dim(snps)[1]*dim(snps)[2])*100


#create and read .lfmm file
write.lfmm(snps, "NER.lfmm")
lfmm.obj = read.lfmm("NER.lfmm")


project.snmf = snmf("NER.lfmm", K = 5, 
                    entropy = TRUE, repetitions = 10,
                    project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.snmf, K = 5))

# Impute the missing genotypes
impute(project.snmf, "NER.lfmm", method = 'mode', K = 5, run = best)

imputed.snps <- read.lfmm("NER.lfmm_imputed.lfmm")


dim(imputed.snps)
imputed.snps[1:10,1:10]
colnames(imputed.snps) <- locNames(NE_LD300)
rownames(imputed.snps) <- indNames(NE_LD300)
# check total % missing data
(sum(is.na(imputed.snps)))/(dim(imputed.snps)[1]*dim(imputed.snps)[2])*100


imputed.snps2 <- 2-imputed.snps

#create indexes for the desired order
x <- order(c(1:ncol(imputed.snps), 1:ncol(imputed.snps2)))

#cbind d1 and d2, interleaving columns with x
dat <- cbind(imputed.snps, imputed.snps2)[,x]
colnames(dat) <- colnames(NE_LD300@tab)


geno <- NE_LD300
geno@tab <- dat
mode(geno@tab) <- "integer"

save.image(file='pop_GEA_MS.RData')

# sort alleles based on splendida major/minor allele (this ensures consistent "direction" of genomic offset among populations)
splendida <- c("splendida01", "splendida02", "splendida11", "splendida12", "splendida17", "splendida35")

sort_alleles.pop <- function(genind, pops) {
  
  tab <- as.data.frame(genind@tab)
  tab$pop <- genind@pop
  tab <- tab[tab$pop %in% pops,]
  tab <- within(tab, rm(pop))
  
  res <- as.data.frame(colSums(tab, na.rm = T))
  colnames(res) <- "allele_count"
  res$allele_name <- rownames(res)
  
  
  my.ord <- do.call(c, lapply(seq(2, nrow(res), by = 2), (function(i){
    c(i-1, i)[order(res$allele_count[(i-1):i], decreasing = T)]})))
  
  
  res1 <- data.frame(res[my.ord,])
  
  genind@tab <- genind@tab[, res1$allele_name]
  
  return(genind)
  
}

geno <- sort_alleles.pop(geno, splendida)


# get pop allele frequencies
genp <- genind2genpop(geno)
af <- makefreq(genp, missing="0")

# get minor allele frequencies
NEaf <- af[,seq(2,ncol(af),2)]
colnames(NEaf) <- locNames(geno)
NEaf[1:10,1:10]


# Re-order env data to pop afs
pop.data <- cbind(unique(NEpops[,2]),env.site.current)
pop.data$popnum <- sprintf("%02d", 1:nrow(pop.data))
pop.data <- pop.data[order(pop.data[,1]),]


# use pop numbers as rownames for allele freq matrix for better plotting
NEaf.popNumbers <- NEaf
rownames(NEaf.popNumbers) <- pop.data$popnum




#############get admixture proportions and define hybrid and pure inds
#####
# get site averages and reorder as per genetic data
admix <- read.table("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/admixture/supervised/NE342MS_LD300_bed.5.Q")
admix <- cbind(NEpops[,2],admix)
colnames(admix) <- c("site", "splendida", "Malanda", "eachamensis", "utcheensis", "Tully")


# calculate mean q values per site
mean_admix <- admix %>% group_by(site) %>% summarize(splendida = mean(splendida), Malanda = mean(Malanda), 
                                                     eachamensis = mean(eachamensis), utcheensis = mean(utcheensis), Tully = mean(Tully))

# check rows match (doesn't work with pop number AF matrix)
stopifnot(all(rownames(NEaf) == mean_admix[,1]))
stopifnot(all(rownames(NEaf) == as.character(pop.data[,1])))
stopifnot(all(mean_admix[,1] == as.character(pop.data[,1])))

meanQmatrix <- as.matrix(mean_admix[,2:6])

popIds <- as.data.frame(colnames(meanQmatrix)[apply(meanQmatrix,1,which.max)])
popIds <- as.data.frame(apply(meanQmatrix,1,which.max))


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





# RDA #
#######################################
##### RDA
# Variable selection and spatial filtering
# PCoA of SNP data
snps.bray <- vegdist(NEaf.popNumbers, method="bray")
snp.pcoa <- cmdscale(snps.bray, k=nrow(NEaf.popNumbers)-1, eig=T, add=T)
plot(snp.pcoa$points[ ,1], snp.pcoa$points[ ,2])
eig <- snp.pcoa$eig/sum(snp.pcoa$eig)
bst <- unname(bstick(length(eig)))

axes <- scores(snp.pcoa)
only <- min(which((eig>bst) == FALSE))
y <- axes[,c(1:only-1)]




# Forward selection of hydroclimatic variables
env_var <- cbind(pop.data$bio_05, pop.data$bio_06, pop.data$bio_18, pop.data$bio_19)
colnames(env_var) <- c("bio_05", "bio_06", "bio_18", "bio_19")

env_var <- as.data.frame(scale(env_var))

# run procedure to remove variables with the highest VIF one at a time until all remaining variables are below 10
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}
keep.env <-vif_func(in_frame=env_var,thresh=10,trace=T)
keep.env  # the retained environmental PCs # "bio_05" "bio_06" "bio_19"

env_var <- env_var[, (colnames(env_var) %in% keep.env)]



# Calculate the full RDA & get adj R2
env.full <- rda(y, env_var)
env.R2full <- RsquareAdj(env.full)$adj.r.squared


# Forward selection using full adj R2 as threshold
sel <- forward.sel(y, as.matrix(env_var), adjR2thresh=env.R2full)
numENV <- nrow(sel)
sigENV <- sort(sel[,2])
ENV.sel <- env_var[,c(sigENV)]
colnames(ENV.sel) <- sel[,1]


envcor <- corr.test(ENV.sel[,c(1:2)], ENV.sel[,c(1:2)])
envcor <- envcor$r

# run procedure to remove variables with the highest VIF one at a time until all remaining variables are below 10
keep.env <-vif_func(in_frame=ENV.sel,thresh=10,trace=T)
keep.env  # the retained environmental PCs

retained.env <- as.data.frame(ENV.sel)
xy <- as.matrix(NEXY[,4:5])


# get MEM vectors to control for spatial genetic structure
snps.bray <- vegdist(NEaf.popNumbers, method="bray")
snpDM <- as.matrix(snps.bray)
memAnalysis <- mgQuick(snpDM, xy, forwardPerm = 10000)
MEM.sel <- memAnalysis$memSelected
colnames(MEM.sel) <- c("MEM1", "MEM3", "MEM4")

# plot spatial representation of retained MEMs
mgMap(xy, memAnalysis$memSelected[, 1:2])
dev.copy2pdf(file = "MEM12.pdf", height = 12, width = 6)
dev.off()
mgMap(xy, memAnalysis$memSelected[, c(1,3)])
dev.copy2pdf(file = "MEM34.pdf", height = 12, width = 6)
dev.off()


# reduced RDA model using the retained environmental variables conditioned on the admixture qmatrix
NEpop.RDA <- vegan::rda(NEaf.popNumbers ~ bio_05 + bio_19 + Condition(MEM.sel), data = retained.env)
NEpop.RDA

# how much genetic variation can be explained by our environmental model?
RsquareAdj(NEpop.RDA)$r.squared

# how much inertia is associated with each axis
screeplot(NEpop.RDA)

# calculate significance of the reduced model, marginal effect of each term and significance of each axis
#   (this can take several hours, depending on your computer)
mod_perm <- anova.cca(NEpop.RDA, nperm=1000, model="reduced", parallel = 4) #test significance of the model
#mod_perm
#margin_perm <- anova.cca(NEpop.RDA, by="margin", nperm=1000,  parallel = 4)#test marginal effect of each individual term in the model
#axis_perm <- anova.cca(NEpop.RDA, by="axis", nperm=1000, parallel = 4)#test significance of each constrained axis


#set colours and labels for plotting
species.pop <- as.data.frame(pop.data[,4])
colnames(species.pop) <- "species"
unique_species <- as.data.frame(sort(unique(as.factor(species.pop$species))))
NE_species <- factor(unique_species[,1])

NE_colors <- ggpubfigs::friendly_pal("ito_seven")[c(2,5,1,3,4)]
scales::show_col(NE_colors)
names(NE_colors) = levels(NE_species)


leg.labs.RDA <- c(expression(italic("M. eachamensis")),
                  expression("Malanda rainbowfish"),
                  expression(italic("M. splendida")),
                  expression("Tully rainbowfish"),
                  expression(italic("M. utcheensis")))




x.lab <- paste0("RDA1 (", paste(round((NEpop.RDA$CA$eig[1]/NEpop.RDA$tot.chi*100),2)),"%)")
y.lab <- paste0("RDA2 (", paste(round((NEpop.RDA$CA$eig[2]/NEpop.RDA$tot.chi*100),2)),"%)")

# plot points
pdf(file = "/Users/chrisbrauer/GoogleDrive/CurrentManuscripts/NER/NCC/reviews/second_revisions/figures/NE_RDA_popMS.pdf", height = 6, width = 6)
par(bg = "#fdf6e3", fg="#073642", mgp=c(1.5,0.5,0), cex.lab=0.8)
pRDAplot <- plot(NEpop.RDA, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.axis=0.8, col="#073642", col.lab="#073642", col.axis="#073642")
text(NEpop.RDA, "bp", choices = c(1, 2), labels = c("bio_05", "bio_19"), col="#073642", cex=0.8)
with(species.pop, text(NEpop.RDA, display = "sites", col = NE_colors[species], pch = 16, cex=1, bg= "transparent"))
with(species.pop, legend("topright", legend = leg.labs.RDA, col=unique(NE_colors), pch=16, pt.cex=1.5, cex=0.8, xpd=1, box.lty = 0, bg= "transparent"))
dev.off()



## obtain list of candidate loci
# calculate coordinates of loci more than 3SD from mean locus scores for each significant RDA axis (in this case RDA1 and RDA2)
locus_scores <- cbind(seq.int(1, length(NEaf.popNumbers[1,]), 1), NEpop.RDA$CCA$v[,1:2])

lims1 <- as.data.frame(mean(locus_scores[,2]) + c(-1, 1) * 3 * sd(locus_scores[,2]))
RDA1.candidates <- locus_scores[ which( locus_scores[,2] < lims1[1,] | locus_scores[,2] > lims1[2,]) , ]

lims2 <- as.data.frame(mean(locus_scores[,3]) + c(-1, 1) * 3 * sd(locus_scores[,3]))
RDA2.candidates <- locus_scores[ which( locus_scores[,3] < lims2[1,] | locus_scores[,3] > lims2[2,]) , ]


# generate final list of unique candidate loci
RDA1_2.candidates <- rbind(RDA1.candidates, RDA2.candidates)
RDA1_2.candidates <- as.data.frame(RDA1_2.candidates[!duplicated(RDA1_2.candidates), ])
candidates <- as.matrix(RDA1_2.candidates[,1])
write.csv(candidates, "NEpopRDA_candidate_SNPsMS211.csv", row.names = FALSE)


# get candidate locus names
candidates <- as.data.frame(unique(genp@loc.fac))[candidates, ]
#subset AF
geno.candidates <- geno[loc=candidates]

genp.candidates <- genind2genpop(geno.candidates)


# get candidate data
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(locus_scores[,2],3)
cand2 <- outliers(locus_scores[,3],3)


#get total number of candidates
ncand <- length(cand1) + length(cand2)
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)


# calculate correlation between candidate snps and environmental variables
snp <- as.data.frame(NEaf.popNumbers)
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp[,nam]
  ENV.sel[i,] <- apply(retained.env,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,ENV.sel)  
head(cand)

length(cand$snp[duplicated(cand$snp)])
env_mat <- cbind(cand$axis, duplicated(cand$snp)) 
table(env_mat[env_mat[,1]==1,2])# none on axis 1
table(env_mat[env_mat[,1]==2,2])# none on axis 2
cand <- cand[!duplicated(cand$snp),]

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,6] <- names(which.max(abs(bar[4:5]))) # gives the variable
  cand[i,7] <- max(abs(bar[4:5]))              # gives the correlation
}

colnames(cand)[6] <- "best_predictor"
colnames(cand)[7] <- "correlation"

table(cand$best_predictor) 

write.csv(cand, "candidate_env_correlations.csv")


setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/GEA/GEA_GV_MS/")
save.image(file='pop_GEA_MS.RData')
#load(file='pop_GEA_MS.RData')




# GV #
#########
# map colours to points
ind_cols <- mapvalues(pop.data$`NEXY$species`, names(NE_colors), unique(NE_colors))


# current
current_var <- cbind(pop.data$bio_05, pop.data$bio_19)
colnames(current_var) <- c("bio_05", "bio_19")
rownames(current_var) <- pop.data$popnum
current_var <- as.data.frame(current_var)


# modified plotting function
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


plotEH.current <- population.shift_cb(EH_var, current_var, title = "EarlyHolocene >> Current", from = "EarlyHolocene", to = "Current", option="PCA")
plotEH.current



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


plotMH.current <- population.shift_cb(MH_var, current_var, title = "MidHolocene >> Current", from = "MidHolocene", to = "Current", option="PCA")
plotMH.current


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




plotEH.MH <- population.shift_cb(EH_var, MH_var, title = "", from = "EarlyHolocene", to = "Current", option="PCA")
plotEH.MH

plotMH.LH <- population.shift_cb(MH_var, LH_var, title = "", from = "EarlyHolocene", to = "Current", option="PCA")
plotMH.LH

plotLH.current <- population.shift_cb(LH_var, current_var, title = "", from = "LateHolocene", to = "Current", option="PCA")
plotLH.current

plotrcp45_2070.current <- population.shift_cb(current_var, rcp45_2070_var, title = "", from = "Current", to = "rcp45_2070", option="PCA")
plotrcp45_2070.current

plotrcp85_2070.current <- population.shift_cb(current_var, rcp85_2070_var, title = "", from = "Current", to = "rcp85_2070", option="PCA")
plotrcp85_2070.current

pdf(file = "EH.MH.pdf", height = 6, width = 6)
plotEH.MH
dev.off()

pdf(file = "MH.LH.pdf", height = 6, width = 6)
plotMH.LH
dev.off()

pdf(file = "LH.current.pdf", height = 6, width = 6)
plotLH.current
dev.off()


pdf(file = "plotrcp45_2070.current.pdf", height = 6, width = 6)
plotrcp45_2070.current
dev.off()

pdf(file = "plotrcp85_2070.current.pdf", height = 6, width = 6)
plotrcp85_2070.current
dev.off()



# rda based
plotLH.current.rda <- population.shift_cb(LH_var, current_var, title = "", from = "LateHolocene", to = "Current", option="RDA")
plotrcp45_2070.current.rda <- population.shift_cb(current_var, rcp45_2070_var, title = "", from = "Current", to = "rcp45_2070", option="RDA")
plotrcp85_2070.current.rda <- population.shift_cb(current_var, rcp85_2070_var, title = "", from = "Current", to = "rcp85_2070", option="RDA")

plotLH.current.rda
plotrcp45_2070.current.rda
plotrcp85_2070.current.rda


# save plots
current.climate.pca <- rda(current_var, scale = TRUE)
x.lab <- paste0("PC1 (", paste(round((current.climate.pca$CA$eig[1]/current.climate.pca$tot.chi*100),2)),"%)")
y.lab <- paste0("PC2 (", paste(round((current.climate.pca$CA$eig[2]/current.climate.pca$tot.chi*100),2)),"%)")

leg.labs <- c(substitute(paste(italic("M. eachamensis"))),
              substitute(paste("Malanda rainbowfish")),
              substitute(paste(italic("M. splendida"))),
              substitute(paste("Tully rainbowfish")),
              substitute(paste(italic("M. utcheensis"))))



# plot points
pdf(file = "current_climate_pop_hulls.pdf", height = 6, width = 6)
pRDAplot <- plot(current.climate.pca, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1)
fit <- envfit(current.climate.pca, current_var, perm = 999)
plot(fit, cex = 0.6)
with(species.pop, text(current.climate.pca, display = "sites", col = NE_colors[species], pch = 16, cex=0.6, bg= "transparent"))
#ordihull(current.climate.pca, group = species.pop$species, col = NE_colors, lwd = 1, lty = 5)
ordihull(current.climate.pca, group = species.pop$species, draw = "polygon", col = NE_colors, alpha = 0.3, border = NA)
#with(species.pop, legend("topleft", legend = leg.labs, col=unique(NE_colors), pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent"))
inset.x <- -0.005
inset.y <- -0.03

for (i in 1:5) {
  # start inset y at 0.02 and place each successive legend element 0.03 below
  inset.y <- inset.y+0.03
  legend("topleft", inset = c(inset.x,inset.y), legend = leg.labs[[i]], col=unique(NE_colors)[i], pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent")
}

dev.off()

LH.climate.pca <- rda(LH_var, scale = TRUE)
x.lab <- paste0("PC1 (", paste(round((LH.climate.pca$CA$eig[1]/LH.climate.pca$tot.chi*100),2)),"%)")
y.lab <- paste0("PC2 (", paste(round((LH.climate.pca$CA$eig[2]/LH.climate.pca$tot.chi*100),2)),"%)")

# plot points
pdf(file = "LH_climate_pop.pdf", height = 6, width = 6)
pRDAplot <- plot(LH.climate.pca, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1, main = "LH climate")
fit <- envfit(LH.climate.pca, LH_var, perm = 999)
plot(fit, cex = 0.6)
with(species.pop, text(current.climate.pca, display = "sites", col = NE_colors[species], pch = 16, cex=0.6, bg= "transparent"))
with(species.pop, legend("topleft", legend = names(NE_colors), col=unique(NE_colors), pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent"))
dev.off()




rcp45_2070.climate.pca <- rda(rcp45_2070_var, scale = TRUE)
x.lab <- paste0("PC1 (", paste(round((rcp45_2070.climate.pca$CA$eig[1]/rcp45_2070.climate.pca$tot.chi*100),2)),"%)")
y.lab <- paste0("PC2 (", paste(round((rcp45_2070.climate.pca$CA$eig[2]/rcp45_2070.climate.pca$tot.chi*100),2)),"%)")

# plot points
pdf(file = "rcp45_2070_climate_pop.pdf", height = 6, width = 6)
pRDAplot <- plot(rcp45_2070.climate.pca, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1, main = "rcp45_2070 climate")
fit <- envfit(rcp45_2070.climate.pca, rcp45_2070_var, perm = 999)
plot(fit, cex = 0.6)
with(species.pop, text(rcp45_2070.climate.pca, display = "sites", col = NE_colors[species], pch = 16, cex=0.6, bg= "transparent"))
with(species.pop, legend("topleft", legend = names(NE_colors), col=unique(NE_colors), pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent"))
dev.off()


rcp85_2070.climate.pca <- rda(rcp85_2070_var, scale = TRUE)
x.lab <- paste0("PC1 (", paste(round((rcp85_2070.climate.pca$CA$eig[1]/rcp85_2070.climate.pca$tot.chi*100),2)),"%)")
y.lab <- paste0("PC2 (", paste(round((rcp85_2070.climate.pca$CA$eig[2]/rcp85_2070.climate.pca$tot.chi*100),2)),"%)")

# plot points
pdf(file = "rcp85_2070_climate_pop.pdf", height = 6, width = 6)
pRDAplot <- plot(rcp85_2070.climate.pca, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1, main = "rcp85_2070 climate")
fit <- envfit(rcp85_2070.climate.pca, rcp85_2070_var, perm = 999)
plot(fit, cex = 0.6)
with(species.pop, text(rcp85_2070.climate.pca, display = "sites", col = NE_colors[species], pch = 16, cex=0.6, bg= "transparent"))
with(species.pop, legend("topleft", legend = names(NE_colors), col=unique(NE_colors), pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent"))
dev.off()



setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/GEA/GEA_GV_MS")
save.image(file='pop_GEA_MS.RData')
load(file='pop_GEA_MS.RData')



# Calibrate the models current to 
genp.count.model <- count.model(genp.candidates, env.data=current_var, ordistep=F, cca.model=F)

genp.pred.baseline <- count.pred(genp.count.model, env.data=current_var)
head(genp.pred.baseline)

genp.freq.model <- freq.model(genp.pred.baseline)
genp.freq.baseline <- freq.pred(genp.freq.model, count.predicted=genp.pred.baseline)
head(genp.freq.baseline)

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
# and remove splendida
retained.pops <- droplevels(as.data.frame(model.fit$Pop[model.fit$GAM.rsq>0.5]))
colnames(retained.pops) <- "pop"
retained.pops$popNumber <- sprintf("%02d", as.numeric(gsub("[[:alpha:]]", "", as.character(retained.pops$pop))))
#retained.pops <- retained.pops[!grepl("splendida",retained.pops$pop),]


#retained.pops <- as.data.frame(model.fit$Pop)
#colnames(retained.pops) <- "pop"
#retained.pops$popNumber <- sprintf("%02d", as.numeric(gsub("[[:alpha:]]", "", as.character(retained.pops$pop))))


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



setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/GEA/GEA_GV_MS")
save.image(file='pop_GEA_MS.RData')
load(file='pop_GEA_MS.RData')


# summarise number of maladapted alleles fixed per population for all pops
af.cand <- makefreq(genp.candidates, missing="0")

# get adaptive allele freqs
adapt.allele.freq <- af.cand[, colnames(af.cand) %in% unique(reduced.genp.freq.rcp85_2070$Allele)]
# count number of loci where adaptive allele is missing from pop
res3 <- as.data.frame(apply(adapt.allele.freq, 1, function(x) sum(x == 0)))
colnames(res3) <- "adapt.allele.missing"
barplot(res3$adapt.allele.missing)

# bar_colors <- c(rep("limegreen", sum(species.pop$species == "eachamensis")),
#                 rep("darkorange1", sum(species.pop$species == "Malanda")),
#                 rep("dodgerblue3", sum(species.pop$species == "splendida")),
#                 rep("deeppink", sum(species.pop$species == "Tully")),
#                 rep("red", sum(species.pop$species == "utcheensis")))          
# 
# hybrids <- rownames(res3) %in% hybrid.pops[,1]
# hybrids <- ifelse(rownames(res3) %in% hybrid.pops[,1], "*", "")
# 
# par(mar = c(7, 4, 2, 2) + 0.2)
# x <- barplot(res3$adapt.allele.missing/211*100, col=bar_colors,
#              ylim=c(0,50), width = 0.9,
#              border="white", xaxt="n",
#              cex.names = 0.6,
#              cex.axis = 0.6,
#              beside=TRUE, xpd=T,axes = F,
#              legend=FALSE,
#              cex.lab=0.8)
# #axis(1, at=x, labels = F)
# axis(2, cex.axis = 0.6, at = seq(0, 50, 10), las=1)
# #rotate 60 degrees (srt = 60)
# text(cex=0.7, x=x-.25, y=-0.5, paste(rownames(res3)), xpd=TRUE, srt=45, adj = c(1,1.8))
# text(x,(res3$adapt.allele.missing/211*100)+2,labels=hybrids, cex=2)
# dev.off()
# 
# 

res3$adapt.allele.missing

rownames(res3) %in% hybrid.pops[,1]
pure.hybrid <- as.data.frame(ifelse(rownames(res3) %in% hybrid.pops[,1], "Hybrid", "Pure"))
pure.hybrid.species <- cbind(res3, species.pop, pure.hybrid)
colnames(pure.hybrid.species) <- c("adapt.allele.missing", "species", "hybrid")
write.csv(pure.hybrid.species, "pure_hybrid_species.csv")

adapt.allele.missing <- aggregate(pure.hybrid.species$adapt.allele.missing ~ species + hybrid, data=pure.hybrid.species, mean)
colnames(adapt.allele.missing) <- c("species", "hybrid", "adapt.allele.missing")
adapt.allele.missing <- adapt.allele.missing[order(adapt.allele.missing$hybrid, decreasing = T), ]
adapt.allele.missing <- adapt.allele.missing[order(adapt.allele.missing$species),]
write.csv(adapt.allele.missing, "adapt_allele_missing.csv")



#NE_colors <- c("limegreen", "darkorange1", "dodgerblue3", "deeppink", "red")
NE_colors <- c("#2aa198", "#cb4b16", "#268bd2", "#d33682", "#6c71c4")
cols <- c("#cb4b16", "#2aa198", "#6c71c4", "#268bd2")
taxa <- c("Malanda", "M.eachamensis", "M. utcheensis", "M.splendida")



bar_colors.agg <- c(rep(NE_colors[1], sum(adapt.allele.missing$species == "eachamensis")),
                    rep(NE_colors[2], sum(adapt.allele.missing$species == "Malanda")),
                    rep(NE_colors[3], sum(adapt.allele.missing$species == "splendida")),
                    rep(NE_colors[4], sum(adapt.allele.missing$species == "Tully")),
                    rep(NE_colors[5], sum(adapt.allele.missing$species == "utcheensis")))


taxa <- c("Malanda", "M. eachamensis", "M. utcheensis", "M.splendida")



leg.labs <- c(substitute(paste("Malanda")),
              substitute(paste(italic("M. eachamensis"))),
              substitute(paste(italic("M. utcheensis"))),
              substitute(paste(italic("M. splendida"))))



pdf(file = "adapt_alleles_absent_percent.pdf", height = 6, width = 6)
par(mar = c(7, 4, 2, 2) + 0.2, bg = "#fdf6e3", fg="#073642")
x <- barplot(adapt.allele.missing$adapt.allele.missing/251*100, col=bar_colors.agg,
             ylim=c(0,30), width = 0.9,
             border="#fdf6e3", xaxt="n",
             cex.names = 0.8,
             cex.axis = 1,
             beside=TRUE, xpd=T,axes = F,
             legend=FALSE,
             cex.lab=0.8)
axis(1, at=x, labels = adapt.allele.missing$hybrid, cex.axis = 1, tick=F, line = -1, col="#073642", col.lab="#073642", col.axis="#073642")
axis(2, lwd = 2, cex.axis = 1.2, at = seq(0, 50, 10), las=1, col="#073642", col.lab="#073642", col.axis="#073642")
title(ylab=expression("% adaptive alleles absent"), line=2, cex.lab=1.2, col="#073642", col.lab="#073642", col.axis="#073642")
#with(species.pop, legend("bottom", legend = names(NE_colors), col=unique(NE_colors), pch=16, pt.cex=1, cex=0.7, xpd=1, box.lty = 0, bg= "transparent"))

# for (i in 1:4) {
#   # start inset y at 0.15 and place each successive legend element 0.03 to the right
#   inset.x <- inset.x+0.15
#   
#   legend("bottomleft", xpd = TRUE, horiz = TRUE, inset = c(inset.x,inset.y), legend = leg.labs[[i]], col=cols[i], 
#          bty = "n", pch = 15, cex = 0.6, pt.cex = 1)
# }

dev.off()



setwd("/Users/chrisbrauer/GoogleDrive/Narrow_endemicRF/bioinfo/genome_call/GEA/GEA_GV_MS")
save.image(file='pop_GEA_MS.RData')
load(file='pop_GEA_MS.RData')
