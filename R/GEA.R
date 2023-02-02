# RDA analysis for Natural hybridisation reduces vulnerability to climate change
library(adegenet)
library(raster)
library(vegan)
library(LEA)
library(packfor)
library(psych)
library(memgene)
library(ggpubfigs)

# define utility functions
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

impute.data <- function(genind) {
  
  alleles <- genind@tab
  snps <- alleles[,seq(1,ncol(alleles),2)]
  colnames(snps) <- locNames(genind)
  
  
  md <- round((sum(is.na(snps)))/(dim(snps)[1]*dim(snps)[2])*100, digits = 2)
  cat(paste0("Total missing data = "),md,"%\n")
  
  write.lfmm(snps, "NER.lfmm")
  lfmm.obj = read.lfmm("NER.lfmm")
  
  project.snmf = snmf("NER.lfmm", K = 5, 
                      entropy = TRUE, repetitions = 10, seed = 81853598, CPU = 4,
                      project = "new")
  
  # select the run with the lowest cross-entropy value
  best = which.min(cross.entropy(project.snmf, K = 5))
  
  # Impute the missing genotypes
  impute(project.snmf, "NER.lfmm", method = 'mode', K = 5, run = best)
  
  imputed.snps <- read.lfmm("NER.lfmm_imputed.lfmm")
  
  
  dim(imputed.snps)
  imputed.snps[1:10,1:10]
  colnames(imputed.snps) <- locNames(genind)
  rownames(imputed.snps) <- indNames(genind)
  # check total % missing data
  (sum(is.na(imputed.snps)))/(dim(imputed.snps)[1]*dim(imputed.snps)[2])*100
  
  
  imputed.snps2 <- 2-imputed.snps
  
  #create indexes for the desired order
  x <- order(c(1:ncol(imputed.snps), 1:ncol(imputed.snps2)))
  
  #cbind d1 and d2, interleaving columns with x
  dat <- cbind(imputed.snps, imputed.snps2)[,x]
  colnames(dat) <- colnames(genind@tab)
  
  
  geno <- genind
  geno@tab <- dat
  mode(geno@tab) <- "integer"
  
  return(geno)
}

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

#####################################################################################
############################ get formatted raster stacks ############################
setwd("/PATH/TO/data/environmental_rasters/")

for (i in c("EH", "MH", "LH", "current", "rcp45_2070", "rcp85_2070")) {
  files <- list.files(path = i, all.files = TRUE, recursive=F, full.names=T, pattern = "\\.grd$")
  rasters <- stack(c(files))
  assign(paste(i,".rasters", sep=""),rasters)
}

####################################################################################
############################ extract environmental data ############################
setwd("/PATH/TO/data")

# get coords and extract env data
NEXY <- read.table("NE_GEA_XY.txt", header = TRUE)
coords <- data.frame(x=NEXY$x,y=NEXY$y)
sites <- SpatialPoints(coords, proj4string = current.rasters@crs)
env.site.current <- raster::extract(current.rasters,sites)
env.site.current <- cbind.data.frame(coordinates(sites),NEXY$species, env.site.current)

#######################################################################################
################################### get genotypes ##########################

NE_LD300 <- read.structure("NE342LD300.stru", n.ind=342, n.loc=13734, onerowperind=FALSE, col.lab=1, col.pop = 1, 
                           col.others=FALSE, row.marknames=1, NA.char="-9", ask=TRUE, quiet=FALSE)
NEpops <- read.table("NE342_popmap.txt", header = FALSE)
NE_LD300$pop <- as.factor(NEpops[,2])
NE_LD300



#######################################################################################
################################### Run analyses ##########################
setwd("/PATH/TO/WD")

# check and impute missing data using snmf ancestry proportions
geno <- impute.data(NE_LD300)

# sort alleles based on splendida major/minor allele (this ensures consistent "direction" of genomic offset relative to splendida)
splendida <- c("splendida01", "splendida02", "splendida11", "splendida12", "splendida17", "splendida35")
geno <- sort_alleles.pop(geno, splendida)


# get pop allele frequencies
genp <- genind2genpop(geno)
af <- makefreq(genp, missing="0")

# get minor allele frequencies
NEaf <- af[,seq(2,ncol(af),2)]
colnames(NEaf) <- locNames(geno)

# Re-order env data to match pop afs
pop.data <- cbind(unique(NEpops[,2]),env.site.current)
pop.data$popnum <- sprintf("%02d", 1:nrow(pop.data))
pop.data <- pop.data[order(pop.data[,1]),]


# use pop numbers as rownames for allele freq matrix for easier plotting
NEaf.popNumbers <- NEaf
rownames(NEaf.popNumbers) <- pop.data$popnum


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

# run vif forward selection procedure to remove variables with the highest VIF one at a time until all remaining variables are below 10
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



# reduced RDA model using the retained environmental variables conditioned on the admixture qmatrix
NEpop.RDA <- vegan::rda(NEaf.popNumbers ~ bio_05 + bio_19 + Condition(MEM.sel), data = retained.env)
NEpop.RDA


# calculate significance of the reduced model
mod_perm <- anova.cca(NEpop.RDA, nperm=1000, model="reduced", parallel = 4)

#set colours and labels for plotting
species.pop <- as.data.frame(pop.data[,4])
colnames(species.pop) <- "species"
unique_species <- as.data.frame(sort(unique(as.factor(species.pop$species))))
NE_species <- factor(unique_species[,1])

NE_colors <- friendly_pal("ito_seven")[c(2,5,1,3,4)]
names(NE_colors) = levels(NE_species)


leg.labs.RDA <- c(expression(italic("M. eachamensis")),
                  expression("Malanda rainbowfish"),
                  expression(italic("M. splendida")),
                  expression("Tully rainbowfish"),
                  expression(italic("M. utcheensis")))




x.lab <- paste0("RDA1 (", paste(round((NEpop.RDA$CA$eig[1]/NEpop.RDA$tot.chi*100),2)),"%)")
y.lab <- paste0("RDA2 (", paste(round((NEpop.RDA$CA$eig[2]/NEpop.RDA$tot.chi*100),2)),"%)")

# plot points
pdf(file = "RDA.pdf", height = 6, width = 6)
par(bg = "#fdf6e3", fg="#073642", mgp=c(1.5,0.5,0), cex.lab=0.8)
pRDAplot <- plot(NEpop.RDA, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.axis=0.8, col="#073642", col.lab="#073642", col.axis="#073642")
text(NEpop.RDA, "bp", choices = c(1, 2), labels = c("bio_05", "bio_19"), col="#073642", cex=0.8)
with(species.pop, text(NEpop.RDA, display = "sites", col = NE_colors[species], pch = 16, cex=1, bg= "transparent"))
with(species.pop, legend("topright", legend = leg.labs.RDA, col=unique(NE_colors), pch=16, pt.cex=1.5, cex=0.8, xpd=1, box.lty = 0, bg= "transparent"))
dev.off()
