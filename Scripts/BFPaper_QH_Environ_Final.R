#*******************************************#
#Pulling Bioclim data and pop coordinates together
#*******************************************#

library(stats); library(dismo); library(maptools); library(rgdal)

#For PCA and correlations

#install_github("vqv/ggbiplot")
library(tibble);library(scales);library(ggbiplot);library(ggplot2)
#install_github("vqv/ggbiplot") # must be installed after devtools is loaded in
library("FactoMineR")
library("factoextra")
library(ggrepel)


#*******************************************#
#Read in data and building matrix
#*******************************************#

setwd("~/Desktop/Quercus/QH_Environ/BF_Paper")
# http://worldclim.org/version2 is where data comes from @ 2.5 min

library(raster)
bioclim2.5 <- getData(name = "worldclim", var = "bio", res = 2.5)

#Extract data for QH Pops
QHloc_Final <- read.csv("BailiePaperCoords.csv", sep=",", header=T, stringsAsFactors = FALSE)
data(wrld_simpl)
climateFinal <- extract(bioclim2.5, QHloc_Final[,1:2]) #Extract climate data for locations

#Renaming columns and rows 
bioclim_names <- c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
colnames(climateFinal) <- bioclim_names
rownames(climateFinal) <-  QHloc_Final$Pops



#*******************************************#

#PCA

#*******************************************#

Region <- QHloc_Final$Region
Collection <- QHloc_Final$Key
Collection_Region <- QHloc_Final$Collection_Region

climateFinal2 <- as.data.frame(climateFinal, row.names = TRUE)

climateFinal.pca2 <- prcomp(climateFinal2, center = TRUE,scale. = TRUE)
climateFinal.pca2  #loadings

#EnvironLoadings <- climateFinal.pca2$rotation
#write.csv(EnvironLoadings, 'EnvironLoadings.csv')
#res.pca <- PCA(climateFinal, graph = FALSE)
#var <- get_pca_var(res.pca)
#write.csv(var$contrib, file = "contributions_Final.csv")
#fviz_pca_var(res.pca, col.var = "black")




#*******************************************#
  
  # Testing for Correlation
  
#*******************************************#
library("Hmisc")

#The output of the function rcorr() is a list containing the following elements : 
#- r : the correlation matrix 
#- n : the matrix of the number of observations used in analyzing each pair of variables 
#- P : the p-values corresponding to the significance levels of correlations.
res <- cor(climateFinal)
round(res,3)
cor(climateFinal, use = "complete.obs")

res2 <- rcorr(as.matrix(climateFinal))
res2
# Extract the correlation coefficients
res2$r
# Extract p-values
res2$P



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

QHcorrmat_varsremoved <- flattenCorrMatrix(res2$r, res2$P)
write.csv(QHcorrmat_varsremoved, file = "corrmat.csv")

symbols_varsremoved <- symnum(res, abbr.colnames = FALSE)
write.csv(symbols_varsremoved, file = "corrmatsymbols_varsremoved.csv")
#symbol codes for above .csv are as follows:  
# 0 = ‘ ’ 
# 0.3 = ‘.’ 
# 0.6 ‘,’ 
# 0.8 ‘+’ 
# 0.9 ‘*’ 
# 0.95 ‘B’
# 1 '1'

#Here are the results: 
#BIO1 = Annual mean temperature	Correlated R2 > 0.9 with	BIO11, BIO10, BIO6
#BIO2 = Annual mean diurnal range	Correlated R2 > 0.9 with	NA
#BIO3 = Isothermality	Correlated R2 > 0.9 with	NA
#BIO4 = Temperature seasonality	Correlated R2 > 0.9 with	BIO7
#BIO5 = Max temperature of warmest month	Correlated R2 > 0.9 with NA
#BIO6 = Min temperature of coldest month	Correlated R2 > 0.9 with	BIO1, BIO11
#BIO7 = Annual temperature range	Correlated R2 > 0.9 with	BIO4
#BIO8 = Mean temperature of wettest quarter	Correlated R2 > 0.9 with	NA
#BIO9 = Mean temperature of driest quarter	Correlated R2 > 0.9 with	NA
#BIO10 = Mean temperature of warmest quarter	Correlated R2 > 0.9 with	BIO1
#BIO11 = Mean temperature of coldest quarter	Correlated R2 > 0.9 with	BIO1, BIO6
#BIO12 = Annual precipitation	Correlated R2 > 0.9 with	BIO13, BIO16
#BIO13 = Precipitation of wettest month	Correlated R2 > 0.9 with	BIO12, BIO16, BIO18
#BIO14 = Precipitation of driest month	Correlated R2 > 0.9 with	BIO17
#BIO15 = Precipitation seasonality	Correlated R2 > 0.9 with	NA
#BIO16 = Precipitation of wettest quarter	Correlated R2 > 0.9 with	BIO12, BIO13, BIO18
#BIO17 = Precipitation of driest quarter	Correlated R2 > 0.9 with	BIO14
#BIO18 = Precipitation of warmest quarter	Correlated R2 > 0.9 with	BIO13, BIO16
#BIO19 = Precipitation of coldest quarter	Correlated R2 > 0.9 with	NA


# Using the info for correlation above, the following variables were removed from the analysis:
# BIO1, BIO6, BIO7, BIO13 , BIO16, BIO17

#Remove correlated variables
climateFinal2 <- as.data.frame(climateFinal2, row.names = TRUE)
climateFinal2 <- climateFinal2[,-c(1,6,7,13,16,17)]
climateFinal.pca2 <- prcomp(climateFinal2, center = TRUE,scale. = TRUE)
climateFinal.pca2  #loadings


#Adding column for regions
climateFinal2 <- cbind(Collection,Region,climateFinal2)
colnames(climateFinal2) <- c("Key", "Region", "BIO2","BIO3","BIO4","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO14","BIO15","BIO17","BIO19")
dim(climateFinal2)






#Run this first to get % variation explained on x and y axes for next plot


p <- ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1, var.scale = 1.5, 
         var.axes = FALSE, labels= NULL, groups= NULL)+
  geom_point(aes(color=Collection, size=0.2),alpha = 1)+
  guides(size=FALSE)+
  geom_density2d(alpha=.5)+
  scale_colour_manual(values = c("blue", "darkgrey","blue", "darkgrey"))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    legend.position = c(0.73, 0.2),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6),
    legend.key.size = unit(1,"line"),
    legend.key = element_rect(fill = "white", color = "white"),
    legend.box.background = element_rect(color="black", size=1),
    legend.text = element_text(size = 12),
    legend.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black"))

p +  guides(fill = guide_legend(override.aes = list(linetype = 0, shape= 16)),
            colour = guide_legend(override.aes = list(linetype=c(0,0),
                                                      shape=c(16,16),
                                                      size=c(3.5))))
#Run Arrow_edits_BFPaper.R for black BIO arrows
jitter <- position_jitter(width = 0.25, height = 0.25)


p2 <- ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1, var.scale = 1.5, 
              var.axes = TRUE, labels= NULL, groups= NULL, alpha = 0)+
  geom_point(position = jitter,aes(color=Collection_Region,shape=Region, size = 1),alpha = 0.85)+
  guides(size=FALSE, shape=FALSE)+
  scale_colour_manual(values = c("blue", "blue", "red", "red"))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black"))

p2$layers <- c(p2$layers, p2$layers[[1]])
p2$layers <- c(p2$layers, p2$layers[[3]])

p2 + theme(
  legend.position = c(0.73, 0.2),
  legend.justification = c("left", "top"),
  legend.box.just = "left",
  legend.margin = margin(6, 6, 6, 6),
  legend.key.size = unit(1,"line"),
  legend.key = element_rect(fill = "white", color = "white"),
  legend.box.background = element_rect(color="black", size=1),
  legend.text = element_text(size = 12),
  legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape= 16)),
         colour = guide_legend(override.aes = list(linetype=c(0,0),
                                                   shape=c(16,17,16,17),
                                                   size=c(4))))


p3 <- ggbiplot(climateFinal.pca2, ellipse = FALSE, obs.scale = 1, var.scale = 1.5, 
               var.axes = FALSE, labels= NULL, groups= NULL)+
  geom_point(aes(color=Collection_Region, size=0.2, shape=Region),alpha = 1)+
  geom_density2d(alpha=.5)+
  guides(size=FALSE, shape=FALSE)+
  scale_colour_manual(values = c("blue", "blue", "darkgrey", "darkgrey"))+
  theme_light(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(fill="transparent", colour = "black"))


p3 + theme(
  legend.position = c(0.1, .95),
  legend.justification = c("left", "top"),
  legend.box.just = "left",
  legend.margin = margin(6, 6, 6, 6),
  legend.key.size = unit(1,"line"),
  legend.key = element_rect(fill = "white", color = "white"),
  legend.box.background = element_rect(color="black", size=1),
  legend.text = element_text(size = 12),
  legend.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape= 16)),
         colour = guide_legend(override.aes = list(linetype=c(0,0),
                                                   shape=c(16,17,16,17),
                                                   size=c(4))))

#Run script Arrow_edits.R to change color of arrows for this graph.

