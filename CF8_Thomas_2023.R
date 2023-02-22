#Code to calculate biomarker indices, correlations, and make plots for Thomas et al., in prep
#Early Holocene Laurentide Ice Sheet retreat influenced summer atmospheric circulation in Baffin Bay

library(tidyr)
library(lipdR) #to read and interact with LiPD data #remotes::install_github("nickmckay/lipdR")
library(geoChronR) #for plotting mostly
library(magrittr); library(FactoMineR); library(factoextra) #we'll be using the magrittr pipe ( %>% ) for simplicity
library(dplyr) #and dplyr for data.frame manipulation
library(ggplot2); library(ggtern) #for plotting
remotes::install_github("nickmckay/compositeR")
library(compositeR) #remotes::install_github("nickmckay/compositeR")
library(foreach) #for parallel processing
library(doParallel)#for parallel processing
library(purrr); library(RColorBrewer)
library(matrixStats); library(cowplot)

#########
#--------------Age model shown in Fig. S3
#If need to make age model (but this is already done and saved in the Thomas.CF8.2023 lipd file)
setwd("/Users/Elizabeth/Documents/MacBacon_2.2")
source("Bacon.R")

C=readLipd()

#Make an age model for the CF8 core using modeled ages of the sand-lacustrine sediment contact from two other Holocene CF8 cores (02CF8 and CF817)
#Tried with thick = 2, but resulted in unrealistically low sed rates just above the boundary at 78.5 cm, so using thick = 1
C = runBacon(C,bacon.thick = 1,d.min = 0, bacon.acc.mean = c(50,20),boundary= 78.5)
C <- loadBaconOutput(C,site.name="Thomas.CF8.2021",bacon.dir = "/Users/Elizabeth/Documents/MacBacon_2.2/Cores")
C <- mapAgeEnsembleToPaleoData(C)

writeLipd(C,"/Users/Elizabeth/Documents/R/GeoChronR/data")


################
#-----------Figure S5
#GDGT concentrations
ages <- selectData(C,var.name = "agemedian")$values 
brGDGT.Ia <- selectData(C,var.name = "peakarea1022brgdgtia")$values
brGDGT.Ib <- selectData(C,var.name = "peakarea1020brgdgtib")$values
brGDGT.Ic <- selectData(C,var.name = "peakarea1018brgdgtic")$values
brGDGT.IIa.5Me <- selectData(C,var.name = "peakarea1036brgdgtiia5me")$values
brGDGT.IIa.6Me <- selectData(C,var.name = "peakarea1036'brgdgtiia6me")$values
brGDGT.IIb.5Me <- selectData(C,var.name = "peakarea1034brgdgtiib5me")$values
brGDGT.IIb.6Me <- selectData(C,var.name = "peakarea1034'brgdgtiib6me")$values
brGDGT.IIc.5Me <- selectData(C,var.name = "peakarea1032brgdgtiic5me")$values
brGDGT.IIc.6Me <- selectData(C,var.name = "peakarea1032'brgdgtiic6me")$values
brGDGT.IIIa.5Me <- selectData(C,var.name = "peakarea1050brgdgtiiia5me")$values
brGDGT.IIIa.6Me <- selectData(C,var.name = "peakarea1050'brgdgtiiia6me")$values
brGDGT.IIIb.5Me <- selectData(C,var.name = "peakarea1048brgdgtiiib5me")$values
brGDGT.IIIb.6Me <- selectData(C,var.name = "peakarea1048'brgdgtiiib6me")$values
brGDGT.IIIc.5Me <- selectData(C,var.name = "peakarea1046brgdgtiiic5me")$values
brGDGT.IIIc.6Me <- selectData(C,var.name = "peakarea1046'brgdgtiiic6me")$values
isoGDGT.0 <- selectData(C, var.name = "peakarea1302gdgt0")$values
isoGDGT.1 <- selectData(C, var.name = "peakarea1300gdgt1")$values
isoGDGT.2 <- selectData(C, var.name = "peakarea1298gdgt2")$values
isoGDGT.3 <- selectData(C, var.name = "peakarea1296gdgt3")$values
isoGDGT.4 <- selectData(C, var.name = "peakarea1292gdgt4")$values
isoGDGT.4isomer <- selectData(C, var.name = "peakarea1292'gdgt4'")$values
IntStd <- selectData(C, var.name = "peakarea744internalstandard")$values
Std_ug <- selectData(C, var.name = "massc46standardadded")$values
SedMass_g <- selectData(C, var.name = "massdrysediment")$values

#to save data so you don't have to upload and choose from the LiPD file again if you restart R
brGDGTdf <- data.frame(IntStd,Std_ug,SedMass_g,brGDGT.Ia,brGDGT.Ib, brGDGT.Ic, brGDGT.IIa.5Me,brGDGT.IIa.6Me,brGDGT.IIb.5Me,brGDGT.IIb.6Me,brGDGT.IIc.5Me,brGDGT.IIc.6Me, brGDGT.IIIa.5Me,brGDGT.IIIa.6Me, brGDGT.IIIb.5Me,brGDGT.IIIb.6Me,brGDGT.IIIc.5Me,brGDGT.IIIc.6Me,ages)
#write.csv(brGDGTdf,"/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_brGDGT_inputdata.csv")

isoGDGTdf <- data.frame(IntStd,Std_ug,SedMass_g,isoGDGT.0,isoGDGT.1,isoGDGT.2,isoGDGT.3,isoGDGT.4,isoGDGT.4isomer,ages)
#write.csv(isoGDGTdf,"/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_isoGDGT_inputdata.csv")

#----start again once have the csv saved
#brGDGTdf <- read.csv("/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_GDGT_inputdata.csv")
#brGDGTdf <-brGDGTdf[,c(2:22)]
#brGDGTdf <- data.frame(ages,IntStd,Std_ug,SedMass_g,isoGDGT.0,isoGDGT.4,brGDGT.Ia,brGDGT.Ib, brGDGT.Ic, brGDGT.IIa.5Me,brGDGT.IIa.6Me,brGDGT.IIb.5Me,brGDGT.IIb.6Me,brGDGT.IIc.5Me,brGDGT.IIc.6Me, brGDGT.IIIa.5Me,brGDGT.IIIa.6Me, brGDGT.IIIb.5Me,brGDGT.IIIb.6Me,brGDGT.IIIc.5Me,brGDGT.IIIc.6Me)
brGDGTdf <- brGDGTdf[rowSums(brGDGTdf[,c(4:18)], na.rm = TRUE)>0,]
ages <- brGDGTdf$ages
IntStd <- brGDGTdf$IntStd
Std_ug <- brGDGTdf$Std_ug
SedMass_g <- brGDGTdf$SedMass_g
brGDGTdf <- brGDGTdf[,c(6:20)]

brGDGTconc <- ((brGDGTdf/IntStd)* Std_ug)/SedMass_g
brGDGT.Total <- rowSums(brGDGTconc, na.rm = TRUE)
brGDGT.Total[brGDGT.Total==0] <- NA
brGDGTfa <- brGDGTdf/rowSums(brGDGTdf, na.rm = TRUE)
brGDGTfa$ages <- ages

brGDGTfaMeth <- data.frame(fIa.meth = brGDGTdf$brGDGT.Ia/rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE),
                           fIIa.5Me.meth = brGDGTdf$brGDGT.IIa.5Me/rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE),
                           fIIIa.5Me.meth = brGDGTdf$brGDGT.IIIa.5Me/rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE),
                           fIIa.6Me.meth = brGDGTdf$brGDGT.IIa.6Me/rowSums(brGDGTdf[,c("brGDGT.IIa.6Me","brGDGT.IIIa.6Me")], na.rm = TRUE),
                           fIIIa.6Me.meth = brGDGTdf$brGDGT.IIIa.6Me/rowSums(brGDGTdf[,c("brGDGT.IIa.6Me","brGDGT.IIIa.6Me")], na.rm = TRUE),
                           fIb.meth = brGDGTdf$brGDGT.Ib/rowSums(brGDGTdf[,c("brGDGT.Ib","brGDGT.IIb.5Me","brGDGT.IIIb.5Me")], na.rm = TRUE),
                           fIIb.5Me.meth = brGDGTdf$brGDGT.IIb.5Me/rowSums(brGDGTdf[,c("brGDGT.Ib","brGDGT.IIb.5Me","brGDGT.IIIb.5Me")], na.rm = TRUE),
                           fIIIb.5Me.meth = brGDGTdf$brGDGT.IIIb.5Me/rowSums(brGDGTdf[,c("brGDGT.Ib","brGDGT.IIb.5Me","brGDGT.IIIb.5Me")], na.rm = TRUE),
                           fIIb.6Me.meth = brGDGTdf$brGDGT.IIb.6Me/rowSums(brGDGTdf[,c("brGDGT.IIb.6Me","brGDGT.IIIb.6Me")], na.rm = TRUE),
                           fIIIb.6Me.meth = brGDGTdf$brGDGT.IIIb.6Me/rowSums(brGDGTdf[,c("brGDGT.IIb.6Me","brGDGT.IIIb.6Me")], na.rm = TRUE),
                           fIc.meth = brGDGTdf$brGDGT.Ic/rowSums(brGDGTdf[,c("brGDGT.Ic","brGDGT.IIc.5Me","brGDGT.IIIc.5Me")], na.rm = TRUE),
                           fIIc.5Me.meth = brGDGTdf$brGDGT.IIc.5Me/rowSums(brGDGTdf[,c("brGDGT.Ic","brGDGT.IIc.5Me","brGDGT.IIIc.5Me")], na.rm = TRUE),
                           fIIIc.5Me.meth = brGDGTdf$brGDGT.IIIc.5Me/rowSums(brGDGTdf[,c("brGDGT.Ic","brGDGT.IIc.5Me","brGDGT.IIIc.5Me")], na.rm = TRUE),
                           fIIc.6Me.meth = brGDGTdf$brGDGT.IIc.6Me/rowSums(brGDGTdf[,c("brGDGT.IIc.6Me","brGDGT.IIIc.6Me")], na.rm = TRUE),
                           fIIIc.6Me.meth = brGDGTdf$brGDGT.IIIc.6Me/rowSums(brGDGTdf[,c("brGDGT.IIc.6Me","brGDGT.IIIc.6Me")], na.rm = TRUE)
)
brGDGTfaMeth[is.na(brGDGTfaMeth)] <- 0

#boxplot of all brGDGTs
par(mar=c(0.5, 1.5, 0.1, 1), mfrow=c(1,1),
    oma = c(2, 4, 2, 4))
boxplot(brGDGTfa[,c(1:15)])

#time series of all brGDGTs
plot(ages[!is.na(brGDGT.Total)],brGDGT.Total[!is.na(brGDGT.Total)], type = "l")

########
#----start again once have the csv saved
#isoGDGTdf <- read.csv("/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_isoGDGT_inputdata.csv")
#isoGDGTdf <-isoGDGTdf[,c(2:22)]
isoGDGTdf <- isoGDGTdf[rowSums(isoGDGTdf[,c(4:9)], na.rm = TRUE)>0,]
ages <- isoGDGTdf$ages
IntStd <- isoGDGTdf$IntStd
Std_ug <- isoGDGTdf$Std_ug
SedMass_g <- isoGDGTdf$SedMass_g
isoGDGTdf <- isoGDGTdf[,c(4:9)]

isoGDGTconc <- ((isoGDGTdf/IntStd)* Std_ug)/SedMass_g
isoGDGT.Total <- rowSums(isoGDGTconc, na.rm = TRUE)
isoGDGT.Total[isoGDGT.Total==0] <- NA
isoGDGTfa <- isoGDGTdf/rowSums(isoGDGTdf, na.rm = TRUE)
isoGDGTfa$ages <- ages

#boxplot of all isoGDGTs
par(mar=c(0.5, 1.5, 0.1, 1), mfrow=c(1,1),
    oma = c(2, 4, 2, 4))
boxplot(isoGDGTfa[,c(1:6)])

#time series of all isoGDGTs
plot(ages[!is.na(isoGDGT.Total)],isoGDGT.Total[!is.na(isoGDGT.Total)], type = "l")


#-------------stats about the values of fractional abundances-----------------

brGDGT_faStats <- data.frame(matrix(nrow=15, ncol=0))

brGDGT_faStats$value <- colnames(brGDGTfa[1:15])
brGDGT_faStats$median <- apply(brGDGTfa[1:15],2,median,na.rm = TRUE)
brGDGT_faStats$stdev <- apply(brGDGTfa[1:15],2,sd,na.rm = TRUE)
brGDGT_faStats$max <- apply(brGDGTfa[1:15],2,max,na.rm = TRUE)
brGDGT_faStats$min <- apply(brGDGTfa[1:15],2,min,na.rm = TRUE)

write.csv(brGDGT_faStats,"/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/CF8_Biomarkerplots_2022/CF8brGDGT_faStats.csv")

################
#-----------Figure S7
#Time series of GDGT indices and Principal Components Analysis

# --------------------- Calculate GDGT indices -----------------------------------------

GDGT0.4 <- isoGDGT.0/isoGDGT.4

GDGT0.4[is.na(GDGT0.4)] <- 0
GDGTindices<-data.frame(GDGT0.4)

GDGTindices$MBT.5me <- rowSums(brGDGTfa[,c("brGDGT.Ia","brGDGT.Ib","brGDGT.Ic")],na.rm = TRUE)/
  (rowSums(brGDGTfa[,c("brGDGT.Ia","brGDGT.Ib","brGDGT.Ic","brGDGT.IIa.5Me","brGDGT.IIb.5Me",
                       "brGDGT.IIc.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE))

GDGTindices$FA_6Me <- rowSums(brGDGTfa[,c("brGDGT.IIa.6Me","brGDGT.IIb.6Me",
                                          "brGDGT.IIc.6Me","brGDGT.IIIa.6Me","brGDGT.IIIb.6Me","brGDGT.IIIc.6Me")], na.rm = TRUE)

GDGTindices$IR <- rowSums(brGDGTfa[,c("brGDGT.IIa.6Me","brGDGT.IIb.6Me","brGDGT.IIc.6Me","brGDGT.IIIa.6Me","brGDGT.IIIb.6Me","brGDGT.IIIc.6Me")],na.rm = TRUE)/
  (rowSums(brGDGTfa[,c("brGDGT.IIa.6Me","brGDGT.IIb.6Me","brGDGT.IIc.6Me","brGDGT.IIIa.6Me","brGDGT.IIIb.6Me","brGDGT.IIIc.6Me",
                       "brGDGT.IIa.5Me","brGDGT.IIb.5Me","brGDGT.IIc.5Me","brGDGT.IIIa.5Me","brGDGT.IIIb.5Me","brGDGT.IIIc.5Me")], na.rm = TRUE))

GDGTindices$BIT <- rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me")],na.rm = TRUE)/
  (rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE)+(isoGDGT.4))

GDGTindices$cren_conc <- ((isoGDGT.0/IntStd)* Std_ug)/SedMass_g

GDGTindices$brGDGTconc <- rowSums(brGDGTconc[],na.rm = TRUE)

GDGTindices$HP5 <- (brGDGTdf[,c("brGDGT.IIIa.5Me")])/
  (rowSums(brGDGTdf[,c("brGDGT.IIa.5Me","brGDGT.IIIa.5Me")], na.rm = TRUE))

GDGTindices$CBT5Me <- -log(rowSums(brGDGTdf[,c("brGDGT.Ib","brGDGT.IIb.5Me")],na.rm = TRUE)/
                             ((rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me")], na.rm = TRUE))))

GDGTindices$fC <- ((rowSums(brGDGTdf[,c("brGDGT.Ib","brGDGT.IIb.5Me","brGDGT.IIIb.5Me","brGDGT.IIb.6Me","brGDGT.IIIb.6Me")],na.rm = TRUE))
                   +2*(rowSums(brGDGTdf[,c("brGDGT.IIc.5Me","brGDGT.IIIc.5Me","brGDGT.IIc.6Me","brGDGT.IIIc.6Me")],na.rm = TRUE)))/
  ((rowSums(brGDGTdf[,c("brGDGT.Ia","brGDGT.IIa.5Me","brGDGT.IIIa.5Me","brGDGT.IIa.6Me","brGDGT.IIIa.6Me","brGDGT.Ib","brGDGT.IIb.5Me","brGDGT.IIIb.5Me","brGDGT.IIb.6Me","brGDGT.IIIb.6Me")],na.rm = TRUE))
   +2*(rowSums(brGDGTdf[,c("brGDGT.IIc.5Me","brGDGT.IIIc.5Me","brGDGT.IIc.6Me","brGDGT.IIIc.6Me")],na.rm = TRUE)))

GDGTindices[GDGTindices == 0] <- NA


# --------------------- PCA of brGDGTs -----------------------------------------

# Scale brGDGT concentrations such that most abundant compound in each sample = 1
brGDGTscaled <- brGDGTdf

for (i in 1:nrow(brGDGTscaled)) {
  max_val <- max(brGDGTscaled[i,], na.rm = T)
  brGDGTscaled[i,] <- brGDGTdf[i,] / max_val
}

brGDGTscaled[is.na(brGDGTscaled)] <- 0

res.pca <- PCA(brGDGTscaled, graph = TRUE)
eig.val <- get_eigenvalue(res.pca)
var <- get_pca_var(res.pca)
ind <- get_pca_ind(res.pca)

pca.res <- data.frame(ages, res.pca$ind$coord[,1],res.pca$ind$coord[,2])


## varimax rotation
#info from here: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/varimax
#and from here: https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
varimaxLoadings <- promax(res.pca$var$coord)
varimaxScores <- scale(res.pca$ind$coord) %*% varimaxLoadings$rotmat

scoresScaled <- varimaxScores
for (i in 1:ncol(scoresScaled)) {
  max_val <- max(scoresScaled[,i], na.rm = T)
  scoresScaled[,i] <- varimaxScores[,i] / max_val
}

CF8varimaxPCloadings <- data.frame(varimaxLoadings$loadings[,1:5])
write.csv(CF8varimaxPCloadings,"/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_varimaxPCLoadings.csv")

#Plot time series of PCA & GDGT indices

midpt1 = mean(c(max(GDGT0.4[!is.na(GDGTindices$GDGT0.4)]),min(GDGT0.4[!is.na(GDGTindices$GDGT0.4)])))
GDGT0.4<-GDGTindices$GDGT0.4

hexa <- ggplot() + geom_line(aes(x = ages, y = brGDGT.TPH$hexa)) + geom_point(aes(x = ages,  y = brGDGT.TPH$hexa, fill = GDGT0.4), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800),
                                                                                                                                 position = 'top') +
  xlab('Age (cal yr BP)') + ylab('Fraction Hexa') +  theme_classic()+theme(legend.position = "none")

MBT <- ggplot() + geom_line(aes(x = ages, y = GDGTindices$MBT.5me)) + geom_point(aes(x = ages,  y = GDGTindices$MBT.5me, fill = GDGT0.4), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  xlab(NULL) + ylab('MBT.5Me') + 
  scale_y_continuous(limits = c(0.07, 0.3))+ theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

varimaxpc1 <- ggplot() + geom_line(aes(x = ages,  y = scoresScaled[,1])) + geom_point(aes(x = ages,  y = scoresScaled[,1], fill = GDGT0.4), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(-2, 1), position = "right")+
  xlab(NULL) + ylab('PC1')+theme_classic() +theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

varimaxpc2 <- ggplot() + geom_line(aes(x = ages,  y = scoresScaled[,2])) + geom_point(aes(x = ages,  y = scoresScaled[,2], fill = GDGT0.4), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  xlab(NULL) + ylab('PC2')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

varimaxpc3 <- ggplot() + geom_line(aes(x = ages,  y = scoresScaled[,3])) + geom_point(aes(x = ages,  y = scoresScaled[,3], fill = GDGT0.4), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(-5, 1), position = "right")+
  xlab(NULL) + ylab('PC3')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

CaldCren <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$GDGT0.4)], y = GDGTindices$GDGT0.4[!is.na(GDGTindices$GDGT0.4)])) + geom_point(aes(x = ages[!is.na(GDGTindices$GDGT0.4)],  y = GDGTindices$GDGT0.4[!is.na(GDGTindices$GDGT0.4)], fill = GDGT0.4[!is.na(GDGTindices$GDGT0.4)]), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(0, 700))+
  xlab(NULL) + ylab('Cald/Cren')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

BITindex <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$BIT)], y = GDGTindices$BIT[!is.na(GDGTindices$BIT)])) + geom_point(aes(x = ages[!is.na(GDGTindices$BIT)],  y = GDGTindices$BIT[!is.na(GDGTindices$BIT)], fill = GDGT0.4[!is.na(GDGTindices$BIT)]), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(0.98, 1), position = "right")+
  xlab(NULL) + ylab('BIT')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")


HP5 <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$HP5)], y = GDGTindices$HP5[!is.na(GDGTindices$HP5)])) + geom_point(aes(x = ages[!is.na(GDGTindices$HP5)],  y = GDGTindices$HP5[!is.na(GDGTindices$HP5)], fill = GDGT0.4[!is.na(GDGTindices$HP5)]), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(0.5, 0.82), position = "right")+
  xlab(NULL) + ylab('HP5')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")


CBT <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$CBT5Me)], y = GDGTindices$CBT5Me[!is.na(GDGTindices$CBT5Me)])) + geom_point(aes(x = ages[!is.na(GDGTindices$CBT5Me)],  y = GDGTindices$CBT5Me[!is.na(GDGTindices$CBT5Me)], fill = GDGT0.4[!is.na(GDGTindices$CBT5Me)]), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(2, 4.6))+
  xlab(NULL) + ylab('CBT5Me')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(), legend.position = "none")

IR <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$IR)], y = GDGTindices$IR[!is.na(GDGTindices$IR)])) + geom_point(aes(x = ages[!is.na(GDGTindices$IR)],  y = GDGTindices$IR[!is.na(GDGTindices$IR)], fill = GDGT0.4[!is.na(GDGTindices$IR)]), col = 'black',  pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(0.05,0.165))+
  xlab('Age (cal yr BP)') + ylab('IR')+theme_classic()+theme(legend.position = "none")

FC <- ggplot() + geom_line(aes(x = ages[!is.na(GDGTindices$fC)], y = GDGTindices$fC[!is.na(GDGTindices$fC)])) + geom_point(aes(x = ages[!is.na(GDGTindices$fC)],  y = GDGTindices$fC[!is.na(GDGTindices$fC)], fill = GDGT0.4[!is.na(GDGTindices$fC)]), col = 'black', pch = 21, size = 2) +
  scale_fill_gradient2(high="#ee964b", mid=  "#f4d06f", low = "#07beb8", midpoint=midpt1, na.value = "grey50") + scale_x_reverse(limits = c(12200, 7800)) +
  scale_y_continuous(limits = c(0,0.06), position = "right")+
  xlab(NULL) + ylab('fC')+theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x = element_blank(),legend.position = "none")

GDGTstackplot<-plot_grid(hexa, varimaxpc1,varimaxpc2,varimaxpc3,MBT,HP5,CaldCren,BITindex,CBT, FC,IR, 
                         align = 'v', ncol = 1,
                         labels = c('A', 'B','C','D','E','F','G','H','I','J','K'), label_x = 0.1, label_y = 1)#+

GDGTstackplot


################
#-----------Figure S6
#brGDGT ternary diagram
brGDGT.TPH <- data.frame(tetra = rowSums(brGDGTfa[,c("brGDGT.Ia","brGDGT.Ib","brGDGT.Ic")],na.rm = TRUE),
                         penta = rowSums(brGDGTfa[,c("brGDGT.IIa.5Me", "brGDGT.IIa.6Me", "brGDGT.IIb.5Me",
                                                     "brGDGT.IIb.6Me", "brGDGT.IIc.5Me", "brGDGT.IIc.6Me")], na.rm = TRUE),
                         hexa = rowSums(brGDGTfa[,c("brGDGT.IIIa.5Me","brGDGT.IIIa.6Me","brGDGT.IIIb.5Me",
                                                    "brGDGT.IIIb.6Me","brGDGT.IIIc.5Me","brGDGT.IIIc.6Me")], na.rm = TRUE))                                                     

brGDGT_tern<-ggtern(data=brGDGT.TPH,aes(x=tetra,y=penta,z=hexa))+
  geom_point(aes(color = GDGTindices$GDGT0.4)) +
  scale_colour_gradientn(colours = c("#FEF07D", "#B1A543", "#642162"))

brGDGT_tern

################
#-----------Figure S8
#Time series of brGDGT-inferred temperature

# --------------------- Calculate temp using various calibrations -----------------------------------------
Russell_MAAT <- -1.21+32.42*GDGTindices$MBT.5me
GDGTtemps<-data.frame(Russell_MAAT)

brGDGT.IIb.6Me <- brGDGTfa$brGDGT.IIb.6Me
brGDGT.IIb.6Me[is.na(brGDGT.IIb.6Me)] <- 0
brGDGT.IIb.5Me <- brGDGTfa$brGDGT.IIb.5Me
brGDGT.IIb.5Me[is.na(brGDGT.IIb.5Me)] <- 0

GDGTtemps$Russell_SFS <- 23.81 - 31.02*brGDGTfa$brGDGT.IIIa.5Me - 41.91*brGDGT.IIb.5Me -51.59*brGDGT.IIb.6Me - 24.70*brGDGTfa$brGDGT.IIa.5Me + 68.80*brGDGTfa$brGDGT.Ib

GDGTtemps$Zhao_Calib <- -1.82 + 56.06*GDGTindices$MBT.5me

brGDGT.IIIa.6Me <- brGDGTfa$brGDGT.IIIa.6Me
brGDGT.IIIa.6Me[is.na(brGDGT.IIIa.6Me)] <- 0
GDGTtemps$Raberg_Full <- -8.06 +37.52*brGDGTfa$brGDGT.Ia-266.83*brGDGTfa$brGDGT.Ib^2+133.42*brGDGTfa$brGDGT.Ib+100.85*brGDGTfa$brGDGT.IIa.6Me^2+
  58.15*brGDGT.IIIa.6Me^2+12.79*brGDGTfa$brGDGT.IIIa.5Me

GDGTtemps$Raberg_Meth <- 92.9 + 63.84*brGDGTfaMeth$fIb.meth^2-130.51*brGDGTfaMeth$fIb.meth-28.77*brGDGTfaMeth$fIIa.5Me.meth^2-72.28*brGDGTfaMeth$fIIb.5Me.meth^2-
  5.88*brGDGTfaMeth$fIIc.5Me.meth^2+20.89*brGDGTfaMeth$fIIIa.5Me.meth^2-40.54*brGDGTfaMeth$fIIIa.5Me.meth-80.47*brGDGTfaMeth$fIIIb.5Me.meth

GDGTtemps$ages <- ages

# --------------Plot temps----------
GDGTtemp_Plot <- ggplot() + geom_line(aes(x = GDGTtemps$ages, y = GDGTtemps$Russell_MAAT), col = 'purple') +
  geom_line(aes(x = GDGTtemps$ages, y = GDGTtemps$Russell_SFS), col = 'black') +
  geom_line(aes(x = GDGTtemps$ages, y = GDGTtemps$Zhao_Calib), col = 'light blue') +
  geom_line(aes(x = GDGTtemps$ages, y = GDGTtemps$Raberg_Full), col = 'red') +
  geom_line(aes(x = GDGTtemps$ages, y = GDGTtemps$Raberg_Meth), col = 'orange') +
  xlim(12200, 7500) +ylim(-2.1, 10) + ylab('Temperature') + theme_minimal()

GDGTtemp_Plot

# --------------Ribbon Plot Zhao-inferred temp----------

C.proxy = GDGTtemps$Zhao_Calib
idx = c(148, 149,150)#these are the three data points excluded bc likely derived from non-lacustrine source
C.proxy <- C.proxy[-idx]
ages <-ages[-idx]
GDGTanalyses <- C$paleoData[[1]]$measurementTable[[1]]$peakArea1050brGDGTIIIa5Me$values
C.ae = selectData(C,var.name = "ageensemble")$values
C.ae = C.ae[!is.na(GDGTanalyses),]#exclude ages for which we do not have GDGT measurements
C.ae = C.ae[-idx,]#exclude ages of samples that had a non-lacustrine source

time.range <- c(7000,12500)

x=ages[!is.na(C.proxy)]
y=C.proxy[!is.na(C.proxy)]

tsPlot <- plotTimeseriesEnsRibbons(X=C.ae,Y=C.proxy)+
  geom_line(mapping = aes(x = x,y = y),color = "red") +
  scale_x_reverse(name = "Age (yr BP)",limits = time.range[order(-time.range)] ) #+ 

tsPlot 

################
#-----------Figure S9
#Fatty acid chain length distributions and hydrogen isotope values

my_palette = c(brewer.pal(11, "PRGn")[c(11,10,9,5,4,3,2)])
my_palette

##### Distribution through time
C20peak <- selectData(C,var.name = "peakAreaC20FattyAcidMethylEster")$values
C22peak <- selectData(C,var.name = "peakAreaC22FattyAcidMethylEster")$values
C24peak <- selectData(C,var.name = "peakAreaC24FattyAcidMethylEster")$values
C26peak <- selectData(C,var.name = "peakAreaC26FattyAcidMethylEster")$values
C28peak <- selectData(C,var.name = "peakAreaC28FattyAcidMethylEster")$values
C30peak <- selectData(C,var.name = "peakAreaC30FattyAcidMethylEster")$values
C32peak <- selectData(C,var.name = "peakAreaC32FattyAcidMethylEster")$values
ages <- selectData(C,var.name = "agemedian")$values 

FA <- data.frame(C20peak,C22peak,C24peak,C26peak,C28peak,C30peak,C32peak)

FAfa <- FA/rowSums(FA, na.rm = FALSE)
FAfa$ages <-ages
FAfa <- FAfa[rowSums(FAfa[,c(1:7)], na.rm = TRUE)>0,]

waxDistPlotdf = FAfa
waxDistPlotdf <- waxDistPlotdf %>%                                   # Apply pivot_longer function
  pivot_longer(!ages, names_to = "compound", values_to = "concentration")

waxDistPlot <- ggplot(data=waxDistPlotdf, aes(x=ages, y=concentration, fill=compound)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=my_palette) +
  scale_x_reverse(c(12500,7500))+
  theme_minimal()+
  xlab('Age (cal yr BP)') + ylab('Fatty Acid Fractional Abundance')
waxDistPlot

#time series of ACL
ACL <- selectData(C,var.name = "ACL20to30")$values
ages <- selectData(C,var.name = "agemedian")$values 
plot(ages[!is.na(ACL)],ACL[!is.na(ACL)], type = "l")

#-------------stats about the values of fractional abundances-----------------

wax_faStats <- data.frame(matrix(nrow=8, ncol=0))

wax_faStats$value <- colnames(FAfa)
wax_faStats$median <- apply(FAfa,2,median,na.rm = TRUE)
wax_faStats$stdev <- apply(FAfa,2,sd,na.rm = TRUE)
wax_faStats$max <- apply(FAfa,2,max,na.rm = TRUE)
wax_faStats$min <- apply(FAfa,2,min,na.rm = TRUE)

write.csv(wax_faStats,"/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/CF8_Biomarkerplots_2022/CF8wax_faStats.csv")


#============modern waxes======================
QPTWaxes <- read.csv("/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/CF8_Biomarkerplots_2022/Hollister2022_ModernPlantWaxDist.csv")
QPTfa <- QPTWaxes[,c(2:8)]/rowSums(QPTWaxes[,c(2:8)],na.rm = FALSE)
QPTfa$plants <- QPTWaxes$X

QPTwaxDistPlotdf <- QPTfa %>%                                   # Apply pivot_longer function
  pivot_longer(!plants, names_to = "compound", values_to = "concentration")

QPTwaxDistPlot <- ggplot(data=QPTwaxDistPlotdf, aes(x=plants, y=concentration, fill=compound)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=my_palette) +
  #scale_x_reverse(c(12500,7500))+
  theme_minimal()+
  xlab('plant taxon') + ylab('Fatty Acid Fractional Abundance')
QPTwaxDistPlot

QPTACL <- QPTWaxes[,c(2:8)]
QPTACL$plants <- QPTWaxes$X
QPTACL[c(1)]<-QPTACL[c(1)]*20
QPTACL[c(2)]<-QPTACL[c(2)]*22
QPTACL[c(3)]<-QPTACL[c(3)]*24
QPTACL[c(4)]<-QPTACL[c(4)]*26
QPTACL[c(5)]<-QPTACL[c(5)]*28
QPTACL[c(6)]<-QPTACL[c(6)]*30
QPTACL[c(7)]<-QPTACL[c(7)]*32
#FAdf <- FAdf[rowSums(FAdf[,c(1:2)], na.rm = TRUE)>0,]
#ACL20to32 = rowSums(ACLdf[,c(1,2,3,4,5,6,7)],na.rm = TRUE)/
# rowSums(FAdf[,c(1,3,5,7,9,11,13)],na.rm = TRUE)
QPTACL$ACL = rowSums(QPTACL[,c(1,2,3,4,5,6,7)],na.rm = TRUE)/
  rowSums(QPTfa[,c(1,2,3,4,5,6,7)],na.rm = TRUE)

#-------------------------- leaf wax d2H data and figures ---------------------------------------
C20d2H <- selectData(C,var.name = "d2HC20")$values
C22d2H <- selectData(C,var.name = "d2HC22")$values
C24d2H <- selectData(C,var.name = "d2HC24")$values
C26d2H <- selectData(C,var.name = "d2HC26")$values
C28d2H <- selectData(C,var.name = "d2HC28")$values
C30d2H <- selectData(C,var.name = "d2HC30")$values
depth <- selectData(C,var.name = "midpointdepth")$values
ages <- selectData(C,var.name = "agemedian")$values

#to save data so you don't have to upload and choose from the LiPD file again if you restart R
d2Hbigdf <- data.frame(ages,depth,C20d2H,C22d2H, C24d2H,C26d2H,C28d2H,C30d2H)
write.csv(d2Hbigdf,"/Users/Elizabeth/Documents/R/GDGTplots/CF8_GDGT_EH/CF8_EH_FAMEd2H_inputdata.csv")
#FA isotope values through time
waxd2HPlotdf = d2Hbigdf[c(1,3,4,5,6,7,8)]
waxd2HPlotdf <- waxd2HPlotdf %>%                                   # Apply pivot_longer function
  pivot_longer(!ages, names_to = "compound", values_to = "concentration")
waxd2HPlotdf <- na.omit(waxd2HPlotdf)

waxd2HPlot <- ggplot(waxd2HPlotdf, aes(x=ages, y=concentration, group=compound)) +
  geom_line(aes(color=compound))+
  #geom_point(aes(color=compound))+
  scale_color_manual(values=my_palette) +
  scale_x_reverse(c(12500,7500))+
  xlab('Age (cal yr BP)') + ylab('Fatty Acid d2H value (per mil VSMOW)')

waxd2HPlot

################
#-----------Table 1, Figure 3, Figure S10
#Determine correlation between brGDGT-inferred temperature and leaf wax d2H
cf8.ens <- selectData(C,var.name = "ageEnsemble")
cf8.agemed <- selectData(C, var.name = "ageMedian")
cf8.d2hc28 <- selectData(C,var.name = "d2hc28")
cf8.d2hc22 <- selectData(C,var.name = "d2hc22")
cf8.temp <- selectData(C,var.name = "ZhaoIWT")
cf8.eps <- selectData(C,var.name = "epsilonc28-c22")


#############

#assign isotope variable
cf8.iso <- cf8.d2hc28 #change to cf8.d2Hc28 OR cf8.d2Hc22 OR cf8.eps

#####plot if you want to double check
plot(cf8.agemed$values[!is.na(cf8.iso$values)],cf8.iso$values[!is.na(cf8.iso$values)], type = "l",col ="blue",xlim = c(12500,7000))#,ylim = c(-250,-200))
par(new=T)
plot(cf8.agemed$values[!is.na(cf8.temp$values)],cf8.temp$values[!is.na(cf8.temp$values)], type = "l",col ="purple",xlim = c(12500,7000))#,ylim = c(15,-30))


#make a dataframe of age, temp, iso
cf8Corrdf = data.frame(cf8.agemed$values,cf8.iso$values,cf8.temp$values)
#then remove all rows with NA
cf8Corrdf <- na.omit(cf8Corrdf)  

#############
#Remove millennial scale variability in the two records

cf8.iso.spline = smooth.spline(x=cf8Corrdf$cf8.agemed.values, y=cf8Corrdf$cf8.iso.values, spar = 0.7, keep.data = TRUE)
#quick plot to check the smooth
plot(cf8.iso.spline$x,cf8.iso.spline$y, type = "l",col ="blue",ylim = c(-300,-200))
par(new=T)
plot(cf8.iso.spline$x,cf8.iso.spline$yin, type = "l",col ="black",ylim = c(-300,-200))

#subtract spline from raw data to yield centennial scale variability
cf8.iso.cent = cf8.iso.spline$yin-cf8.iso.spline$y
plot(cf8.iso.spline$x,cf8.iso.cent, type = "l",col ="blue")#,ylim = c(-250,-200))

####now spline for T
cf8.temp.spline = smooth.spline(x=cf8Corrdf$cf8.agemed.values, y=cf8Corrdf$cf8.temp.values, spar = 0.7, keep.data = TRUE)
#quick plot to check the smooth
plot(cf8.temp.spline$x,cf8.temp.spline$y, type = "l",col ="purple",ylim = c(3,9))
par(new=T)
plot(cf8.temp.spline$x,cf8.temp.spline$yin, type = "l",col ="black",ylim = c(3,9))

#subtract spline from raw data to yield centennial scale variability
cf8.temp.cent = cf8.temp.spline$yin-cf8.temp.spline$y
plot(cf8.temp.spline$x,cf8.temp.cent, type = "l",col ="purple")#,ylim = c(-250,-200))
par(new=T)
plot(cf8.iso.spline$x,cf8.iso.cent, type = "l",col ="blue",ylim = c(15,-30))

#############
#Correlation between iso and GDGT-inferred temp, full record, from same samples, unfiltered + filtered millennial and centennial

resPears.unfilt<-cor.test(cf8Corrdf$cf8.temp.values,cf8Corrdf$cf8.iso.values, method="pearson")
resPears.mill<-cor.test(cf8.temp.spline$y,cf8.iso.spline$y, method="pearson")
resPears.cent<-cor.test(cf8.temp.cent,cf8.iso.cent, method="pearson")
#resPears.unfilt
#resPears.mill
#resPears.cent

#############
#Correlation between iso and GDGT-inferred temp, from same samples, unfiltered + filtered millennial and centennial
brkPt <- c(20:105)
corrdf <- data.frame(matrix(ncol =13, nrow = length(brkPt)))
colnames(corrdf)<-c('brkPtAge','r.early','pvalue.early','r.late','pvalue.late',
                    'r.early.mill','pvalue.early.mill','r.late.mill','pvalue.late.mill',
                    'r.early.cent','pvalue.early.cent','r.late.cent','pvalue.late.cent')

for(i in brkPt) {
  cf8.temp.late <- cf8Corrdf$cf8.temp.values[1:i]
  cf8.temp.early <- cf8Corrdf$cf8.temp.values[(i-1):112]
  cf8.iso.late <- cf8Corrdf$cf8.iso.values[1:i]
  cf8.iso.early <- cf8Corrdf$cf8.iso.values[(i-1):112]
  resPears.unfilt.early<-cor.test(cf8.temp.early,cf8.iso.early, method="pearson")
  resPears.unfilt.late<-cor.test(cf8.temp.late,cf8.iso.late, method="pearson")
  corrdf$brkPtAge[i-19] <- cf8Corrdf$cf8.agemed.values[i]
  corrdf$r.early[i-19] <- resPears.unfilt.early$estimate
  corrdf$pvalue.early[i-19] <- resPears.unfilt.early$p.value
  corrdf$r.late[i-19] <- resPears.unfilt.late$estimate
  corrdf$pvalue.late[i-19] <- resPears.unfilt.late$p.value
  
  cf8.temp.late.mill <- cf8.temp.spline$y[1:i]
  cf8.temp.early.mill <- cf8.temp.spline$y[(i-1):112]
  cf8.iso.late.mill <- cf8.iso.spline$y[1:i]
  cf8.iso.early.mill <- cf8.iso.spline$y[(i-1):112]
  resPears.unfilt.early.mill<-cor.test(cf8.temp.early.mill,cf8.iso.early.mill, method="pearson")
  resPears.unfilt.late.mill<-cor.test(cf8.temp.late.mill,cf8.iso.late.mill, method="pearson")
  corrdf$r.early.mill[i-19] <- resPears.unfilt.early.mill$estimate
  corrdf$pvalue.early.mill[i-19] <- resPears.unfilt.early.mill$p.value
  corrdf$r.late.mill[i-19] <- resPears.unfilt.late.mill$estimate
  corrdf$pvalue.late.mill[i-19] <- resPears.unfilt.late.mill$p.value
  
  cf8.temp.late.cent <- cf8.temp.cent[1:i]
  cf8.temp.early.cent <- cf8.temp.cent[(i-1):112]
  cf8.iso.late.cent <- cf8.iso.cent[1:i]
  cf8.iso.early.cent <- cf8.iso.cent[(i-1):112]
  resPears.unfilt.early.cent<-cor.test(cf8.temp.early.cent,cf8.iso.early.cent, method="pearson")
  resPears.unfilt.late.cent<-cor.test(cf8.temp.late.cent,cf8.iso.late.cent, method="pearson")
  corrdf$r.early.cent[i-19] <- resPears.unfilt.early.cent$estimate
  corrdf$pvalue.early.cent[i-19] <- resPears.unfilt.early.cent$p.value
  corrdf$r.late.cent[i-19] <- resPears.unfilt.late.cent$estimate
  corrdf$pvalue.late.cent[i-19] <- resPears.unfilt.late.cent$p.value
  
}




###################
#find time window when all correlations are strongest, use this time range as the breakpoint, calculate some age stats
#cutoffs for C28
#brkPt <- subset(corrdf,(pvalue.early.mill<0.001 & pvalue.late.mill<0.001 & pvalue.early.cent<0.001))

#cutoffs for C22
brkPt <- corrdf[65:66,]

brkPtAgeEns <- c(cf8.ens$values[104,],cf8.ens$values[106,])
brkPtAgeUnc <- summary(brkPtAgeEns)
brkPtAgeUnc
brkPtAgesd <- sd(brkPtAgeEns)
brkPtAgesd

###################
#down-sample the record to make sure the correlation isn't due to different resolution
#make a dataframe of age, temp, iso binned at regular intervals, so not higher resolution in upper part of record
cf8Biniso = geoChronR::bin(cf8Corrdf$cf8.agemed.values,cf8Corrdf$cf8.iso.values,bin.vec = seq(7000,12500,by=50),bin.fun = mean)
cf8BinT = geoChronR::bin(cf8Corrdf$cf8.agemed.values,cf8Corrdf$cf8.temp.values,bin.vec = seq(7000,12500,by=50),bin.fun = mean)
cf8Bin <- data.frame(matrix(ncol =3, nrow = length(cf8Biniso$x)))
colnames(cf8Bin)<-c('Age','cf8.iso','cf8.temp')
cf8Bin$Age <- cf8Biniso$x
cf8Bin$cf8.iso <- cf8Biniso$y
cf8Bin$cf8.temp <- cf8BinT$y

#then remove all rows with NA
cf8Bin <- na.omit(cf8Bin)  

#then do correlation again!
resPears.unfilt.bin<-cor.test(cf8Bin$cf8.temp,cf8Bin$cf8.iso, method="pearson")

brkPtbin <- c(10:60)
corrdfbin <- data.frame(matrix(ncol =5, nrow = length(brkPtbin)))
colnames(corrdfbin)<-c('brkPtAge','r.early','pvalue.early','r.late','pvalue.late')

for(i in brkPtbin) {
  cf8.temp.late.bin <- cf8Bin$cf8.temp[1:i]
  cf8.temp.early.bin <- cf8Bin$cf8.temp[(i-1):71]
  cf8.iso.late.bin <- cf8Bin$cf8.iso[1:i]
  cf8.iso.early.bin <- cf8Bin$cf8.iso[(i-1):71]
  resPears.unfilt.early.bin<-cor.test(cf8.temp.early.bin,cf8.iso.early.bin, method="pearson")
  resPears.unfilt.late.bin<-cor.test(cf8.temp.late.bin,cf8.iso.late.bin, method="pearson")
  corrdfbin$brkPtAge[i-9] <- cf8Bin$Age[i]
  corrdfbin$r.early[i-9] <- resPears.unfilt.early.bin$estimate
  corrdfbin$pvalue.early[i-9] <- resPears.unfilt.early.bin$p.value
  corrdfbin$r.late[i-9] <- resPears.unfilt.late.bin$estimate
  corrdfbin$pvalue.late[i-9] <- resPears.unfilt.late.bin$p.value
  
}

###########build a table of the correlation coefficient values

corcoeffdf <- data.frame(matrix(ncol = 6, nrow=4))
colnames(corcoeffdf) <- c('Full record r','Full record p','Before 9.8 ka r','Before 9.8 ka p','After 9.8 ka r','After 9.8 ka p')
rownames(corcoeffdf) <- c('Raw record','Millennial','Centennial','Binned record')
corcoeffdf[1,1] <- round(resPears.unfilt$estimate,digits = 2)
corcoeffdf[1,2] <- ifelse(round(resPears.unfilt$p.value,digits = 3)==0,"<0.001",(round(resPears.unfilt$p.value,digits = 2)))
corcoeffdf[1,3] <- round(brkPt[1,2], digits = 2)
corcoeffdf[1,4] <- ifelse(round(brkPt[1,3],digits = 3)==0,"<0.001",(round(brkPt[1,3],digits = 2)))
corcoeffdf[1,5] <- round(brkPt[1,4],digits = 2)
corcoeffdf[1,6] <- ifelse(round(brkPt[1,5],digits = 3)==0,"<0.001",(round(brkPt[1,5],digits = 2)))
corcoeffdf[2,1] <- round(resPears.mill$estimate,digits = 2)
corcoeffdf[2,2] <- ifelse(round(resPears.mill$p.value,digits = 3)==0,"<0.001",(round(resPears.mill$p.value,digits = 2)))
corcoeffdf[2,3] <- round(brkPt[1,6],digits = 2)
corcoeffdf[2,4] <- ifelse(round(brkPt[1,7],digits = 3)==0,"<0.001",(round(brkPt[1,7],digits = 2)))
corcoeffdf[2,5] <- round(brkPt[1,8],digits = 2)
corcoeffdf[2,6] <- ifelse(round(brkPt[1,9],digits = 3)==0,"<0.001",(round(brkPt[1,9],digits = 2)))
corcoeffdf[3,1] <- round(resPears.cent$estimate,digits = 2)
corcoeffdf[3,2] <- ifelse(round(resPears.cent$p.value,digits = 3)==0,"<0.001",(round(resPears.cent$p.value,digits = 2)))
corcoeffdf[3,3] <- round(brkPt[1,10],digits = 2)
corcoeffdf[3,4] <- ifelse(round(brkPt[1,11],digits = 3)==0,"<0.001",(round(brkPt[1,11],digits = 3)))
corcoeffdf[3,5] <- round(brkPt[1,12],digits = 2)
corcoeffdf[3,6] <- ifelse(round(brkPt[1,13],digits = 3)==0,"<0.001",(round(brkPt[1,13],digits = 2)))
corcoeffdf[4,1] <- round(resPears.unfilt.bin$estimate,digits = 2)
corcoeffdf[4,2] <- ifelse(round(resPears.unfilt.bin$p.value,digits = 3)==0,"<0.001",(round(resPears.unfilt.bin$p.value,digits = 2)))
corcoeffdf[4,3] <- round(corrdfbin[35,2],digits = 2)
corcoeffdf[4,4] <- ifelse(round(corrdfbin[35,3],digits = 3)==0,"<0.001",(round(corrdfbin[35,3],digits = 2)))
corcoeffdf[4,5] <- round(corrdfbin[35,4],digits = 2)
corcoeffdf[4,6] <- ifelse(round(corrdfbin[35,5],digits = 3)==0,"<0.001",(round(corrdfbin[35,5],digits = 2)))

#change the file name when writing csv for different isotope time series
#write.csv(corcoeffdf,file = "/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/CF8_T_eps_corrTable.csv")
##########Make plots####

plot(corrdf$brkPtAge,corrdf$r.early, type = "l",col ="blue",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdf$brkPtAge,corrdf$r.early.mill, type = "l",lty = 2,col ="blue",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdf$brkPtAge,corrdf$r.early.cent, type = "l",lty = 3,col ="blue",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdf$brkPtAge,corrdf$r.late, type = "l",col ="purple",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdf$brkPtAge,corrdf$r.late.mill, type = "l",lty = 2,col ="purple",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdf$brkPtAge,corrdf$r.late.cent, type = "l",lty = 3,col ="purple",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdfbin$brkPtAge,corrdfbin$r.early, type = "l",lwd=2,col ="blue",xlim = c(12500,7000),ylim = c(-1,1))
par(new=T)
plot(corrdfbin$brkPtAge,corrdfbin$r.late, type = "l",lwd=2,col ="purple",xlim = c(12500,7000),ylim = c(-1,1))


plot(corrdf$brkPtAge,corrdf$pvalue.early, type = "l",col ="blue",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdf$brkPtAge,corrdf$pvalue.early.mill, type = "l",lty = 2,col ="blue",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdf$brkPtAge,corrdf$pvalue.early.cent, type = "l",lty = 3,col ="blue",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdf$brkPtAge,corrdf$pvalue.late, type = "l",col ="purple",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdf$brkPtAge,corrdf$pvalue.late.mill, type = "l",lty = 2,col ="purple",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdf$brkPtAge,corrdf$pvalue.late.cent, type = "l",lty = 3,col ="purple",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdfbin$brkPtAge,corrdfbin$pvalue.early, type = "l",lwd=2,col ="blue",xlim = c(12500,7000),ylim = c(0,0.01))
par(new=T)
plot(corrdfbin$brkPtAge,corrdfbin$pvalue.late, type = "l",lwd=2,col ="purple",xlim = c(12500,7000),ylim = c(0,0.01))




################
#-----------Figure S12
#Determine temperature at the two main summer precipitation moisture sources: Eastern and Western North America. There are no available high-resolution Holocene data for the North Pacific, so we do not consider that source in these calculations.

#First import the data, filter according to region and resolution, and make a composite. Later, we'll make a stack plot and maps of the time series used in the composites.
#D <- readLipd("https://lipdverse.org/Temp12k/1_0_2/Temp12k1_0_2.zip") #this is the online version of temp12k Lipdverse
D <- readLipd("/Users/Elizabeth/Documents/R/Temperature12k/ScientificDataAnalysis/lipdFilesWithEnsembles") #this is the  version of temp12k Lipdverse with age uncertainty

TS <- #extractTs(D) %>% #extract to lipd-ts
  as.lipdTsTibble(D) %>% # and the then to lipd-ts-tibble for filtering
  #filter(between(geo_longitude,-105,-50)) %>%  #only NE NA longitudes
  #filter(between(geo_latitude,35,55)) %>%  #between 35 and 55 N
  filter(between(geo_longitude,-175,-90)) %>%  #only W NA longitudes
  filter(between(geo_latitude,55,80)) %>%  #between 55 and 80 N
  #filter(between(geo_longitude,-180,-125)| between(geo_longitude,160,180)) %>%  #only Pacific longitudes
  #filter(between(geo_latitude,0,55)) %>%  #between 20 and 50 N
  #this pacific region has NO temp records at <200 yr median resolution! :(, I tried filtering without the resolution criterion...still only one coastal record
  filter(interpretation1_variable == "T") %>% #only variables sensitive temperature
  filter(paleoData_medianRes12k < 200) %>% #only time series at highres
  #filter(paleoData_medianRes12k < 50) %>% #only time series at highest res
  filter(interpretation1_seasonalityGeneral == "summer+") %>% #only summer proxies
  as.lipdTs() #back to TS for compositeR



#bin the TS
binvec <-  seq(5800, to = 12500,by = 50)#so the composite overlaps with the 6 ka time slice
binYears <- rowMeans(cbind(binvec[-1],binvec[-length(binvec)]))


#setup ensemble
nens <- 20#this is the number of ensembles, which is used to estimate the uncertainty


registerDoParallel(2)#this is to use multiple computer cores to run this in parallel

ensOut <- foreach(i = 1:nens) %dopar% {
  tc <- compositeEnsembles(TS,
                           binvec,
                           stanFun = standardizeMeanIteratively,#how to standardize data, "standardize over interval" command has the 'bundle of sticks effect' where interval chose arbitrarily affects the uncertainty structure, so Nick has standardizeoverRandomInterval, StandardizeMeanIteratively finds mean value for each record, trying to find any solution that has the most agreement for all records
                           binFun = sampleEnsembleThenBinTS,#simpleBinTs does not incorporate age model uncertainty; sampleEnsembleThenBinTS to incorporate age and proxy uncertainty, to get age ensembles...get dataset names
                           ageVar = "age",
                           alignInterpDirection = FALSE,#we need this if there are interpretation directions, ie with isotope data, only works if interpretations are in there and correct
                           spread = TRUE,#for low res datasets, this spreads the info from a given datapoint to cover the distance halfway to the next value, before we bin it. This deals with aliasing
                           duration = 3000,
                           searchRange = c(5800,12500),
                           normalizeVariance = FALSE)#this is whether we take a z-score before we composite or not (do this if using multiple different proxies)
  
  return(list(composite = tc$composite,count = tc$count))
}


#thisComposite <-  as.matrix(purrr::map_dfc(ensOut,extract,"composite"))
thisComposite <-  as.matrix(purrr::map_dfc(ensOut,"composite"))


compPlot <- plotTimeseriesEnsRibbons(X = binYears,Y = thisComposite)+
  scale_x_reverse(name = "Year (BP)",oob = scales::squish)+
  scale_y_continuous(name = "Summer Temperature (°C)",oob = scales::squish)+
  theme_bw()+
  ggtitle("Temperature sensitive")

compPlot

medianComposite <- rowMedians(thisComposite)
thisComposite <- as.data.frame(thisComposite)
thisComposite$ages <- binYears
thisComposite$median <- medianComposite

#To convert to absolute temps, I estimated the mean temperature of the warmest month in each source region (WNA and ENA) at 6 ka using data from Bartlein et al., 2010
#Fig 4 panel shows modern temps: WNA MTWA is ~10°C, ENA MTWA is ~25°C
#Fig 6 panel shows 6 ka temps: WNA MTWA is ~ -5 to +5°C compared to modern, ENA MTWA is ~ -5 to +5°C compared to modern, but mostly warming in the Great Lakes region, south of the middle Holocene Laurentide
#So, call WNA the same temp as modern, convert to an absolute value of 10°C for WNA MTWA and call ENA slightly warmer than modern, convert to an absolute value of 25.5°C for ENA MTWA
thisComposite$absolute <- thisComposite$median+10 #convert to abs temp for WNA
#thisComposite$absolute <- thisComposite$median+25.5 #convert to abs temp for ENA

write.csv(thisComposite,file = "/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/ManuscriptCalcs/ENASource200yrRes/ENAcomposite.csv")
write.csv(thisComposite,file = "/Users/Elizabeth/Documents/R/GeoChronR/07CF8L4_2021/ManuscriptCalcs/WNAsource200yrRes/WNAcomposite.csv")

#Make a stack plot
tidyData <- tidyTs(TS,age.var = "age")%>% 
  arrange(geo_latitude) #and sort by latitude

myPlot <- plotTimeseriesStack(tidyData, time.var = "age", color.var =  "archiveType", #create another plot, storing the output
                              lab.size = 2,
                              color.ramp = "black",#if we specify one color, all the records will get that
                              fill.alpha = 0, #make the fill transparent
                              lab.space = 3, #add it a bit more space between the axis labels and ticks
                              line.size = 0.4) + #make the lines a bit thinner 
  scale_x_reverse()

myPlot + annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 9100, xmax = 9600, ymin = 0, ymax = length(unique(tidyData$dataSetName))+1,alpha = 0.2) +
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 7800, xmax = 8300, ymin = 0, ymax = length(unique(tidyData$dataSetName))+1,alpha = 0.2) +
  annotate(geom = "rect", colour = NA, fill = "yellow", xmin = 10150, xmax = 10250, ymin = 0, ymax = length(unique(tidyData$dataSetName))+1,alpha = 0.2) +
  annotate(geom = "rect", colour = NA, fill = "green", xmin = 10450, xmax = 10700, ymin = 0, ymax = length(unique(tidyData$dataSetName))+1,alpha = 0.2)



#Then make a map
tm <- mapTs(tidyData,
            global = FALSE,
            color = "archiveType", #color by the archiveType
            size = 4
)

print(tm)

##############
