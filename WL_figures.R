# Code for creating figures for Waterlogged paper - lessons learned.
########
rm(list = ls())
ls()
# Get data with function
source("/Users/user/SpecScripts/WLgetdata_function.R")
dir <- '/Users/user/Dropbox/PhD Work/Thesis chapters/Citizen Science Chapter/Figures'

WLdata <- getWLdata(project = "WL", 
                    spectro.direct = "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_spectrodata", 
                    fppar.directory = "/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_spectrodata/WL_fppar", 
                    pathlength = 33)
WLdata <- as.data.frame(WLdata)

## Get data for our own data - from DBP project
source("/Users/user/SpecScripts/DBP_datacompile_function.R")
DBPdata <- getDBPdata(save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata',
                      CM.pre.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_pre_CM_PARAFAC",
                      pre.directory <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMSRaleigh")
DBPdata <- DBPdata[!rownames(DBPdata) %in% "119", ] #take out null last row

# if the columns are factors, convert to numeric
indx <- sapply(DBPdata, is.factor)
DBPdata[indx] <- lapply(DBPdata[indx], function(x) as.numeric(as.character(x)))

indx <- sapply(WLdata, is.factor)
WLdata[indx] <- lapply(WLdata[indx], function(x) as.numeric(as.character(x)))

# PLot DOC
#WL DOC
#plot(WLdata$DOCcorr, WLdata$NO3.N)
#plot(DBPdata$NPOC_DOC_uncorrected, DBPdata$NO3)
# 99 -102 in DBP are quite high? Something going wrong in the correction phase
# Use uncorrected for now. Only affects DOC measurements from the spectro::lyzer
# Note also that there are Nas within certain of the FI measurements... why?

#required packages
library(tidyr)
library(dplyr)
library(plyr)
library(stringr)
library(stringi)
library(ggplot2)
library(pacman)
library(ggplot2)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc)

# convert DOC into the right format for density/histogram plots
varmerge <- function(WLdata_var, DBPdata_var, varname){
  temp <- data.frame(t(cbind(t(WLdata_var), t(DBPdata_var))))
  code <- data.frame(t(cbind(t(WLdata$sample), t(DBPdata$sample))))
  temp$code <- code
  temp$code <- sapply(temp$code, function(x) stri_sub(x, 1, -5)) # get code column such that only "WL or DBp" is in the column
  colnames(temp)[1] <- varname
  return(temp)
}

# Variables to examine
DOC <- varmerge(WLdata$DOCcorr, DBPdata$NPOC_DOC_uncorrected, varname = "DOC_mgL")
NO3 <- varmerge(WLdata$NO3.N, DBPdata$NO3, varname = "NO3mgL")
SUVA <- varmerge(WLdata$SUVA1, DBPdata$SUVA, varname = "SUVA")
e2e3 <- varmerge(WLdata$e2e3, DBPdata$e2e3.dec, varname = "e2e3")
e4e6 <- varmerge(WLdata$e4e6, DBPdata$e4e6.dec, varname = "e4e6")
SR <- varmerge(WLdata$SR, DBPdata$SR.dec, varname = "SR")
FI <- varmerge(WLdata$FI, DBPdata$FI, varname = "FI")
HIX <- varmerge(WLdata$HIX_ohno_area, DBPdata$HIX_ohno_area, varname = "HIX_ohno")
BIX <- varmerge(WLdata$FrI, DBPdata$FrI, varname = "BIX")
peakA <- varmerge(WLdata$peakA, DBPdata$peakA, varname = "peakA")
peakC <- varmerge(WLdata$peakC, DBPdata$peakC, varname = "peakC")
peakB <- varmerge(WLdata$peakB, DBPdata$peakB, varname = "peakB")
peakT <- varmerge(WLdata$peakT, DBPdata$peakT, varname = "peakT")
peakT.C <- varmerge(WLdata$peakt.peakC, DBPdata$peakt.peakC, varname = "peakTC")
redox <- varmerge(WLdata$redox, DBPdata$redox, varname = "redox")
perprotein <- varmerge(WLdata$perprotein, DBPdata$perprotein, varname = "perprotein")

# Do density plots
# http://stackoverflow.com/questions/6939136/how-to-overlay-density-plots-in-r
########### colour blind colour palettes
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Default theme
theme = theme_set(theme_bw() + 
                    theme(strip.text = element_text(size=14),
                          axis.text=element_text(size=14, color = "black"), 
                          axis.title = element_text(size=14), legend.text = element_text(size=14),
                          plot.title = element_text(size=14)))

# DOC
DOCplot <- ggplot(DOC, aes(DOC_mgL, fill = code)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[2]),
                    name="Data\nSource",
                    breaks=c("DBP", "WL"),
                    labels=c(" Scientist", " Citizen")) + 
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.8, .85)) +
  geom_vline(aes(xintercept=mean(WLdata$DOCcorr, na.rm=T)),   # Ignore NA values for mean
             color=cbPalette[2], linetype="dashed", size=1) +
  geom_vline(aes(xintercept=mean(DBPdata$NPOC_DOC_uncorrected, na.rm=T)),   # Ignore NA values for mean
             color=cbPalette[6], linetype="dashed", size=1) +
  ggtitle("DOC Concentration") +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="DOC Concentration (mg/L)")

# NO3
NO3plot <- ggplot(NO3, aes(NO3mgL, fill = code)) + geom_density(alpha = 0.2) +
  scale_fill_manual(values=c(cbPalette[6], cbPalette[2]),
                    name="Data\nSource",
                    breaks=c("DBP", "WL"),
                    labels=c(" Scientist", " Citizen")) +
  theme(legend.title=element_blank()) +
  theme(legend.position=c(.8, .85)) +
  geom_vline(aes(xintercept=mean(WLdata$NO3.N, na.rm=T)),   # Ignore NA values for mean
             color=cbPalette[2], linetype="dashed", size=1) +
  geom_vline(aes(xintercept=mean(DBPdata$NO3, na.rm=T)),   # Ignore NA values for mean
             color=cbPalette[6], linetype="dashed", size=1) + 
  ggtitle(expression('NO'[3]*'- Concentration')) +
  ylab(expression('')) +
  xlab(bquote('NO'[3]*'- Concentration (mg/L)')) 

pdf(file=paste0(dir,"/WLFigures_DOCNO3.pdf"), width = 11, height = 8.5)
grid.arrange(DOCplot,NO3plot,ncol=2)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

### Plot scatterplot - DOC and NO3 concentrations
# bind NO3 and DOC together
DN <- cbind(DOC, NO3, SUVA, perprotein, FI)

# Plot as scatter plot to add points
pdf(file=paste0(dir,"/WLscatter_NO3DOC.pdf"), width = 11, height = 8.5)

ggplot(DN, aes(x = NO3mgL, y = DOC_mgL, color = code)) + geom_point() +
  scale_color_manual(values=c(cbPalette[6], cbPalette[3]),
                     name="Data\nSource",
                     breaks=c("DBP", "WL"),
                     labels=c(" Scientist"," Citizen")) +
  theme(legend.position=c(.8, .85)) +
  xlim(0, 3) +
  scale_shape_manual(values=c(5,2)) +
  ylab("DOC (mg/L)") +
  xlab(expression('NO'[3]*'- Concentration (mg/L)')) +
  #geom_point(aes(size = FI)) + 
  #scale_size_continuous(name = FI, range = c(0.5,3))
  geom_vline(aes(xintercept=mean(DN$NO3mgL, na.rm=T)),    # Add mean NO3 for both datasets
             color=cbPalette[3], linetype="solid", size=1) + 
  geom_vline(aes(xintercept=as.numeric(quantile(DN$NO3mgL, 0.25, na.rm = TRUE))),    # Add 25% NO3 for both datasets
             color=cbPalette[3], linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=as.numeric(quantile(DN$NO3mgL, 0.75, na.rm = TRUE))),    # Add 75% NO3 for both datasets
             color=cbPalette[3], linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=mean(DN$DOC_mgL, na.rm=T)),   # Add mean DOC for both datasets
             color=cbPalette[1], linetype="solid", size=1) +
  geom_hline(aes(yintercept=as.numeric(quantile(DN$DOC_mgL, 0.25, na.rm = TRUE))),   # Add 25% DOC for both datasets
             color=cbPalette[1], linetype="dashed", size=1) +
  geom_hline(aes(yintercept=as.numeric(quantile(DN$DOC_mgL, 0.75, na.rm = TRUE))),   # Add 75% DOC for both datasets
             color=cbPalette[1], linetype="dashed", size=1) 

# need to deal with NAs in FI...
dev.off() 

#########
### Figure 3 - specific stories
# Took out size = FI - too much information.

pdf(file=paste0(dir,"/WLstories.pdf"), width = 11, height = 8.5)
# plot mean, 1st and 3rd quartile for DOC and NO3 
DOCNO3 <- ggplot(data = DN, aes(x = NO3mgL, y = DOC_mgL, color = "white")) + geom_point() +
  scale_color_manual(values=c("white", "white")) +
  theme(legend.position="none") +
  xlim(0, 3) +
  scale_shape_manual(values=c(5,2)) +
  ylab("DOC (mg/L)") +
  xlab(expression('NO'[3]*'- Concentration (mg/L)')) +
  #geom_point(aes(size = FI)) + 
  #scale_size_continuous(name = FI, range = c(0.5,3))
  geom_vline(aes(xintercept=mean(DN$NO3mgL, na.rm=T)),    # Add mean NO3 for both datasets
             color=cbPalette[3], linetype="solid", size=1) + 
  geom_vline(aes(xintercept=as.numeric(quantile(DN$NO3mgL, 0.25, na.rm = TRUE))),    # Add 25% NO3 for both datasets
             color=cbPalette[3], linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=as.numeric(quantile(DN$NO3mgL, 0.75, na.rm = TRUE))),    # Add 75% NO3 for both datasets
             color=cbPalette[3], linetype="dashed", size=1) + 
  geom_hline(aes(yintercept=mean(DN$DOC_mgL, na.rm=T)),   # Add mean DOC for both datasets
             color=cbPalette[1], linetype="solid", size=1) +
  geom_hline(aes(yintercept=as.numeric(quantile(DN$DOC_mgL, 0.25, na.rm = TRUE))),   # Add 25% DOC for both datasets
           color=cbPalette[1], linetype="dashed", size=1) +
  geom_hline(aes(yintercept=as.numeric(quantile(DN$DOC_mgL, 0.75, na.rm = TRUE))),   # Add 75% DOC for both datasets
           color=cbPalette[1], linetype="dashed", size=1)


### Evergreen - show points for Evergreen samples in the context of DOC and NO3 histrograms
# Evergreen samples = 46-56, 80-84, 121-124
sample <- c("WL0046", "WL0047", "WL0048", "WL0049", "WL0050", "WL0051", "WL0052",
            "WL0053", "WL0054", "WL0055", "WL0056", "WL0080", "WL0081", "WL0083",
            "WL0084", "WL0121", "WL0122", "WL0123")
EGsamples <- sample[sample %in% WLdata$sample]
#match samples that have data
ES <- subset(WLdata, WLdata$sample == EGsamples)
ES <- WLdata[which(WLdata$sample %in% EGsamples),]
# plot
St1 <-DOCNO3 +
  geom_point(data=ES, aes(x=ES$NO3.N, y=ES$DOCcorr, size = 1), colour=cbbPalette[2]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 1')) +
  theme(legend.position="none")

### Silver Creek Streamkeepers
# Prior to derailment: 157, 158, 165
# Down = 27, 164
up <- c("WL0157", "WL0158", "WL0165")
down <- c("WL0027", "WL0164")
all <- c("WL0157", "WL0158", "WL0165", "WL0027", "WL0164")
# match samples that have data
up <- WLdata[which(WLdata$sample %in% up),]
down <- WLdata[which(WLdata$sample %in% down),]
St2 <- DOCNO3 +
  geom_point(data=up, aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[2])  +
  geom_point(data=down, aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[7]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 2')) 

#### Bowen Island samples 218-221. 221 = control
BI <- c("WL0218", "WL0219", "WL0220", "WL0221")
BIs <- WLdata[which(WLdata$sample %in% BI),]
St3 <- DOCNO3 +
  geom_point(data=BIs, aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[7]) +
  geom_point(data=BIs[4,], aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[2]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 3')) 

### North Shore Streamkeepers
NSS <- c("WL0090", "WL0091", "WL0092", "WL0093", "WL0094", "WL0061", "WL0062", "WL0063", "WL0064", "WL0065",
        "WL0067", "WL0068", "WL0069", "WL0093", "WL0094")
NSS <- WLdata[which(WLdata$sample %in% NSS),]

St4 <- DOCNO3 +
  geom_point(data=NSS, aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[2]) +
  geom_point(data=NSS[6,], aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[7]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 4')) 

### North Shore  Zo Ann sampples. High DOC 186
ZA <- c("WL0036", "WL0038","WL0058",  "WL0074", "WL0075", "WL0076", "WL0077", "WL0078", "WL0079", 
        "WL0179", "WL0180", "WL0181", "WL0182", "WL0183", "WL0184", "WL0185", "WL0186", "WL0187", "WL0188" )
ZA <- WLdata[which(WLdata$sample %in% ZA),] 

St5 <- DOCNO3 +
  geom_point(data=ZA, aes(x = NO3.N, y = DOCcorr, size = 1), colour = cbbPalette[2]) +
  geom_point(data=ZA[8,], aes(x = NO3.N, y = DOCcorr, size = 1), colour = cbbPalette[7]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 5')) 

### Mt Polley Spill
MP <- c("WL0040")
MP <- WLdata[which(WLdata$sample %in% MP),] 
St6 <- DOCNO3 +
  geom_point(data=MP, aes(x=NO3.N, y=DOCcorr, size = 1), colour = cbbPalette[2]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 6')) 

### Lawson Creek Samples
#66, 67 - PWP
#102, 110, 130 - West Van Streamkeepers
LC <- c("WL0066", "WL0067", "WL0102", "WL0110", "WL0130")
LC <- WLdata[which(WLdata$sample %in% LC),] 
St7 <- DOCNO3 +
  geom_point(data=LC, aes(x=NO3.N, y=DOCcorr, size = 1), colour = cbbPalette[2]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 7')) 
# missing 66?

### Byrne Creek = Paul C and 
BC <- c("WL0155", "WL0193", "WL0194","WL0195","WL0197","WL0200")
BC <- WLdata[which(WLdata$sample %in% BC),] 
St8 <- DOCNO3 +
  geom_point(data=BC, aes(x=NO3.N, y=DOCcorr, size = 1), colour=cbbPalette[2]) +
  theme(legend.position="none") +
  ggtitle(expression('Story 8'))

# arrange in grid for saving as pdf
grid.arrange(St1, St2, St3, St4, St5, St6, ncol=2)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit((2/3), "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("D", x=unit(0.5, "npc"), y=unit((2/3), "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("E", x=unit(0, "npc")+ unit(2,"mm"), y=unit((1/3), "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("F", x=unit(0.5, "npc"), y=unit((1/3), "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("G", x=unit(0, "npc")+ unit(2,"mm"), y=unit((1/4), "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
#grid.text("H", x=unit(0.5, "npc"), y=unit(1/4, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
dev.off()

######## 
# Supplemental figures
# plot the two in histograms

#http://www.r-bloggers.com/how-to-make-a-histogram-with-ggplot2/

ggplot(data=DBPdata, aes(DBPdata$NPOC_DOC_uncorrected)) +
  geom_histogram(aes(y =..density..),
                 breaks=seq(0, 30, by = 2),
                 col="black",
                 fill="red",
                 alpha = .2) +
  geom_density(col="red") +
  labs(title="Our Data - [DOC]") +
  labs(x="DOC (mg/L)", y="Count") +
  geom_vline(aes(xintercept=mean(DBPdata$NPOC_DOC_uncorrected, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

ggplot(data=WLdata, aes(WLdata$DOCcorr)) +
  geom_histogram(aes(y =..density..),
                 breaks=seq(0, 30, by = 2),
                 col="black",
                 fill="blue",
                 alpha = .2) +
  geom_density(col="blue") +
  labs(title="Citizen Data - [DOC]") +
  labs(x="DOC (mg/L)", y="Count") +
  geom_vline(aes(xintercept=mean(WLdata$DOCcorr, na.rm=T)),   # Ignore NA values for mean
             color="blue", linetype="dashed", size=1)

#now make your lovely plot - CDF
# https://r-dir.com/blog/2014/03/cdfs-in-r.html
# http://stackoverflow.com/questions/6989015/ecdf-on-the-same-plot-using-ggplot2

ggplot(DOC, aes(x = DOC_mgL)) + 
  stat_ecdf(aes(group = code, colour = code)) 

# Statistical test for significant differences - table
# T-tests
DOCT <- t.test(DOC_mgL ~ code, data = DOC)
ggplot(DOC, aes(x=code, y=DOC_mgL, fill=code)) + geom_boxplot()

NO3T <- t.test(NO3mgL ~ code, data = NO3)
ggplot(NO3, aes(x=code, y=NO3mgL, fill=code)) + geom_boxplot()

SUVAt <- t.test(SUVA ~ code, data = SUVA)
ggplot(SUVA, aes(x=code, y=SUVA, fill=code)) + geom_boxplot()

e2e3T <- t.test(e2e3 ~ code, data = e2e3)
ggplot(e2e3, aes(x=code, y=e2e3, fill=code)) + geom_boxplot()

e4e6T <- t.test(e4e6 ~ code, data = e4e6)
ggplot(e4e6, aes(x=code, y=e4e6, fill=code)) + geom_boxplot()

FIT <- t.test(FI ~ code, data = FI)
ggplot(FI, aes(x=code, y=FI, fill=code)) + geom_boxplot()

HIXt <- t.test(HIX_ohno ~ code, data = HIX)
ggplot(HIX, aes(x=code, y=HIX_ohno, fill=code)) + geom_boxplot()

BIXt <- t.test(BIX ~ code, data = BIX)
ggplot(BIX, aes(x=code, y=BIX, fill=code)) + geom_boxplot()

t.test(peakA ~ code, data = peakA)
ggplot(peakA, aes(x=code, y=peakA, fill=code)) + geom_boxplot()

t.test(peakC ~ code, data = peakC)
ggplot(peakC, aes(x=code, y=peakC, fill=code)) + geom_boxplot()

t.test(peakT ~ code, data = peakT)
ggplot(peakT, aes(x=code, y=peakT, fill=code)) + geom_boxplot()

t.test(peakB ~ code, data = peakB)
ggplot(peakB, aes(x=code, y=peakB, fill=code)) + geom_boxplot()

t.test(redox ~ code, data = redox)
ggplot(redox, aes(x=code, y=redox, fill=code)) + geom_boxplot()

t.test(perprotein ~ code, data = perprotein)
ggplot(perprotein, aes(x=code, y=perprotein, fill=code)) + geom_boxplot()

test <- cbind(DOC, NO3, SUVA, e2e3, e4e6, FI, HIX, BIX, peakA, peakC, peakT, peakB, redox, perprotein)
test <- test[, !duplicated(colnames(test))]
means <- ddply(test, "code", colwise(mean), na.rm = TRUE)
stdevd <- ddply(test, "code", colwise(sd), na.rm = TRUE) # calculate st dev for all of the datasets 

sterrorDBP <- stdevd[1,-1]/sqrt(dim(subset(test, code =="DBP"))[1])
sterrorWL <- stdevd[2,-1]/sqrt(dim(subset(test, code =="WL"))[1])