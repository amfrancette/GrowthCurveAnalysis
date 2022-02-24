BiocManager::install("tidyverse")
BiocManager::install("ggplot2")
BiocManager::install("lubridate")
BiocManager::install("inflection")


library(tidyverse)
library(ggplot2)
library(lubridate) 
library(inflection) 
# library(ggpattern)

setwd("~/Documents/GitHub/GrowthCurveAnalysis/GC_Data/")

# Set colorscheme
{
  
  #unique(GC_Data_Longer$Mutant)
  cols <- c("WT" =  "grey60", 
            "Paf1-AID" = "#F5E733", "Paf1-AID-osTIR" = "#D1C52A", "paf1D" = "#9B921F", "Blank" =  "grey",
            "Rtf1-AID" = "#55CA61", "Rtf1-AID-osTIR" = "#409B4A",  "rtf1D" = "#2B6A32",
            "htz1D_Rtf1-AID" = "#48C9B0", "htz1D_Rtf1-AID-osTIR" = "#17A589",  "htz1D" = "#117864",
            "Ctr9-AID" = "#F7AE3F", "Ctr9-AID-osTIR" = "#CC9239",  "ctr9D" = "#A67831",
            "Cdc73-AID" = "#549EFC","Cdc73-AID-osTIR" = "#4586D9",  "cdc73D" = "#3D6394",
            "Leo1-AID" = "#EE46E9","Leo1-AID-osTIR" = "#CE3ECA",  "leo1D" = "#8E2A8B",
            "1" = "black", "2" = "red", "3" = "blue", "4" = "green" )  
  media <- c("YPD" = 1, "YPD_Aux" = 2, "YPD_Veh" = 3)
  sampleOrder <- c("Wt", "Paf1", "Ctr9", "Rtf1", 
                      "Cdc73", "Leo1", 
                   "WT", "Paf1-AID", "Paf1-AID-osTIR", "paf1D",
                   "Ctr9-AID", "Ctr9-AID-osTIR",  "ctr9D" ,
                   "Rtf1-AID", "Rtf1-AID-osTIR",  "rtf1D" ,
                   "Cdc73-AID","Cdc73-AID-osTIR",  "cdc73D" ,
                   "Leo1-AID","Leo1-AID-osTIR",  "leo1D" ,
                   "htz1D_Rtf1-AID", "htz1D_Rtf1-AID-osTIR" ,  "htz1D", "Blank")
}


# Load Data and Data Keys
# Data keys have a field that connects to the well number of the plate, 
# and stores other metadata such as unique observation_index, media, mutant,
# genetic target, and so on. It takes a minute to fill out, but it goes faster with
# a template, and makes everything work pretty smoothly
# make note of any wells that were dropped from the original Data 

{
# Read in data key csv files here. 
GC_Data_Key_1 <- read.csv("DataKeys/CtCdGCDataKey_001.csv", header = T, sep = ",")
GC_Data_Key_2 <- read.csv("DataKeys/PaRtGCDataKey_002.csv", header = T, sep = ",")
GC_Data_Key_3 <- read.csv("DataKeys/PaRtGCDataKey_003.csv", header = T, sep = ",")
GC_Data_Key_4 <- read.csv("DataKeys/CtLeGCDataKey_004.csv", header = T, sep = ",")
GC_Data_Key_5 <- read.csv("DataKeys/RtRthGCDataKey_005.csv", header = T, sep = ",")
GC_Data_Key_6 <- read.csv("DataKeys/RthLeGCDataKey_006.csv", header = T, sep = ",")
GC_Data_Key_7 <- read.csv("DataKeys/CtCdGCDataKey_007.csv", header = T, sep = ",")

# Read in data csv files here.
GC_Data_1 <- read.csv("DataTables/CtCdGCData_001.csv", header = T, sep = ",")
GC_Data_2 <- read.csv("DataTables/PaRtGCData_002.csv", header = T, sep = ",")
  # Note Wells G1 and G6 were dropped from original data of Plate 2 as they never reached saturation. Clearly a technical issue. 
GC_Data_3 <- read.csv("DataTables/PaRtGCData_003.csv", header = T, sep = ",")
GC_Data_4 <- read.csv("DataTables/CtLeGCData_004.csv", header = T, sep = ",")
GC_Data_5 <- read.csv("DataTables/RtRthGCData_005.csv", header = T, sep = ",")
GC_Data_6 <- read.csv("DataTables/RthLeGCData_006.csv", header = T, sep = ",")
GC_Data_7 <- read.csv("DataTables/CtCdGCData_007.csv", header = T, sep = ",")

# this for loop finds the objects matching the pattern "GC_Data_*" and trims them 
# leaving just the index. This index is used to perform a few operations on the data frames
# including making an observation index that appends plate and well ids, making a condition
# variable that combines mutant and media parameters, and making a bioRep index to be able to average 
# technical replicates. 

# this just cleans out the combined data frame of all plates if you wish to re-run the script in the same directory
rm(GC_Data_Longer)
for (i in str_remove(ls(pattern="GC_Data_Key_", all.names = TRUE), "GC_Data_Key_")) {
  print(paste("Plate ", i))
  # the get(eval()) piece allows me to interpret a string as a variable, I use this
  # to refer to variables with varying plate indices
  tempDataKey <- get(eval(paste0("GC_Data_Key_", i)))
  tempDataTable <- get(eval(paste0("GC_Data_", i)))
  tempDataKey$Observation_Index <- paste(tempDataKey$Plate_Num, "_", tempDataKey$Well, sep = '')
  tempDataKey$Condition <- paste(tempDataKey$Mutant, "_", tempDataKey$Media, sep = '')
  tempDataKey$BioRep_Index <-  paste(tempDataKey$Condition, "_", tempDataKey$BioRep, sep = '')
  # "longifies" data table to provide a "Well" column for merging with datakeys
  tempDataTable <- pivot_longer(tempDataTable,
                 -Time,
                 names_to = "Well",
                 values_to = "OD600")
  tempDataTable_Longer <-  merge(tempDataTable,
                                 tempDataKey,
          by.x = 'Well',
          by.y = 'Well')
  
  # this if else statement lets me make the variable GC_Data_Longer if it doesn't
  # already exist. If it does, then I append the curent temptable to the growing 
  # combined data frame
  if (exists("GC_Data_Longer")) {
    print(dim(GC_Data_Longer))
    GC_Data_Longer <- rbind(GC_Data_Longer, tempDataTable_Longer)
    print(dim(GC_Data_Longer))
  }  else {
    GC_Data_Longer <- tempDataTable_Longer
    print(dim(GC_Data_Longer))
  }
}

}

# removes htz1 samples from dataset
GC_Data_Longer <- GC_Data_Longer[!grepl("htz1",GC_Data_Longer$Target),]
GC_Data_Longer <- GC_Data_Longer[!grepl("YPD_Veh",GC_Data_Longer$Media),]

# Ties the genetic mutant to the biorep to keep track of which pairwise conditions came from the same colony/culture
GC_Data_Longer$Mutant_BioRep <- paste0(GC_Data_Longer$Mutant, "_" ,GC_Data_Longer$BioRep)

# Numerifies Time column and calculates mean by techincial and biological replicate
{
  GC_Data_Longer$Time <- (as.numeric(hms(GC_Data_Longer$Time)))/(60*60)
  # aggregate by time for everything but Observation Index to average multiple tech reps
  GC_Data_Long <- aggregate(OD600 ~ Time + Plate_Num + Condition + Target + Mutant + Media + BioRep_Index, data = GC_Data_Longer, mean)
  GC_Data_Long$sdByTechRep <- aggregate(OD600 ~ Time + Plate_Num + Condition + Target + Mutant + Media + BioRep_Index, data = GC_Data_Longer, sd)$OD600
  GC_Data_Long$sdByTechRep[is.na(GC_Data_Long$sdByTechRep)] <- 0
  # aggregate by time for everything but BioRep Index to average multiple biological reps
  GC_Data_Reduced <- aggregate(OD600 ~ Time + Condition + Target + Mutant + Media , data = GC_Data_Long, mean)
  GC_Data_Reduced$sdByBioRep <- aggregate(OD600 ~ Time + Condition + Target + Mutant + Media , data = GC_Data_Long, sd)$OD600
  GC_Data_Reduced$ConditionByPlate <- paste(GC_Data_Reduced$Condition, GC_Data_Reduced$Plate_Num)
  GC_Data_Reduced$sdByBioRep[is.na(GC_Data_Reduced$sdByBioRep)] <- 0
}

attach(GC_Data_Long)

# Keeps an untrimmed version for posterity
GC_Data_Longer_Untrimmed <- GC_Data_Longer
# removes Rtf1-AID BioRep2 bc of wierdness
GC_Data_Longer <- GC_Data_Longer[GC_Data_Longer$Mutant_BioRep != "Rtf1-AID_2",]
GC_Data_Longer <- GC_Data_Longer[GC_Data_Longer$Mutant_BioRep != "Leo1-AID-osTIR_3",]

{
  ggplot(GC_Data_Longer[GC_Data_Longer$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = .7 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols) +  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 3)
  
  ggplot(GC_Data_Long[GC_Data_Longer$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols) +  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 3)
  
  ggplot(GC_Data_Reduced[GC_Data_Reduced$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = ConditionByPlate)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) + geom_ribbon(aes(ymin=(OD600 - sdByBioRep), ymax=(OD600 + sdByBioRep), alpha = .05), linetype="blank", alpha=0.2)+
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols)+  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Target, ncol = 3)

}


# ifelse(GC_Data_Reduced$Target == "Rtf1", "Rtf1isTarget", "Rtf1isnotTarget")

# Extract Slopes!
{
  # calculates inflection points for combined BioReps
  {
    inflectionPointMatrix <- as.data.frame(GC_Data_Long[GC_Data_Long$Mutant != "Blank" & GC_Data_Long$Time > 1.9 & GC_Data_Long$Time < 20, ] %>%
                                             group_by(BioRep_Index, Media) %>%
                                             do({
                                               y = as.vector(.$OD600)
                                               x = as.vector(.$Time)
                                               cc = check_curve(x,y)
                                               ip = (bese(x,y,cc$index))
                                               data.frame(curveType = cc$ctype, inflectionPoint = ip$iplast)
                                             }))
    NoIPDatasets <- inflectionPointMatrix[rowSums(is.na(inflectionPointMatrix)) != 0,]
    inflectionPointMatrix <- inflectionPointMatrix[rowSums(is.na(inflectionPointMatrix)) == 0,]
  }
  
  # prints out if there are some datasets without an Inflection point
  NoIPDatasets
  
  # plots IP-less curves for visual inspection
  ggplot(GC_Data_Long[GC_Data_Long$BioRep_Index == NoIPDatasets$BioRep_Index,], aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols) +  xlim(0,25) +  ylim(0,2.5)
  
  # visualizing the pieces of curves that will be 
  # transforming Datasets to align with midpoint at T = 14
  GC_Data_Long_IP <- merge(GC_Data_Long, inflectionPointMatrix[,c(1,4)], by.x = 'BioRep_Index', by.y = 'BioRep_Index')
  GC_Data_Long_IP$Time = (GC_Data_Long_IP$Time - GC_Data_Long_IP$inflectionPoint)
  GC_Data_Long_IP$BioRep <- str_sub(GC_Data_Long_IP$BioRep_Index, -1)
  
  # plotting the curves shifted to align with midpoints and with a grey rectangle highlighting the line segment used for extracting slopes
  IPShifted_plot_log2 <- ggplot(GC_Data_Long_IP[GC_Data_Long_IP$Target != "Blank",], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = .5 ) + theme_bw(base_size = 8) +
    scale_colour_manual(values = cols) +  xlim(-14,14) +  ylim(-4,2.5) + facet_wrap(~Mutant, ncol = 3) + geom_vline(xintercept = 0 ) +
    annotate("rect", xmin = -4, xmax = -1, ymin = -4, ymax = 2.5,
             alpha = .2)
  ggsave(width = 10, height = 20, "../res/IPShifted_plot_log2.pdf", IPShifted_plot_log2)
  
  IPShifted_plot <- ggplot(GC_Data_Long_IP[GC_Data_Long_IP$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = .5 ) + theme_bw(base_size = 8) +
    scale_colour_manual(values = cols) +  xlim(-14,14) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 3) + geom_vline(xintercept = 0 ) +
    annotate("rect", xmin = -4, xmax = -1, ymin = 0, ymax = 2.5,
             alpha = .2)
  ggsave(width = 10, height = 20, "../res/IPShifted_plot.pdf", IPShifted_plot)
  

  GC_Long_IP <-  merge(GC_Data_Long, inflectionPointMatrix[,-2], by.x = 'BioRep_Index',  by.y = 'BioRep_Index' )
  Data_subset <- GC_Long_IP[(GC_Long_IP$Time > (GC_Long_IP$inflectionPoint - 3.48) & GC_Long_IP$Time < (GC_Long_IP$inflectionPoint - 0.48)),]
  
  # # Uncomment the lines below to check to see how your cut curves look
  # full traces of untransformed data
  allData_plot <- ggplot(GC_Data_Long[GC_Data_Long$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = .5 ) + theme_bw(base_size = 8) + geom_ribbon(aes(ymin=(OD600 - sdByTechRep), ymax=(OD600 + sdByTechRep), alpha = .2), linetype="blank", alpha=0.2)+
    scale_colour_manual(values = cols) +  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 3)
  ggsave(width = 10, height = 20, "../res/allData_plot.pdf", allData_plot)
  
  

        
        
  # full traces of log2 transformed data
  allData_plot_log2 <- ggplot(GC_Data_Long[GC_Data_Long$Target != "Blank",], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = .5 ) + theme_bw(base_size = 8) +
    scale_colour_manual(values = cols) +  xlim(0,25) +  ylim(-2.5,2) + facet_wrap(~Mutant, ncol = 3)
  ggsave(width = 10, height = 20, "../res/allData_plot_log2.pdf", allData_plot_log2)
  
  allData_Reduced_plot <-  ggplot(GC_Data_Reduced[GC_Data_Reduced$Target != "Blank",], aes(x=Time, y=OD600, color = `Mutant`, group = ConditionByPlate)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + 
    theme_bw(base_size = 15) + theme(legend.position  = "none") +
    geom_ribbon(aes(ymin=(OD600 - sdByBioRep),
                    ymax=(OD600 + sdByBioRep), 
                    alpha = .05), linetype="F1", 
                alpha=0.2, size = .3) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols)+  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 6)
  ggsave(width = 10, height = 5, "../res/allData_Reduced_plot.pdf", allData_Reduced_plot)
  
  # cut traces from untransformed data
  cutTrace_plot <- ggplot(Data_subset, aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols)  +  xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant, ncol = 3)
  ggsave(width = 10, height = 20, "../res/cutTrace_plot.pdf", cutTrace_plot)
  
  # cut traces from log2 transformed data
  cutTrace_plot_log2 <- ggplot(Data_subset, aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols)  +  xlim(0,25) +  ylim(-2.5,2.5) + facet_wrap(~Mutant, ncol = 3)
  ggsave(width = 10, height = 20, "../res/cutTrace_plot_log2.pdf", cutTrace_plot_log2)

  # making plate-specific plots to view tech rep data
  for (i in unique(GC_Data_Longer$Plate_Num)) {
    Plot_path <-  paste("../res/DataByPlate/PlatePlot_", i, ".pdf", sep = "")
    print(Plot_path)
    tempPlot <-  ggplot(GC_Data_Longer[GC_Data_Longer$Plate_Num == as.character(i) ,], aes(x=Time, y=OD600, color = as.factor(`BioRep`), group = `Observation_Index`)) +
      geom_line(aes(linetype = `Media`), alpha=0.8, size = 1) + theme_bw(base_size = 16) +
      ggtitle(paste0("Data from Plate ", i)) +
      xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant)
    ggsave(width = 10, height = 8, Plot_path, tempPlot)
  }
  
  # making plots by target
  for (i in unique(GC_Data_Long$Target)) {
    Plot_path <-  paste("../res/DataByTargetGene/TargetPlot_", i, ".pdf", sep = "")
    print(Plot_path)
    tempPlot <- ggplot(GC_Data_Long[GC_Data_Long$Target == i,], aes(x=Time, y=OD600, color = Mutant, group = `BioRep_Index`)) +
      geom_line(aes(linetype = `Media`), alpha=0.8, size = 1) + theme_bw(base_size = 16) +
      ggtitle(paste0(i, "Data")) + scale_color_manual(values = cols) +
      xlim(0,25) +  ylim(0,2.5) + facet_wrap(~Mutant)
    ggsave(width =20, height = 5, Plot_path, tempPlot)
  }

  slopes <- as.data.frame(Data_subset %>% 
                            group_by(BioRep_Index, Mutant, Media, Target) %>% 
                            do({
                              mod = lm(log2(OD600) ~ Time, data = .)
                              data.frame(Intercept = coef(mod)[1],
                                         Slope = coef(mod)[2])
                            }))
  slopes$DT <- (1/slopes$Slope)
  
  DTMatrix <- as.data.frame(cbind(Mutant = slopes$Mutant, 
                                  BioRep_Index = str_remove( str_remove(str_remove(slopes$BioRep_Index, "_YPD"), "_Aux"), 
                                                             "_Veh"), 
                                  Media = slopes$Media,
                                  Target = slopes$Target,
                                  DT = slopes$DT))
}


DTMatrix$BioRep <- str_sub(DTMatrix$BioRep_Index, -1)
DTMatrix$DT <- as.numeric(DTMatrix$DT)

DTMatrix <- DTMatrix[order(DTMatrix$BioRep_Index),]
DTMatrix$Target <- factor(DTMatrix$Target,
                          levels = sampleOrder)
DTMatrix$Mutant <- factor(DTMatrix$Mutant,
                          levels = sampleOrder)
attach(DTMatrix)

DTMatrix
DTMatrix_avg <- DTMatrix %>% 
  group_by(Mutant, Media, Target) %>% 
  summarise(meanDT = mean(DT), sdDT = sd(DT) )
DT_plot <- ggplot(DTMatrix_avg[DTMatrix_avg$Mutant != "htz1D_Rtf1-AID-osTIR",], aes(x=Mutant, y=meanDT, group=Media, fill = Mutant, color = Media)) + 
  geom_col( width=0.7, position=position_dodge(), alpha = 1) +
  geom_errorbar(aes(ymin=meanDT-sdDT, ymax=meanDT+sdDT), width=.2,
                position=position_dodge(.7)) +
  scale_fill_manual(values = cols) + 
  geom_point(data = DTMatrix[DTMatrix$Mutant != "htz1D_Rtf1-AID-osTIR",], aes(x=Mutant, y=DT, group=Media, fill = Mutant, color = Media, shape = as.factor(Media)),
             position=position_dodge(.7) ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("grey10","grey10","grey10")) + theme(legend.position = "none")
DT_plot
ggsave(width = 10, height = 6, "../res/DT_plot.pdf", DT_plot)
# looking at data in a bio-rep paired manner
attach(DTMatrix)



#### BELOW THIS POINT IS FOR TESTING #####


diff <- DTMatrix[(Media == "YPD"),]$DT - DTMatrix[Media == "YPD_Aux",]$DT
DTDiffMatrix <- cbind(DTDiffMatrix, diff, YPD_DT = DTMatrix[Media == "YPD",]$DT)

# DTDiffMatrix$percentdiff <- DTDiffMatrix$diff/DTDiffMatrix
DTDiffMatrix_avg <- DTDiffMatrix %>% 
  group_by(Mutant) %>% 
  summarise(
    meanDTDiff = mean(diff),
    sdDT = sd(diff)
  )

DT_diff_plot <- ggplot(DTDiffMatrix_avg, aes(x=Mutant, y=meanDTDiff, group=Mutant, fill = Mutant, color = Mutant)) + 
  geom_col( width=.9, position="dodge", alpha = 1, color = "black") +
  geom_errorbar(aes(ymin=meanDTDiff-sdDT, ymax=meanDTDiff+sdDT), width=.2,
                position=position_dodge(.7), color = "black") +
  scale_fill_manual(values = cols) + theme_bw(base_size = 10) +
  geom_point(data = DTDiffMatrix, aes(x=Mutant, y=diff),
             position=position_dodge(.7), color = "Black") 
DT_diff_plot
ggsave(width = 20, height = 6, "../res/DT_diff_plot.pdf", DT_diff_plot)


ggplot(DTDiffMatrix, aes(x=Mutant, y=diff, group=Mutant, color = Mutant)) + 
  geom_point(position=position_dodge(.7) ) +
  scale_color_manual(values = cols)


ggplot(GC_Data_Long_IP[GC_Data_Long_IP$Mutant == "Rtf1-AID-osTIR",], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
  geom_line(aes(linetype = `Media`), alpha=0.8, size = .5 ) + theme_bw(base_size = 8) +
  scale_colour_manual(values = cols) +  xlim(-14,14) +  ylim(-3,2.5) + facet_wrap(~BioRep, ncol = 3) + geom_vline(xintercept = 0 ) +
  annotate("rect", xmin = -4, xmax = -1, ymin = -3, ymax = 2.5,
           alpha = .2)
GC_Data_Longer_Untrimmed <- merge(GC_Data_Longer_Untrimmed, inflectionPointMatrix[,c(1,4)], by.x = 'BioRep_Index', by.y = 'BioRep_Index')
GC_Data_Longer_Untrimmed$IP_Time <-  (GC_Data_Longer_Untrimmed$Time - GC_Data_Longer_Untrimmed$inflectionPoint)

ggplot(GC_Data_Longer_Untrimmed[GC_Data_Longer_Untrimmed$Mutant == "Rtf1-AID",], aes(x=Time, y=log2(OD600), color = as.factor(`TechRep`), group = `Observation_Index`)) +  
  geom_line(aes(linetype = `Media`), alpha=0.5, size = 1 ) + theme_bw(base_size = 8) +
  scale_colour_manual(values = cols) +
  xlim(-6,5) +  ylim(-3,2.5) + facet_wrap(~Mutant_BioRep, ncol = 3) + geom_vline(xintercept = 0 ) +
  annotate("rect", xmin = -4, xmax = -1, ymin = -3, ymax = 2.5,
           alpha = .2)
ggplot(GC_Data_Longer_Untrimmed[GC_Data_Longer_Untrimmed$Mutant == "Rtf1-AID",], aes(x=Time, y=log2(OD600), color = as.factor(`TechRep`), group = `Observation_Index`)) +  
  geom_line(aes(linetype = `Media`), alpha=0.5, size = 1 ) + theme_bw(base_size = 8) +
  scale_colour_manual(values = cols) +
  xlim(-6,5) +  ylim(-3,2.5) + facet_wrap(~Mutant_BioRep, ncol = 3) + geom_vline(xintercept = 0 ) +
  annotate("rect", xmin = -4, xmax = -1, ymin = -3, ymax = 2.5,
           alpha = .2)

ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == c("ctr9D", "htz1D_Rtf1-AID-osTIR") ,], aes(x=Time, y=OD600, color = Mutant, group = `Observation_Index`)) +
  geom_line(aes(linetype = `Media`), alpha=0.8, size = .5) + theme_bw(base_size = 8) +
  ggtitle(paste0("Data from Plate ", 4)) +
  xlim(0,25) +  ylim(0,2.5) +  scale_colour_manual(values = cols)


DTMatrix_avg2 <- DTMatrix_avg[DTMatrix_avg$Media != "YPD_Veh",]
DTMatrix2 <- DTMatrix[DTMatrix$Media != "YPD_Veh",]

?linetype
ltype <- c("YPD" = 1, "YPD_Aux" = 14)


DT_plot <- ggplot(DTMatrix_avg2, aes(x=Mutant, y=meanDT, group=Media, fill = Mutant, color = Media, linetype = Media)) + 
  geom_col( width=0.7, position=position_dodge(), alpha = 1) +
  geom_errorbar(aes(ymin=meanDT-sdDT, ymax=meanDT+sdDT), width=.2,
                position=position_dodge(.7)) +
  scale_fill_manual(values = cols) + 
  geom_point(data = DTMatrix2, aes(x=Mutant, y=DT, group=Media, fill = Mutant, color = Media, shape = as.factor(Media)),
             position=position_dodge(.7) ) +
  scale_linetype_manual(values = ltype) +
  scale_color_manual(values = c("grey10","grey10","grey10")) + scale_y_continuous(breaks=seq(0,6,0.5)) + 
  geom_hline(yintercept = as.numeric(DTMatrix_avg[DTMatrix_avg$Mutant == "WT" & DTMatrix_avg$Media == "YPD", 4])) 
DT_plot

ggsave(width = 20, height = 6, "../res/DT_plotNoVeh.pdf", DT_plot)
?geom_col
?scale_linetype_manual()

