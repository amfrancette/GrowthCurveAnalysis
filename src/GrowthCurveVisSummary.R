library(tidyverse)
library(dplyr)
library(ggplot2)
library(mgcv)  
library(lubridate)  
library(inflection)  
library(mosaic)


setwd("/Users/amf198/Documents/ProjectNotes/Project_Paf1C_Depletion_4tU/Paf1CAIDPhenotyping/GrowthCurves/GC_Data/")
# Set colorscheme
{
  
  unique(GC_Data_Longer$Mutant)
  cols <- c("WT" =  "black", "Paf1-AID-osTIR" = "goldenrod2", "Paf1-AID osTIR" = "lightgoldenrod4",
            "paf1D" = "darkgoldenrod4", "Blank" =  "grey",
            "Rtf1-AID-osTIR" = "deeppink", "Rtf1-AID" = "coral", "rtf1D" = "red",
            "Ctr9-AID-osTIR" = "orange", "Ctr9-AID" = "sienna1", "ctr9D" = "orangered",
            "Cdc73-AID-osTIR" = "blue", "Cdc73-AID" = "lightblue", "cdc73D" = "navy")    
}

# Load Data and Data Keys
# Data keys have a field that connects to the well number of the plate, 
# and stores other metadata such as unique observation_index, media, mutant,
# genetic target, and so on. It takes a minute to fill out, but it goes faster with
# a template, and makes everything work pretty smoothly
# make note of any wells that were dropped from the original Data 

GC_Data_Key_001 <- read.csv("DataKeys/CtCdGCDataKey_001.csv", header = T, sep = ",")
GC_Data_Key_001$Observation_Index <- paste(GC_Data_Key_001$Plate_Num, "_", GC_Data_Key_001$Well, sep = '')
GC_Data_Key_001$Condition <-  paste(GC_Data_Key_001$Mutant, "_", GC_Data_Key_001$Media, sep = '')
GC_Data_Key_001$BioRep_Index <-  paste(GC_Data_Key_001$Condition, "_", GC_Data_Key_001$BioRep, sep = '')
GC_Data_001 <- read.csv("DataTables/CtCdGCData_001.csv", header = T, sep = ",")
GC_Data_001 <- pivot_longer(GC_Data_001, -Time, names_to = "Well", values_to = "OD600" )
GC_Data_001_Longer <- merge(GC_Data_001, GC_Data_Key_001, by.x = 'Well', by.y = 'Well')

GC_Data_Key_002 <- read.csv("DataKeys/PaRtGCDataKey_002.csv", header = T, sep = ",")
GC_Data_Key_002$Observation_Index <- paste(GC_Data_Key_002$Plate_Num, "_", GC_Data_Key_002$Well, sep = '')
GC_Data_Key_002$Condition <-  paste(GC_Data_Key_002$Mutant, "_", GC_Data_Key_002$Media, sep = '')
GC_Data_Key_002$BioRep_Index <-  paste(GC_Data_Key_002$Condition, "_", GC_Data_Key_002$BioRep, sep = '')
GC_Data_002 <- read.csv("DataTables/PaRtGCData_002.csv", header = T, sep = ",")
GC_Data_002 <- pivot_longer(GC_Data_002, -Time, names_to = "Well", values_to = "OD600" )
GC_Data_002_Longer <-merge(GC_Data_002, GC_Data_Key_002, by.x = 'Well', by.y = 'Well')
# Note Wells G1 and G6 were dropped from original data as they never reached saturation. Clearly a technical issue. 

GC_Data_Key_003 <- read.csv("DataKeys/PaRtGCDataKey_003.csv", header = T, sep = ",")
GC_Data_Key_003$Observation_Index <- paste(GC_Data_Key_003$Plate_Num, "_", GC_Data_Key_003$Well, sep = '')
GC_Data_Key_003$Condition <-  paste(GC_Data_Key_003$Mutant, "_", GC_Data_Key_003$Media, sep = '')
GC_Data_Key_003$BioRep_Index <-  paste(GC_Data_Key_003$Condition, "_", GC_Data_Key_003$BioRep, sep = '')
GC_Data_003 <- read.csv("DataTables/PaRtGCData_003.csv", header = T, sep = ",")
GC_Data_003 <- pivot_longer(GC_Data_003, -Time, names_to = "Well", values_to = "OD600" )
GC_Data_003_Longer <-merge(GC_Data_003, GC_Data_Key_003, by.x = 'Well', by.y = 'Well')

GC_Data_Longer <- rbind(GC_Data_001_Longer, GC_Data_002_Longer, GC_Data_003_Longer)
GC_Data_Longer

# sets plate ID to examine if you are into that kind of thing
Plate = (GC_Data_Longer$Plate_Num == "2")
Plate1 = (GC_Data_Longer$Plate_Num == "1")
Plate2 = (GC_Data_Longer$Plate_Num == "2")
Plate3 = (GC_Data_Longer$Plate_Num == "3")
IncludeWT <- (GC_Data_Longer$Mutant == as.character("WT"))

# Numerifies Time column and calculates mean by techincial and biological replicate
{
  GC_Data_Longer$Time <-(as.numeric(hms(GC_Data_Longer$Time)))/(60*60)
  # aggregate by time for everything but Observation Index to average multiple tech reps
  GC_Data_Long <- aggregate(OD600 ~ Time + Plate_Num + Condition + Target + Mutant + Media + BioRep_Index, data = GC_Data_Longer, mean)
  
  # aggregate by time for everything but BioRep Index to average multiple biological reps
  GC_Data_Reduced <- aggregate(OD600 ~ Time + Plate_Num + Condition + Target + Mutant + Media , data = GC_Data_Long, mean)
  GC_Data_Reduced$ConditionByPlate <- paste(GC_Data_Reduced$Condition, GC_Data_Reduced$Plate_Num)
}



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
    scale_colour_manual(values = cols) +  xlim(0,24) +  ylim(0,2.5)
  
  GC_Long_IP <-  merge(GC_Data_Long, inflectionPointMatrix[,-2], by.x = 'BioRep_Index',  by.y = 'BioRep_Index' )
  Data_subset <- GC_Long_IP[(GC_Long_IP$Time > (GC_Long_IP$inflectionPoint - 3.48) & GC_Long_IP$Time < (GC_Long_IP$inflectionPoint - 0.48)),]

  # # Uncomment the lines below to check to see how your cut curves look
  # full traces of untransformed data
  ggplot(GC_Data_Long, aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols) +  xlim(0,24) +  ylim(0,2.5)
  
  # cut traces from untransformed data
  ggplot(Data_subset, aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols)  +  xlim(0,24) +  ylim(0,2.5)
  
  # full traces of log2 transformed data
  ggplot(GC_Data_Long[GC_Data_Long$BioRep_Index == "rtf1D_YPD_Aux_4",], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols)
  
  # cut traces of log2 transformed data
  ggplot(Data_subset[Data_subset$Mutant == "paf1D",], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols)  +  xlim(0,24) 
  
  ggplot(Data_subset[Data_subset$Mutant == as.character("Ctr9-AID-osTIR"),], aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
    geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
    scale_colour_manual(values = cols) +  xlim(0,24) +  ylim(-3,1)
  
  slopes <- as.data.frame(Data_subset %>% 
                            group_by(BioRep_Index, Media) %>% 
                            do({
                              mod = lm(log2(OD600) ~ Time, data = .)
                              data.frame(Intercept = coef(mod)[1],
                                         Slope = coef(mod)[2])
                            }))
  slopes$DT <- (1/slopes$Slope)
  DTMatrix <- as.data.frame(cbind(BioRep = slopes[slopes$Media == "YPD",c(1)], YPD_DT = slopes[slopes$Media == "YPD",c(5)],
                                  Aux_DT = slopes[slopes$Media == "YPD_Aux",5]))
  DTMatrix$DifDT <- (as.numeric(DTMatrix$Aux_DT) - as.numeric(DTMatrix$YPD_DT))
  write.csv(DTMatrix, file = "DTMatrix_summary.csv")
}

# constructing matricies for DT comparisons by Tag Presence and Condition

DTMatrix

ggplot(testdata, aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
  geom_line(aes(linetype = `Media`), alpha=0.8, size = 1 ) + theme_bw(base_size = 15) +
  scale_colour_manual(values = cols)  +  xlim(0,24) 

# modeling troubleshooting
{
  GC_Curve_Fits <- as.data.frame(GC_Data_Long %>% 
                                   group_by(BioRep_Index) %>% 
                                   do({
                                     mod = lm(log2(OD600) ~ Time + I(Time^2) + I(Time^3) + I(Time^4) + I(Time^5) + I(Time^6), data = .)
                                     modFun <- makeFun(mod)
                                     data.frame(Intercept = coef(mod)[1],
                                                Coef1 = coef(mod)[2],
                                                Coef2 = coef(mod)[3],
                                                Coef3 = coef(mod)[4],
                                                Coef4 = coef(mod)[5],
                                                Coef5 = coef(mod)[6],
                                                Coef6 = coef(mod)[7])
                                   }))
  
  GC_Curve_Fit_Derivative <- as.data.frame(GC_Curve_Fits$BioRep_Index)
  GC_Curve_Fit_Derivative$Intercept <- (GC_Curve_Fits$Coef1)
  colnames(GC_Curve_Fit_Derivative) <- c("BioRep_Index", "Intercept")
  GC_Curve_Fit_Derivative$Coef1 <- (2 * GC_Curve_Fits$Coef2)
  GC_Curve_Fit_Derivative$Coef2 <- (3 * GC_Curve_Fits$Coef3)
  GC_Curve_Fit_Derivative$Coef3 <- (4 * GC_Curve_Fits$Coef4)
  GC_Curve_Fit_Derivative$Coef4 <- (5 * GC_Curve_Fits$Coef5)
  GC_Curve_Fit_Derivative$Coef5 <- (6 * GC_Curve_Fits$Coef6)
  
  # as.data.frame(GC_Curve_Fit_Derivative %>% 
  #                 group_by(BioRep_Index) %>% 
  #                 do({
  #                   curve
  #                   data.frame(Intercept = coef(mod)[1],
  #                              Coef1 = coef(mod)[2],
  #                              Coef2 = coef(mod)[3],
  #                              Coef3 = coef(mod)[4],
  #                              Coef4 = coef(mod)[5],
  #                              Coef5 = coef(mod)[6],
  #                              Coef6 = coef(mod)[7])
  #                 }))  

coeffs <- GC_Curve_Fits[GC_Curve_Fits$BioRep_Index == "WT_YPD_1",]
fun1 <- function(x) {(coeffs$Intercept + (coeffs$Coef1 * x) + (coeffs$Coef2 * x^2) + 
                        (coeffs$Coef3 * x^3) + (coeffs$Coef4 * x^4) + (coeffs$Coef5 * x^5) +
                        (coeffs$Coef6 * x^6))}

Dcoeffs <- GC_Curve_Fit_Derivative[GC_Curve_Fit_Derivative$BioRep_Index == "WT_YPD_1",]
fun4 <- function(x) {(Dcoeffs$Intercept + (Dcoeffs$Coef1 * x) + (Dcoeffs$Coef2 * x^2) + 
                (Dcoeffs$Coef3 * x^3) + (Dcoeffs$Coef4 * x^4) + (Dcoeffs$Coef5 * x^5))}

curve(2^(fun1(x)), xlim = c(.5,24), ylim = c(-1,2)) 
curve(fun1(x), xlim = c(.5,24)) 
curve(fun4(x), xlim = c(.5,24)) 

  list(GC_Curve_Fits %>% 
                  group_by(BioRep_Index) %>% 
                  do({
                    funs = function(x){(GC_Curve_Fits$Intercept + (GC_Curve_Fits$Coef1 * x) + (GC_Curve_Fits$Coef2 * x^2) + 
                                          (GC_Curve_Fits$Coef3 * x^3) + (GC_Curve_Fits$Coef4 * x^4) + (GC_Curve_Fits$Coef5 * x^5) +
                                          (GC_Curve_Fits$Coef6 * x^6))}
                    as.data.frame(curveFunction = funs)
                  }))
  

  GC_Curve_Fits[GC_Curve_Fits$BioRep_Index == "rtf1D_YPD_3",]
  curve(fun1(x), xlim = c(0.5,24))
  
  fun2 <- function(x) (makeFun(coeffs$Intercept + (coeffs$Coef1 * x) + (coeffs$Coef2 * x^2) + 
                                 (coeffs$Coef3 * x^3) + (coeffs$Coef4 * x^4) + (coeffs$Coef5 * x^5) +
                                 (coeffs$Coef6 * x^6) ~ x))
  
  eval(fun2(x), 1:24)
  ?eval()
  
  coeffs <- GC_Curve_Fits[GC_Curve_Fits$BioRep_Index == "rtf1D_YPD_3",]
  fun1 <- function(x) {(coeffs$Intercept + (coeffs$Coef1 * x) + (coeffs$Coef2 * x^2) + 
                          (coeffs$Coef3 * x^3) + (coeffs$Coef4 * x^4) + (coeffs$Coef5 * x^5) +
                          (coeffs$Coef6 * x^6))}
  
  fun3 <- expression((coeffs$Intercept) + (coeffs$Coef1 * x) + (coeffs$Coef2 * x^2) + 
                       (coeffs$Coef3 * x^3) + (coeffs$Coef4 * x^4) + (coeffs$Coef5 * x^5) +
                       (coeffs$Coef6 * x^6))
  
  eval(fun3, 1:24)
  ?eval
  class(fun3)
  D(fun3, 'x')
  ## Example of fitting 6 coeff polynomial model to GC data
  {# Select a testdata subset
    testdata <- GC_Data_Long[GC_Data_Long$BioRep_Index == "rtf1D_YPD_3",]
    #log2 transform data
    testdata$OD600 <- log2(testdata$OD600)
    # plot 
    plot(testdata$Time,testdata$OD600, col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3, ylim = c(-3,2), xlim = c(-3,25)) 
    curve(makeFun(testmodel)(x), add = TRUE)
    deriv(as.expression(makeFun(testmodel)), "Time", func = TRUE)
    
    testmodel <- lm(OD600 ~ Time + I(Time^2) + I(Time^3) + I(Time^4) + I(Time^5) + I(Time^6), data = testdata)
    plot(makeFun(testmodel), xlim = c(0.5,24), add = T)
    myPredict <- predict(testmodel)
    summary(testmodel)$adj.r.squared
    ix <- sort(testdata$Time, index.return=T)$ix
    lines(testdata$Time[ix], myPredict[ix], col=2, lwd=2 )  
    con <- names(coef(testmodel))
  }
  }



# GC Pipeline Example
{
  # Normal GC example for well A1
  ggplot(GC_Data_Long[GC_Data_Long$BioRep_Index == "001_A1" | GC_Data_Longer$Observation_Index == "001_D1" ,], aes(x=Time, y=OD600, group = `Well`)) +  
    geom_line(aes(color = `Mutant`), alpha=0.8, size = 2) + theme_bw() + scale_x_continuous(breaks=seq(0,24,1)) + 
    scale_y_continuous(breaks=seq(0,2.4,.2)) + scale_color_manual(values = cols) +
    geom_vline(xintercept = c(2,8.486111,11.48))
  
  # Log GC example for well A1
  ggplot(GC_Data_Longer[GC_Data_Longer$Observation_Index == "001_A1" | GC_Data_Longer$Observation_Index == "001_D1" ,], aes(x=Time, y=log2(OD600), group = `Well`)) +  
    geom_line(aes(color = `Mutant`), alpha=0.8, size = 2) + theme_bw() + scale_x_continuous(breaks=seq(0,24,1)) + scale_color_manual(values = cols) +
    geom_vline(xintercept = c(2,8.486111,11.48))
}

ggplot(GC_Data_Long, aes(x=Time, y=OD600, color = `Target`, group = `BioRep_Index`, shape = `Media`)) +  geom_line(aes(), alpha=0.8, size = 1) + theme_bw() 
ggplot(GC_Data_Long, aes(x=Time, y=log2(OD600), color = `Target`, group = `BioRep_Index`)) +  geom_line(aes(), alpha=0.8, size = 1) + theme_bw() 
ggplot(GC_Data_Longer, aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  geom_line(aes(), alpha=0.8, size = 1) + theme_bw() + scale_colour_manual(values = cols)



# Paf1 GCs
{
  # Checking Tech Reps
  {  
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Paf1-AID-osTIR") & GC_Data_Longer$Plate_Num == "3",],
           aes(x=Time, y=OD600, color = `BioRep_Index`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() 
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Paf1-AID") | 
                            GC_Data_Longer$Mutant == "WT" & GC_Data_Longer$Plate_Num == "3",],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("paf1D") | 
                            GC_Data_Longer$Mutant == "WT" & GC_Data_Longer$Plate_Num == "3",],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Examining Combined Tech Reps
  {  
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Paf1-AID-osTIR") | 
                          GC_Data_Long$Mutant == as.character("WT") & GC_Data_Long$Plate_Num == Plate,],
           aes(x=Time, y=log2(OD600), color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[ GC_Data_Long$Plate_Num == 3 & (GC_Data_Long$Mutant == as.character("Rtf1-AID-osTIR") | 
                                                          GC_Data_Long$Mutant == as.character("WT")) ,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[ GC_Data_Long$Plate_Num == 3 & (GC_Data_Long$Mutant == as.character("rtf1D") | 
                                                          GC_Data_Long$Mutant == as.character("WT")) ,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Group photo
  {
    ggplot(GC_Data_Reduced[GC_Data_Reduced$Target == as.character("Paf1") | 
                             GC_Data_Reduced$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `ConditionByPlate`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Target == as.character("Paf1") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
}

# Rtf1 GCs
{
  # Checking Tech Reps
  {  
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Rtf1-AID-osTIR") | 
                            IncludeWT & Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Rtf1-AID") | 
                            IncludeWT & Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw()  +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("rtf1D") | 
                            IncludeWT & Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Examining Combined Tech Reps
  {  
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Rtf1-AID-osTIR") | 
                          GC_Data_Long$Mutant == as.character("WT") & GC_Data_Long$Plate_Num == Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Rtf1-AID") | 
                          GC_Data_Long$Mutant == as.character("WT") & GC_Data_Long$Plate_Num == Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("rtf1D") | 
                          GC_Data_Long$Mutant == as.character("WT") & GC_Data_Long$Plate_Num == Plate,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Group photo
  {
    ggplot(GC_Data_Reduced[GC_Data_Reduced$Target == as.character("Rtf1") | 
                             GC_Data_Reduced$Mutant == as.character("WT"),],
           aes(x=Time, y=log2(OD600), color = `Mutant`, group = `ConditionByPlate`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    
    ggplot(GC_Data_Long[GC_Data_Long$Target == as.character("Rtf1") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
}

# Ctr9 GCs
{
  # Checking Tech Reps
  {  
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Ctr9-AID-osTIR") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Ctr9-AID") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("ctr9D") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Examining Combined Tech Reps
  {  
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Ctr9-AID-osTIR") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Ctr9-AID") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("ctr9D") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Group photo
  {
    ggplot(GC_Data_Reduced[GC_Data_Reduced$Target == as.character("Ctr9") | 
                             GC_Data_Reduced$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `ConditionByPlate`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
}

# Cdc73 GCs
{
  # Checking Tech Reps
  {  
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Cdc73-AID-osTIR") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("Cdc73-AID") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Longer[GC_Data_Longer$Mutant == as.character("cdc73D") | 
                            IncludeWT,],
           aes(x=Time, y=OD600, color = `Mutant`, group = `Observation_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Examining Combined Tech Reps
  {  
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Cdc73-AID-osTIR") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("Cdc73-AID") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
    
    ggplot(GC_Data_Long[GC_Data_Long$Mutant == as.character("cdc73D") | 
                          GC_Data_Long$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `BioRep_Index`)) +  
      geom_line(aes(linetype = `Media`), alpha=.5, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
  
  # Group photo
  {
    ggplot(GC_Data_Reduced[GC_Data_Reduced$Target == as.character("Cdc73") | 
                             GC_Data_Reduced$Mutant == as.character("WT"),],
           aes(x=Time, y=OD600, color = `Mutant`, group = `ConditionByPlate`)) +  
      geom_line(aes(linetype = `Media`), alpha=.9, size = 1) + theme_bw() +
      scale_colour_manual(values = cols)
  }
}
