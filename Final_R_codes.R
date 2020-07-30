library(grid) # for rasterGrob
library(dplyr) # for join
library(tidyr) # for complete()
library(ggplot2) # for all plot
library(cowplot) # for ggdraw
library(ggpubr) # for compare_means
library(stringr) # for string manipulation
library(rcompanion) # for bootstrap CI for median - groupwiseMedian

# Time = vector of times
# relative2Time = Reference time point to treat as 0 hours (in this case, start of dusk), format H:MM
RelativeTime <- function(Time, relative2Time){
  Splittime <- strsplit(as.character(Time),":")
  RTime <- strsplit(relative2Time, ":")
  RelTimeMin <- 60*as.numeric(RTime[[1]][1]) + as.numeric(RTime[[1]][2]) #Relative time in minutes from 00:00, e.g. 16:30 will be 16*60+30 = 990 minutes from 00:00
  RelTimeMin2 <- rep(0,length(Splittime))
  RelativeT <- rep(0,length(Splittime))
  for (i in 1:length(Splittime)){
    RelTimeMin2[i] <- 60*as.numeric(Splittime[[i]][1]) +  as.numeric(Splittime[[i]][2])
    RelativeT[i] <- RelTimeMin2[i] - RelTimeMin
  }
  return(RelativeT)
}


# Must remove empty cells first
TransformCallus1 <- function(dataframe, Callus="Callus"){
  NewCallus <- ifelse(as.character(dataframe[,Callus]) %in% c("<<25","10"), 5,
                      ifelse(as.character(dataframe[,Callus]) %in% c("<25", "25","25+"), 20,
                             ifelse(as.character(dataframe[,Callus]) %in% c("25++","50","50+"), 50,
                                    ifelse(as.character(dataframe[,Callus]) %in% c("75","775"), 70,
                                           ifelse(as.character(dataframe[,Callus]) %in% c("75+","75++","95"), 90, 
                                                  ifelse(as.character(dataframe[,Callus]) %in% c("100","100+"),100,0 ))))))
  return(NewCallus)
}

# Must remove empty cells first
TransformCallus2 <- function(dataframe, Callus="Callus"){
  NewCallus <- ifelse(as.character(dataframe[,Callus]) %in% c("<<25","10","<25", "25","25+","25++","50","50+","75","775","75+","75++"), "Some yellow",
                      ifelse(as.character(dataframe[,Callus]) %in% c("95","100","100+"),"Yellow","Brown"))
  return(NewCallus)
}

# Based on the original ordinal characterisation
TransformCallus3 <- function(dataframe, Callus="Callus"){
  NewCallus <- ifelse(as.character(dataframe[,Callus]) %in% c("<<25","10","<25"), 'bbbY',
                      ifelse(as.character(dataframe[,Callus]) %in% c("25","25+"), 'bbYY',
                             ifelse(as.character(dataframe[,Callus]) %in% c("25++","50","50+"), 'bbYY',
                                    ifelse(as.character(dataframe[,Callus]) %in% c("75","775"), 'bYYY',
                                           ifelse(as.character(dataframe[,Callus]) %in% c("75+","75++","95"), 'bYYY', 
                                                  ifelse(as.character(dataframe[,Callus]) %in% c("100","100+"),'YYYY',
                                                         ifelse(is.null(dataframe$Line) == F ,ifelse(as.character(dataframe$Line) %in% c("S06"),"YYYY","bbbb"),"bbbb" )))))))
  return(NewCallus)
}


# Must remove empty cells first
TransformMating <- function(dataframe, Time="RelTime"){
  NewMatingTime <- ifelse(dataframe[,Time] >= 0, "Dusk-try",
                          ifelse(dataframe[,Time] >= -3, "Intermediate", "Day-neo"))
  return(NewMatingTime)
}

Expand_df2 <- function(dataframe, Name){
  A <- dataframe[0,]
  for (i in 1:length(dataframe[,Name])){
    A <- rbind(A,dataframe[rep(i,dataframe[,Name][i]),])
  }
  return(A)
}

SE <- function(x){
  return(sd(x)/sqrt(length(x)))
}

CI <- function(x){
  return(qnorm(0.975)*SE(x))
}


# Expand a dataset of frequencies (with column named 'Mating')
# data frame contains columns 'Line', 'Time' and 'Mating'
Expand_df <- function(dataframe, Name="Mating"){
  for (i in 1:length(dataframe[,Name])){
    if (dataframe[,Name][i] > 1){
      dataframe <- rbind(dataframe, data.frame(Line=rep(dataframe$Line[i], dataframe$Mating[i]-1),
                                               Time=rep(dataframe$Time[i], dataframe$Mating[i]-1),
                                               Mating=rep(dataframe$Mating[i], dataframe$Mating[i]-1)))  
    }
  }
  return(dataframe)
}


### splitvar is a string denoting the variable to which proportion are to be calculated within this variable, for facet_grid in ggplot
### xvar is the string denoting the variable/group that will be plotted on the x-axis
### df is the data frame containing the raw data to be counted and proportion evaluated
### Grouping is a vector of string names of the variables to group by... this must include the variable in splitvar and xvar
### To compute the count and proportion across groupings (nested within one splitvar)... e.g. if groupby c('A','B','C')... all combinations of 'A','B','C' will be generated
### If say 'B' is the splitvar, then proportion is calculated wihin each level of 'B', the proportions for each of 'B' should sum to 1
nest_prop <- function(df,Grouping,xvar,splitvar=""){
  library(dplyr)
  xvar <<- xvar ### for some reason, ddply needs to access xvar from outside the function, thus we need to force xvar to be a global variable
  DF1 <- ddply(df,Grouping,summarise, N = length(eval(parse(text=xvar))))
  #DF1 <- ddply(df,Grouping,summarise, N = eval(parse(text=paste('length(eval(parse(text=',xvar,')))'))) )
  if (splitvar != ""){
    for (i in 1:length(levels(df[,splitvar]))){
      if (i == 1){
        A <- DF1[DF1[,splitvar]==levels(df[,splitvar])[i],]
        A$Proportion <- A$N/sum(A$N)
      }
      else{
        B <- DF1[DF1[,splitvar]==levels(df[,splitvar])[i],]
        B$Proportion <- B$N/sum(B$N)
        A <- rbind(A,B)
      }
    }
  }
  else {
    A <- DF1
  }
  #print(A)
  return(as.data.frame(eval(parse(text=paste('complete(A,', paste(Grouping,collapse=','),')',sep=''))) )) # complete(A,...)
}

### Allow for replicates
### Works only if all groups have the same splitvar
nest_prop2 <- function(df,Grouping,xvar,splitvar,Rep){
  library(plyr)
  xvar <<- xvar
  DF1 <- ddply(df,c(Grouping,Rep),summarise, N = length(eval(parse(text=xvar))))
  if (splitvar != ""){
    for (i in 1:length(levels(df[,splitvar]))){
      if (i == 1){
        A <- DF1[DF1[,splitvar]==levels(df[,splitvar])[i],]
        A[,splitvar] <- factor(A[,splitvar], levels = levels(df[,splitvar])[i]) ### So that it contains only one level
        A <- as.data.frame(eval(parse(text=paste('complete(A,', paste(c(Grouping,Rep) ,collapse=','),')',sep=''))) )
      }
      else{
        B <- DF1[DF1[,splitvar]==levels(df[,splitvar])[i],]
        B[,splitvar] <- factor(B[,splitvar], levels = levels(df[,splitvar])[i]) ### So that it contains only one level
        B <- as.data.frame(eval(parse(text=paste('complete(B,', paste(c(Grouping,Rep) ,collapse=','),')',sep=''))) )
        A <- rbind(A,B)
      }
    }
  }
  else {
    A <- DF1
  }
  A$Proportion <- NA
  for (i in 1:length(A$N)){
    if (splitvar != ''){
      A$Proportion2[i] <- A$N[i]/sum(A$N[A[,splitvar]==A[,splitvar][i] & A[,Rep]==A[,Rep][i]], na.rm = T)
    }
    else {
      A$Proportion2[i] <- A$N[i]/sum(A$N[A[,Rep]==A[,Rep][i]], na.rm = T)
    }
  }
  A$Proportion2 <- ifelse(is.na(A$Proportion2)==T,0,A$Proportion2)
  #print(A)
  #print(str(A))
  A[,xvar] <- factor(A[,xvar], levels = levels(df[,xvar]))
  DF2 <- ddply(A, Grouping, summarise, Proportion = mean(Proportion2), P_sd = sd(Proportion2), N_group = length(Proportion2))
  #print(DF2)
  DF2$sdlow <- ifelse(DF2$Proportion - DF2$P_sd < 0, DF2$Proportion, DF2$P_sd)
  DF2$sdupper <- ifelse(DF2$Proportion + DF2$P_sd > 1, 1-DF2$Proportion, DF2$P_sd)
  DF2[DF2 == 0] <- NA
  DF2 <- as.data.frame(eval(parse(text=paste('complete(DF2,', paste(Grouping,collapse=','),')',sep=''))) )
  DF2$N <- nest_prop(df,Grouping,xvar,splitvar='')$N
  DF2$se <- DF2$P_sd/sqrt(DF2$N_group)
  DF2$CI <- qnorm(0.975)*DF2$se
  DF2$CI_lower <- ifelse(DF2$Proportion - DF2$CI < 0, DF2$Proportion, DF2$CI)
  DF2$CI_upper <- ifelse(DF2$Proportion + DF2$CI > 1, 1-DF2$Proportion, DF2$CI)
  return(DF2)
}


# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}


# Backcross generation data
Mating <- read.csv(".../Data_June-July2015.csv")
Mating$Time[Mating$Time == ""] <- "17:00"
Mating$RelTime <- RelativeTime(Mating$Time, "16:30") # in minutes
Mating$RelTime2 <- Mating$RelTime/60 # in hours (continuous)
Mating$RelTime3 <- floor(Mating$RelTime2*2)/2 # Rounded down to half hour
Mating1 <- subset(Mating, Mating$Callus != "")
Mating1$Callus1 <- TransformCallus1(Mating1,"Callus")
Mating1$Callus2 <- TransformCallus2(Mating1,"Callus")
Mating1$Callus3 <- TransformCallus3(Mating1,"Callus")
Mating1$SimpleTime <- TransformMating(Mating1,"RelTime2")
Mating1$SimpleTime <- factor(Mating1$SimpleTime, levels = c("Day-neo","Intermediate","Dusk-try"))
Mating1$Callus3 <- factor(Mating1$Callus3, levels = c("bbbb","bbbY","bbYY","bYYY","YYYY")) ### Reorder the callus

New_Callus <- c("bbbb","bbbY","bbYY","bYYY","YYYY")

### Table of simplified phenotype
Mating1$Start <- ifelse(Mating1$Filter == "Yes" & Mating1$Dusk == "Day", "Day", ifelse(Mating1$Filter == "No" & Mating1$Dusk == "Dusk", "Day","Dusk"))
Callus_count <- aggregate(Mating1$Callus3, by=list(Mating1$Callus3,Mating1$Sex), FUN=length)

# F1 data
F1_mating <- read.csv(".../F1_Mating_April2015.csv")
F1_callus <- read.csv(".../F1_callus.csv")
F1_callus$Class2 <- ifelse(F1_callus$Class == 4, "try",
                           ifelse(F1_callus$Class == 0, "neo", F1_callus$Class))
F1_callus$Class <- New_Callus[as.numeric(F1_callus$Class)+1]
F1_mating$Time2 <- ifelse(F1_mating$Time == 1, "Day-neo",
                          ifelse(F1_mating$Time == 4, "Dusk-try", "Intermediate"))

# F25 data
# Profile of mating of males (before selection)
Male25 <- read.csv(".../Males_mating_time.csv")
#length(subset(Male25$Type, Male25$Type =="Drowned")) # 16 out of 388 are not mating data
M_25 <- subset(Male25, Male25$Type != "Drowned")
#length(M_25$Line)

M_25$Time <- as.character(M_25$Time) #Not factor, otherwise the time will be converted into order
M_25$Time <- ifelse(M_25$Time == "", ifelse(M_25$Remarks == "Dusk", "17:30",""), M_25$Time) # All records as dusk placed at mid-point 17:30

# Create a list to see how many times particular individual had mated over the two days
Freq_mating <- aggregate(M_25$Concatenate_ID, by = list(M_25$Concatenate_ID), FUN = length)
Twice <- subset(Freq_mating, Freq_mating$x == 2)

Twice_M <- M_25[M_25$Concatenate_ID %in% Twice$Group.1,] #Only keep individuals that have two data points
Once_M <- M_25[!(M_25$Concatenate_ID %in% Twice$Group.1),]

# output first/second occurrence of each ID
Twice_M$Occurrence <- 0
A <- c()
for (i in 1:length(Twice_M$Concatenate_ID)) {
  if (as.character(Twice_M$Concatenate_ID[i]) %in% A){
    Twice_M$Occurrence[i] <- 2
    A[i] <- as.character(Twice_M$Concatenate_ID[i])
  }
  else {
    Twice_M$Occurrence[i] <- 1
    A[i] <- as.character(Twice_M$Concatenate_ID[i])
  }
}

First_M <- subset(Twice_M, Twice_M$Occurrence == 1)
First_M <- rbind(Once_M,Twice_M[,!(names(Twice_M)=="Occurrence")])
Second_M <- subset(Twice_M, Twice_M$Occurrence == 2)
names(Second_M)[7] <- "Time2"
Both_M <- dplyr::full_join(First_M, Second_M[,c("Concatenate_ID","Time2")])

Both_M$RelTime <- RelativeTime(Both_M$Time, "17:00")
Both_M$RelTime2 <- Both_M$RelTime/60 # in hours (continuous)
Both_M$RelTime3 <- floor(Both_M$RelTime2*2)/2 # Rounded down to half hour
Both_M$RelTime4 <- floor(Both_M$RelTime2*4)/4 # Rounded down to quarter hour

Both_M$RelTime12 <- RelativeTime(Both_M$Time2, "17:00")
Both_M$RelTime22 <- Both_M$RelTime12/60 # in hours (continuous)
Both_M$RelTime32 <- floor(Both_M$RelTime22*2)/2 # Rounded down to half hour
Both_M$RelTime42 <- floor(Both_M$RelTime22*4)/4 # Rounded down to quarter hour

Both_M$EarlyR <- pmin(Both_M$RelTime3, Both_M$RelTime32, na.rm = T) 
Both_M$LateR <- pmax(Both_M$RelTime3, Both_M$RelTime32, na.rm = F) 

Both_M$EarlyR1 <- pmin(Both_M$RelTime2, Both_M$RelTime22, na.rm = T) 
Both_M$LateR1 <- pmax(Both_M$RelTime2, Both_M$RelTime22, na.rm = F) 

Both_M$EarlyR2 <- pmin(Both_M$RelTime4, Both_M$RelTime42, na.rm = T) 
Both_M$LateR2 <- pmax(Both_M$RelTime4, Both_M$RelTime42, na.rm = F)



F25Mat <- Both_M[Both_M$Line != "S06",]

length(F25Mat[F25Mat$EarlyR1>=0 & is.na(F25Mat$EarlyR1) != T,]$EarlyR1)
length(F25Mat[as.numeric(F25Mat$EarlyR1)>=0 & as.numeric(F25Mat$LateR1)>=0 & is.na(F25Mat$LateR1) != T,]$EarlyR1)


# Might not be used
# F25 males after selection... and re-mating with hybrid females
Male25DayDusk <- read.csv("C:/Users/yea017/Dropbox/CSIRO/Research related/Mating_time/Mating_timeF25-26/F25_SelectedDayDusk.csv")
Male25DayDusk$RelTime <- RelativeTime(Male25DayDusk$Time, "17:00")
Male25DayDusk$RelTime2 <- Male25DayDusk$RelTime/60 # in hours (continuous)
Male25DayDusk$RelTime3 <- floor(Male25DayDusk$RelTime2*2)/2 # Rounded down to half hour
Male25DayDusk$RelTime4 <- floor(Male25DayDusk$RelTime2*4)/4 # Rounded down to quarter hour

## F25 callus colour with mating group
F25_callus <- read.csv(".../F25_callus.csv")
F25_callus <- subset(F25_callus, F25_callus$Updated_callus !="")
F25_callus$Callus1 <- TransformCallus1(F25_callus,"Updated_callus")
F25_callus$Callus2 <- TransformCallus2(F25_callus,"Updated_callus")
F25_callus$Callus3 <- TransformCallus3(F25_callus,"Updated_callus")

## Might not be used because no heritability estimates
### Female25 mating time
Fem25 <- read.csv(".../FemaleF25.csv")
Fem25$RelTime <- RelativeTime(Fem25$Time, "17:00")
Fem25$RelTime2 <- Fem25$RelTime/60


# F26 data (all of the males two matings, including post selectiong mating - no number identifier)
Male26Day <- read.csv(".../F26_SelectedDay.csv")
Male26Day$Callus3 <- Male26Day$Callus - 1
Male26Day$RelTime <- RelativeTime(Male26Day$Time_Start, "17:00")
Male26Day$RelTime2 <- Male26Day$RelTime/60 # in hours (continuous)
Male26Day$RelTime3 <- floor(Male26Day$RelTime2*2)/2 # Rounded down to half hour

Male26Day$RelTime12 <- RelativeTime(Male26Day$Time_Start2, "17:00")
Male26Day$RelTime22 <- Male26Day$RelTime12/60 # in hours (continuous)
Male26Day$RelTime32 <- floor(Male26Day$RelTime22*2)/2 # Rounded down to half hour

Male26Day$EarlyR <- pmin(Male26Day$RelTime3, Male26Day$RelTime32, na.rm = T) 
Male26Day$LateR <- pmax(Male26Day$RelTime3, Male26Day$RelTime32, na.rm = F) 


Male26Day$EarlyR1 <- pmin(Male26Day$RelTime2, Male26Day$RelTime22, na.rm = T) 
Male26Day$LateR1 <- pmax(Male26Day$RelTime2, Male26Day$RelTime22, na.rm = F) 

#Male26Day <- Male26Day[is.na(Male26Day$RelTime) == F,]

Male26Day1 <- subset(Male26Day, is.na(Male26Day$EarlyR1)== F)
Male26Early2 <- subset(Male26Day1, Male26Day1$Line %in% c("1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","2_1","2_2","2_3","2_4","2_5","2_6","2_7","3_1","3_2","3_3"))
Male26Early2$Line <- factor(Male26Early2$Line, levels = c("1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","2_1","2_2","2_3","2_4","2_5","2_6","2_7","3_1","3_2","3_3") )

# F26 Late day selected
Male26Late <- read.csv(".../F26_LateDay_Dusk.csv")
Male26Late$RelTime <- RelativeTime(Male26Late$Time_Start, "17:00")
Male26Late$RelTime2 <- Male26Late$RelTime/60 # in hours (continuous)
Male26Late$RelTime3 <- floor(Male26Late$RelTime2*4)/4 # Rounded down to quarter hour
Male26Dusk <- subset(Male26Late, Male26Late$Line == "Dusk") ### Only dusk selected lines

# Callus frequency F26 - to compare callus frequency between selected and unselected lines
General26 <- read.csv(".../F26_General_Data.csv")
General26$Cal1 <- ifelse(General26$Callus >= 4, "Yellow",General26$Callus)
Gen26Y <- subset(General26, General26$Cal1 == "Yellow")



# F35 data
Final35 <- read.csv(".../F35_experiment.csv")
Final35$Time1 <- ifelse(Final35$Dusk1 == "Yes","17:00","")
Final35$RelTime1 <- RelativeTime(Final35$Time1, "16:30")
Final35$RelTime1B <- Final35$RelTime1/60
Final35$RelTime1C <- floor(Final35$RelTime1B*2)/2

Final35$RelDay1 <- RelativeTime(Final35$Day1, "16:30")
Final35$RelDay1B <- Final35$RelDay1/60
Final35$RelDay1C <- floor(Final35$RelDay1B*2)/2

Final35$RelDay2 <- RelativeTime(Final35$Day2, "16:30")
Final35$RelDay2B <- Final35$RelDay2/60
Final35$RelDay2C <- floor(Final35$RelDay2B*2)/2

Final35$RelDay3 <- RelativeTime(Final35$Day3, "16:30")
Final35$RelDay3B <- Final35$RelDay3/60
Final35$RelDay3C <- floor(Final35$RelDay3B*2)/2


Final35$Time2 <- ifelse(Final35$Dusk2 != "", "17:00","")
Final35$RelTime2 <- RelativeTime(Final35$Time2, "16:30")
Final35$RelTime2B <- Final35$RelTime2/60
Final35$RelTime2C <- floor(Final35$RelTime2B*2)/2
Final35$MultiMate <- 5 - (is.na(Final35$RelTime1C) + is.na(Final35$RelTime2C) + is.na(Final35$RelDay1C) +  is.na(Final35$RelDay2C) + is.na(Final35$RelDay3C))
Final35$EarlyR <- pmin(Final35$RelTime1C, Final35$RelTime2C, Final35$RelDay1C,  Final35$RelDay2C, Final35$RelDay3C, na.rm = T)
Final35$LateR <- ifelse(Final35$MultiMate > 1,pmax(Final35$RelTime1C, Final35$RelTime2C, Final35$RelDay1C,  Final35$RelDay2C, Final35$RelDay3C, na.rm = T),NA)
Final35$EarlyR2 <- pmin(Final35$RelTime1B, Final35$RelTime2B, Final35$RelDay1B,  Final35$RelDay2B, Final35$RelDay3B, na.rm = T)
Final35$LateR2 <- ifelse(Final35$MultiMate > 1,pmax(Final35$RelTime1B, Final35$RelTime2B, Final35$RelDay1B,  Final35$RelDay2B, Final35$RelDay3B, na.rm = T),NA)


Final35$Range <- Final35$LateR - Final35$EarlyR
Final35$Range2 <- Final35$LateR2 - Final35$EarlyR2
Final35$SimpleTime <- TransformMating(Final35,"EarlyR")
Final35$SimpleTime2 <- TransformMating(Final35,"LateR")

Dusk35 <- subset(Final35, Final35$Start == "Dusk")
Day35 <- subset(Final35, Final35$Start == "DAY")


Day35M <- subset(Day35 ,Day35$Sex == "Male")
Day35F <- subset(Day35, Day35$Sex == "Female")
Dusk35M <- subset(Dusk35, Dusk35$Sex == "Male")
Dusk35F <- subset(Dusk35, Dusk35$Sex == "Female")


Final35$Callus1 <- TransformCallus1(Final35,"Yellow")
Final35$Callus2 <- TransformCallus2(Final35,"Yellow")
Final35$Callus3 <- as.factor(TransformCallus3(Final35,"Yellow"))
Final35$Callus3 <- factor(Final35$Callus3, levels=c("bbbb","bbbY","bbYY","bYYY","YYYY"))




# Create full early mating time dataset
F25_35_mating <- data.frame(Line = c(as.character(Both_M$Line),as.character(Male26Day$Line), substring(Mass28$Line,2), as.character(Mass34$Line), substring(as.character(Day35M$Line),2), substring(as.character(Dusk35M$Line),2)),
                            Start = c(rep(NA,372+800+91+126+420)),
                            Generation = c(rep(25,372),rep(26,800),rep(28,91),rep(34,126),rep(35,420)),
                            Early = c(Both_M$EarlyR,Male26Day$EarlyR,Mass28$RelTime3,Mass34$RelTime3,Day35M$EarlyR,Dusk35M$EarlyR),
                            Late = c(Both_M$LateR,Male26Day$LateR,rep(NA,91+126),Day35M$LateR,Dusk35M$LateR),
                            Early2 = c(Both_M$EarlyR1,Male26Day$EarlyR1,Mass28$RelTime3,Mass34$RelTime3,Day35M$EarlyR2,Dusk35M$EarlyR2),
                            Late2 = c(Both_M$LateR1,Male26Day$LateR1,rep(NA,91+126),Day35M$LateR2,Dusk35M$LateR2))

F25_35_mating$SimpleTime <- TransformMating(F25_35_mating,"Early")
F25_35_mating$Line2 <- ifelse(F25_35_mating$Generation == 25, substring(F25_35_mating$Line,2),
                              ifelse(F25_35_mating$Generation == 26,ifelse(F25_35_mating$Line %in% c("2_1","2_2","2_3","2_4","2_5","2_6","2_7"), 2,
                                                                           as.character(F25_35_mating$Line)),as.character(F25_35_mating$Line)))

F25_35_Twice <- subset(F25_35_mating, is.na(F25_35_mating$Late) == F & F25_35_mating$Line2 %in% c("1","2","1_1","1_5") )
F25_35_Twice$Range <- F25_35_Twice$Late - F25_35_Twice$Early
F25_35_Twice$Range2 <- F25_35_Twice$Late2 - F25_35_Twice$Early2


library(tiff)

neo0 <- rasterGrob(readTIFF(".../BC(Neo)_177_0.tif"),interpolate=T)
neo1 <- rasterGrob(readTIFF(".../BC(Neo)_035_25.tif"),interpolate=T)
neo2 <- rasterGrob(readTIFF(".../BC(S06)_051_50.tif"),interpolate=T)
neo3 <- rasterGrob(readTIFF(".../BC(Neo)_029_85.tif"),interpolate=T)
neo4 <- rasterGrob(readTIFF(".../BC(Neo)_063_100.tif"),interpolate=T)



# F1
F1_mating$Time2 <- factor(F1_mating$Time2, levels = c("Day-neo","Intermediate","Dusk-try"))
F1_callus$Class <- factor(F1_callus$Class, levels = c("bbbb","bbbY","bbYY","bYYY","YYYY"))

try1 <- subset(F1_mating,F1_mating$Line %in% c("try"))
neo1A <- subset(F1_mating,F1_mating$Line %in% c("neo"))
F1_1 <- subset(F1_mating,F1_mating$Line %in% c("F1(neo)","F1(try)"))
B1_1 <- subset(Mating1,Mating1$Start=="Day" & Mating1$Mate2 %in% c('neo','try'))
B1_1$Mate2 <- factor(B1_1$Mate2, levels=c('neo','try')) 
B1_2 <- subset(Mating1,Mating1$Start=="Dusk" & Mating1$Mate2 %in% c('neo','try'))
B1_2$Mate2 <- factor(B1_2$Mate2, levels=c('neo','try')) 

B1_All <- subset(Mating1,Mating1$Mate2 %in% c('neo','try'))
B1_All$Mate2 <- factor(B1_All$Mate2, levels=c('neo','try')) 


try_graph <- nest_prop2(try1,c("Sex","Mate2","Time2"),"Time2","Sex","Rep")
neo_graph <- nest_prop2(neo1A,c("Sex","Mate2","Time2"),"Time2","Sex","Rep")
F1_graph <- nest_prop2(F1_1,c("Sex","Mate2","Time2"),"Time2","Sex","Rep")
B1_graph <- nest_prop2(B1_1,c("Sex","Mate2","SimpleTime"),"SimpleTime","Sex","Rep")
B1_Dusk_graph <- nest_prop2(B1_2,c("Sex","Mate2","SimpleTime"),"SimpleTime","Sex","Rep")
B1_Day_graph <- nest_prop2(B1_1,c("Sex","Mate2","SimpleTime"),"SimpleTime","Sex","Rep")
B1_All_graph <- nest_prop2(B1_All,c("Sex","Mate2","SimpleTime"),"SimpleTime","Sex","Rep")


ggplot(B1_Dusk_graph, aes(y = Proportion,x = SimpleTime,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none")+ ylim(0,1.1)
ggplot(B1_All_graph, aes(y = Proportion,x = SimpleTime,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none")+ ylim(0,1.1)


neo_FMNT <- nest_prop2(neo1A,c("Sex","Mate2"),"Sex","Sex","Rep")
try_FMNT <- nest_prop2(try1,c("Sex","Mate2"),"Sex","Sex","Rep")
F1_FMNT <- nest_prop2(F1_1,c("Sex","Mate2"),"Sex","Sex","Rep")
B1_FMNT <- nest_prop2(B1_All,c("Sex","Mate2"),"Sex","Sex","Rep")

B1_All$SexS <- interaction(B1_All$Sex, B1_All$Start)
B1_All$SexSD <- interaction(B1_All$Sex, B1_All$Start, B1_All$Dusk)
nest_prop2(B1_All,c("SexS","SimpleTime"),"SimpleTime","SexS","Rep")
B1_DDFMNT <- nest_prop2(B1_All,c("SexSD","Mate2"),"Mate2","SexSD","Rep")

# Male - neo vs try female mate
t.test2(try_FMNT$Proportion[3],try_FMNT$Proportion[4],try_FMNT$P_sd[3],try_FMNT$P_sd[4],
        try_FMNT$N[3],try_FMNT$N[4],0,equal.variance = T) # No evidence of difference in mating to either parentals in B. tryoni males
t.test2(F1_FMNT$Proportion[3],F1_FMNT$Proportion[4],F1_FMNT$P_sd[3],F1_FMNT$P_sd[4],
        F1_FMNT$N[3],F1_FMNT$N[4],0,equal.variance = T) # No evidence for mating more with one parental 
t.test2(B1_FMNT$Proportion[3],B1_FMNT$Proportion[4],B1_FMNT$P_sd[3],B1_FMNT$P_sd[4],
        B1_FMNT$N[3],B1_FMNT$N[4],0,equal.variance = T) # Evidence for mating preferentially with neo females

t.test2(B1_DDFMNT$Proportion[15],B1_DDFMNT$Proportion[16],B1_DDFMNT$P_sd[15],B1_DDFMNT$P_sd[16],
        B1_DDFMNT$N[15],B1_DDFMNT$N[16],0,equal.variance = T) # Evidence for mating preferentially with neo females in dusk males if start during dusk
t.test2(B1_DDFMNT$Proportion[7],B1_DDFMNT$Proportion[8],B1_DDFMNT$P_sd[7],B1_DDFMNT$P_sd[8],
        B1_DDFMNT$N[7],B1_DDFMNT$N[8],0,equal.variance = F) # Evidence for mating preferentially with try females in dusk males if start during the day



# Female - neo vs try male mate
t.test2(neo_FMNT$Proportion[1],neo_FMNT$Proportion[2],neo_FMNT$P_sd[1],neo_FMNT$P_sd[2],
        neo_FMNT$N[1],neo_FMNT$N[2],0,equal.variance = F)


t.test2(B1_FMNT$Proportion[1],B1_FMNT$Proportion[2],B1_FMNT$P_sd[1],B1_FMNT$P_sd[2],
        B1_FMNT$N[1],B1_FMNT$N[2],0,equal.variance = T)


#Supplementary observations

# Discordant mating - if starting at dusk, less mating to tryoni than when starting at day
# More mating to neo
######## NEED TO CHANGE EVERYTHING TO FACTOR
#nest_prop2(B1_All,c("Sex","Mate2","Start"),"Start","Sex","Rep")



# Female
wilcox.test(subset(try1, try1$Sex=="Female")$Time, subset(try1, try1$Sex=="Male")$Time, conf.int = T)
wilcox.test(subset(neo1A, neo1A$Sex=="Female")$Time, subset(neo1A, neo1A$Sex=="Male")$Time)
wilcox.test(subset(F1_1, F1_1$Sex=="Female")$Time, subset(F1_1, F1_1$Sex=="Male")$Time, conf.int = T)
wilcox.test(subset(B1_All, B1_All$Sex=="Female")$RelTime, subset(B1_All, B1_All$Sex=="Male")$RelTime)

B1_calDay <- nest_prop2(B1_1,c("Sex","Mate2","Callus3"),"Callus3","Sex","Rep")
B1_calDusk <- nest_prop2(B1_2,c("Sex","Mate2","Callus3"),"Callus3","Sex","Rep")
B1_calAll <- nest_prop2(B1_All,c("Sex","Mate2","Callus3"),"Callus3","Sex","Rep")

B1_All$SCal <- interaction(B1_All$Sex,B1_All$Callus3)
B1_ScalAll <- nest_prop2(B1_All,c("SCal","Mate2"),"Mate2","SCal","Rep")
B1_ScalAll
t.test2(B1_ScalAll$Proportion[1],B1_ScalAll$Proportion[2],B1_ScalAll$P_sd[1],B1_ScalAll$P_sd[2],B1_ScalAll$N[1],B1_ScalAll$N[2],
        0, equal.variance = T)


B1_calAll2 <- nest_prop2(B1_All,c("Sex","Callus3"),"Callus3","Sex","Rep")

# proportion of callus colour in overall B1
prop.table(table(B1_All$Callus3))


### Manually create table - counting total number of callus and mating... IMPORTANT... callus and mating variance not clear in males
levels(interaction(Mating1$Dusk,Mating1$Start))
nest_prop2(Mating1, c("Start","Sex","Callus3"),"Callus3","Sex","Rep")
aggregate(Mating1$SimpleTime, by=list(Mating1$SimpleTime, Mating1$Sex, Mating1$Start), FUN=length)
aggregate(Mating1$Callus3, by=list(Mating1$Callus3, Mating1$Sex, Mating1$Start), FUN=length)


aggregate(Mating1$Callus3, by=list(Mating1$Callus3, Mating1$SimpleTime, Mating1$Sex, Mating1$Start), FUN=length)
ggplot(Mating1, aes(x=SimpleTime,y=..prop.., group=1)) + geom_bar() + facet_grid(rows=vars(Start), cols=vars(Sex))
ggplot(Mating1, aes(x=Callus3,y=..prop.., group=1)) + geom_bar() + facet_grid(rows=vars(Start), cols=vars(Sex))


### Split by day and dusk mating (day = day + intermediate)
B1_DAY <- subset(Mating1, Mating1$Mate2 %in% c('neo','try') & Mating1$SimpleTime == "Day-neo")
B1_DAY$Mate2 <- factor(B1_DAY$Mate2, levels=c('neo','try'))
B1_DAY$Callus4 <- as.numeric(str_count(B1_DAY$Callus3,"Y"))
B1_DAY$MateTime <- "Day" 

B1_Int <- subset(Mating1, Mating1$Mate2 %in% c('neo','try') & Mating1$SimpleTime == "Intermediate")
B1_Int$Mate2 <- factor(B1_Int$Mate2, levels=c('neo','try'))
B1_Int$Callus4 <- as.numeric(str_count(B1_Int$Callus3,"Y"))
B1_Int$MateTime <- "Day"

B1_DUSK <- subset(Mating1, Mating1$Mate2 %in% c('neo','try') & Mating1$SimpleTime == "Dusk-try")
B1_DUSK$Mate2 <- factor(B1_DUSK$Mate2, levels=c('neo','try'))
B1_DUSK$Callus4 <- as.numeric(str_count(B1_DUSK$Callus3,"Y"))
B1_DUSK$MateTime <- "Dusk"



B1_calMate <- rbind(nest_prop2(B1_DAY,c("Sex","Callus3"),"Callus3","Sex","Rep"),nest_prop2(B1_Int,c("Sex","Callus3"),"Callus3","Sex","Rep"),nest_prop2(B1_DUSK,c("Sex","Callus3"),"Callus3","Sex","Rep") )
B1_calMate$Dusk <- c(rep("Day-neo",10),rep("Intermediate",10),rep("Dusk-try",10))
B1_calMate$Recom <- ifelse(B1_calMate$Dusk == "Day-neo" & B1_calMate$Callus3 %in% c("bYYY","YYYY"),
                           "Recombinant", ifelse(B1_calMate$Dusk == "Dusk-try" & B1_calMate$Callus3 %in% c("bbbb"), "Recombinant", 
                                                 ifelse(B1_calMate$Dusk == "Day-neo" & B1_calMate$Callus3 %in% c("bbbb"),"Parental",
                                                        ifelse(B1_calMate$Dusk == "Dusk-try" & B1_calMate$Callus3 %in% c("bYYY","YYYY"),"Parental","Other")))) 


B1_DAY2 <- subset(Mating1, Mating1$Mate2 %in% c('neo','try') & Mating1$Dusk == "Day" & Mating1$Start == "Dusk")
B1_DAY2$Mate2 <- factor(B1_DAY2$Mate2, levels=c('neo','try'))
B1_DAY2$Callus4 <- as.numeric(str_count(B1_DAY2$Callus3, "Y"))
B1_DAY2$MateTime <- "Day"

B1_DUSK2 <- subset(Mating1, Mating1$Mate2 %in% c('neo','try') & Mating1$Dusk == "Dusk" & Mating1$Start == "Dusk")
B1_DUSK2$Mate2 <- factor(B1_DUSK2$Mate2, levels=c('neo','try'))
B1_DUSK2$Callus4 <- as.numeric(str_count(B1_DUSK2$Callus3, "Y"))
B1_DUSK2$MateTime <- "Dusk"
B1_calMate2 <- rbind(nest_prop2(B1_DAY2,c("Sex","Callus3"),"Callus3","Sex","Rep"),nest_prop2(B1_DUSK2,c("Sex","Callus3"),"Callus3","Sex","Rep") )
B1_calMate2$Dusk <- c(rep("Day",10),rep("Dusk",10))


wilcox.test(subset(B1_DAY, B1_DAY$Sex=="Female")$Callus4, subset(B1_DUSK, B1_DUSK$Sex=="Female")$Callus4)
wilcox.test(subset(B1_DAY, B1_DAY$Sex=="Male")$Callus4, subset(B1_DUSK, B1_DUSK$Sex=="Male")$Callus4)

### Proportion of callus for each mating class
prop.table(table(B1_DAY$Callus3))


### This chisq test for linkage section
chisq.test(subset(rbind(B1_DAY,B1_DUSK), rbind(B1_DAY,B1_DUSK)$Sex=="Female")$Dusk,subset(rbind(B1_DAY,B1_DUSK),rbind(B1_DAY,B1_DUSK)$Sex=="Female")$Callus3)
chisq.test(subset(rbind(B1_DAY,B1_DUSK), rbind(B1_DAY,B1_DUSK)$Sex=="Male")$Dusk,subset(rbind(B1_DAY,B1_DUSK),rbind(B1_DAY,B1_DUSK)$Sex=="Male")$Callus3)

wilcox.test(subset(B1_DAY, B1_DAY$Sex=="Female")$Callus4, subset(B1_DAY, B1_DAY$Sex=="Male")$Callus4)
wilcox.test(subset(B1_DUSK, B1_DUSK$Sex=="Female")$Callus4, subset(B1_DUSK, B1_DUSK$Sex=="Male")$Callus4)
chisq.test(subset(rbind(B1_DAY,B1_DUSK), rbind(B1_DAY,B1_DUSK)$Dusk=="Day")$Sex,subset(rbind(B1_DAY,B1_DUSK),rbind(B1_DAY,B1_DUSK)$Dusk=="Day")$Callus3)
chisq.test(subset(rbind(B1_DAY,B1_DUSK), rbind(B1_DAY,B1_DUSK)$Dusk=="Dusk")$Sex,subset(rbind(B1_DAY,B1_DUSK),rbind(B1_DAY,B1_DUSK)$Dusk=="Dusk")$Callus3)



ggplot(try_graph, aes(y = Proportion,x = Time2,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N),stat='identity', vjust = -1,position=position_dodge(width = 0.5)) + geom_errorbar(aes(ymin=Proportion-P_sd,ymax = Proportion+P_sd), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Count", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1)


F1_mateplot <- ggplot(F1_graph, aes(y = Proportion,x = Time2,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1.1)
ggplot(subset(F1_mating,F1_mating$Line %in% c("neo")), aes(x = Time2)) + geom_bar(aes(y = ..count..,fill = Mate2)) + facet_grid(cols=vars(Sex))  + labs(y= "Count", x = "Mating time profile") + guides(fill=guide_legend(title="Mated to"))
try_mate <- ggplot(try_graph, aes(y = Proportion,x = Time2,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+ifelse(is.na(CI_upper)==T,0,CI_upper) ),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1.1)
neo_mate <- ggplot(neo_graph, aes(y = Proportion,x = Time2,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1.1)


Cal_top <- ggplot() + xlim(0,1) + ylim(0,0.5) + 
  annotation_custom(neo0, xmin=0,xmax =0.18,ymin=0.08, ymax=0.40) + 
  annotation_custom(neo1, xmin=0.2,xmax =0.38,ymin=0.08, ymax=0.40) + 
  annotation_custom(neo2, xmin=0.4,xmax =0.58,ymin=0.08, ymax=0.40) + 
  annotation_custom(neo3, xmin=0.6,xmax =0.78,ymin=0.08, ymax=0.40) + 
  annotation_custom(neo4, xmin=0.8,xmax =0.98,ymin=0.08, ymax=0.40) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

F1_calplot <- ggplot(F1_callus, aes(x = Class)) + geom_bar(aes(y = ..prop.., group = 1)) + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE) + ylim(0,0.7) + ggtitle("F1 mixed sexes, N = 40") + 
  theme(plot.title = element_text(size = 10))


# Backcross
B1_day <- ggplot(B1_graph, aes(y = Proportion,x = SimpleTime,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1.1)
B1_dusk <- ggplot(B1_Dusk_graph, aes(y = Proportion,x = SimpleTime,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,1.1)
B1_all <- ggplot(B1_All_graph, aes(y = Proportion,x = SimpleTime,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none")+ ylim(0,1.1)
ggplot(subset(Mating1,Mating1$Start=="Dusk"), aes(x = SimpleTime)) + geom_bar(aes(y = ..count..,fill=Mate2)) + facet_grid(cols=vars(Sex)) + labs(y= "Proportion", x = "Mating time profile")# + ylim(0,1)
ggplot(Mating1, aes(x = Callus3)) + geom_bar(aes(y = ..prop.., group = 1)) + labs(y= "Proportion", x = "Callus") + ylim(0,1) + ggtitle("Backcross, n = 3068")
length(subset(Mating1,Mating1$Start=="Day")$Callus3)

B1_dayCal <- ggplot(B1_calDay, aes(y = Proportion,x = Callus3,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,0.75)
B1_duskCal <- ggplot(B1_calDusk, aes(y = Proportion,x = Callus3,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y = Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none",axis.title.x=element_blank())+ ylim(0,0.75)
Sexlab <- c(paste("Female, N =",sum(subset(B1_calAll$N, B1_calAll$Sex == "Female"))), paste("Male, N =", sum(subset(B1_calAll$N, B1_calAll$Sex == "Male"))))
names(Sexlab) <- c("Female", "Male")
B1_allCal <- ggplot(B1_calAll) + geom_bar(aes(y = Proportion,x = Callus3,fill = Mate2),stat='identity', position='dodge') + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper, x=Callus3, fill=Mate2), width=0.2, position=position_dodge(width = 1)) + 
  facet_grid(~Sex, labeller = labeller(Sex = Sexlab)) + scale_fill_grey()  + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + theme(legend.position = "none")+ ylim(0,0.7) #+ geom_text(aes(label=N, y = Proportion+CI_upper, x=Callus3, fill=Mate2),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) 

B1_plain <- ggplot(B1_All, aes(x = Callus3)) + geom_bar(aes(y = ..prop.., group = 1)) + labs(y= "Proportion", x = "Callus") + ylim(0,1) + ggtitle("Backcross, n = 3057")



MatB1 <- subset(Mating1, Mating1$Mate2 %in% c("neo","try"))
MatB1$Mate2 <- factor(as.factor(MatB1$Mate2), levels = c("neo","try"))



legend1 <- get_legend(
  # create some space to the left of the legend
  ggplot(F1_graph, aes(y = Proportion,x = Time2,fill = Mate2)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N),stat='identity', vjust = -0.15,position=position_dodge(width = 0.5)) + geom_errorbar(aes(ymin=Proportion-P_sd,ymax = Proportion+P_sd), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Count", x = "Mating time profile") + scale_x_discrete(drop=FALSE) + guides(fill=guide_legend(title="Mated to")) + ylim(0,1) + guides(fill=guide_legend(title="Mated to")) + theme(legend.box.margin = margin(0, 0, 0, 12))
)

### Figure 2 - try 1500 (W) by 800 (H) (10inch x8inch)
jpeg(".../Figure2.jpg",width = 5800, height = 4000, res = 600)
ggdraw() +
  #draw_label("Mating time profile", fontface='bold', x = 0.5, y = 0.97)+
  draw_plot(neo_mate, x = 0.02, y = 0.62, width = .48, height = 0.28) +
  draw_plot(try_mate, x = .52, y = 0.62, width = .48, height = 0.28) +
  draw_plot(F1_mateplot, x = 0.25, y = 0.32, width = 0.5, height = 0.28) +
  draw_plot(plot_grid(NULL,legend1), x = 0.1) +
  draw_plot(B1_all, x = 0.25, y = 0, width = 0.5, height = 0.3) +
  draw_plot_label(label = c("P","Bactrocera neohumeralis", "Bactrocera tryoni", "F1","B1"), size = 15,
                  x = c(0,0.05,0.6,0,0), y = c(0.85,0.94, 0.94, 0.55, 0.25), fontface = c("bold",rep("bold.italic",2), rep("bold",2)))
dev.off()
###
groupwiseMedian(RelTime2 ~ 1, data = subset(B1_All, B1_All$Sex == "Male"), conf = 0.95, R = 5000, percentile = TRUE, bca = FALSE, basic = FALSE, normal = FALSE, wilcox = FALSE, digits = 3)
groupwiseMedian(RelTime2 ~ 1, data = subset(B1_All, B1_All$Sex == "Female"), conf = 0.95, R = 5000, percentile = TRUE, bca = FALSE, basic = FALSE, normal = FALSE, wilcox = FALSE, digits = 3)
# Supplementary 1
# Compare dusk start vs day start


# Figure 3
jpeg(".../Figure3.jpg",width = 4800, height = 4000, res = 600)
ggdraw() +
  #draw_label("Callus colour profile", fontface='bold', x = 0.5, y = 0.97)+
  draw_plot(Cal_top, x = 0.01, y = 0.7, width = .99, height = 0.28) +
  #draw_plot(try_mate, x = .525, y = 0.62, width = .475, height = 0.28) +
  draw_plot(F1_calplot, x = 0.25, y = 0.35, width = 0.5, height = 0.30) +
  draw_plot(plot_grid(NULL,legend1), x = 0.2, y = -0.35) +
  draw_plot(B1_allCal, x = 0.1, y = 0, width = 0.8, height = 0.3) +
  draw_plot_label(label = c("bbbb","bbbY","bbYY","bYYY","YYYY", "F1","B1"), size = 15,
                  x = c(0.08,0.26,0.44,0.62,0.79, 0,0), y = c(rep(0.75,5), 0.55, 0.25))
dev.off()
# Supplementary figure for female vs male callus colour
ggplot(B1_calAll2, aes(y = Proportion,x = Callus3)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y = Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey()  + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE) + ylim(0,0.75)



# Figure 4
jpeg(".../Figure4.jpg",width = 4000, height = 4000, res = 600)
B1_calMate$Dusk <- factor(B1_calMate$Dusk, levels = c("Day-neo","Intermediate","Dusk-try"))
B1_calMate$Perc <- paste(as.character(format(round(B1_calMate$Proportion*100,1)), nsmall=1),"%",sep="")
B1_plotcM <- ggplot(B1_calMate, aes(y = Proportion,x = Callus3, fill=Recom)) + geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex),rows = vars(Dusk)) + theme(legend.position = "none") +
  scale_fill_manual("legend", values = c("Other" = "dimgray", "Parental" = "Black", "Recombinant" = "Grey"))  + labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE)  + geom_text(aes(label=Perc, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + ylim(0,0.75)
B1_plotcM2 <- ggplot(B1_calMate, aes(y = Proportion,x = Callus3, fill=Dusk)) + geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1)) + facet_grid(cols=vars(Sex)) + scale_fill_grey() +
  labs(y= "Proportion", x = "Callus class") + scale_x_discrete(drop=FALSE)  + geom_text(aes(label=Perc, y=Proportion+CI_upper),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + ylim(0,0.75)

ggdraw() +
  draw_plot(B1_plotcM, x = 0, y = 0, width = 1, height = 1) +
  draw_plot_label(label = c("N = 532","N = 620","N = 498","N = 408","N = 152","N = 847"), size = 12,
                  x = c(rep(0.65,3),rep(0.20,3)), y = rep(c(0.95,0.65,0.35),2))
dev.off()

#plot_grid(plot1A, legend1)



# F25 only


Count_MatF25<- nest_prop(Both_M[Both_M$Line != "S06",],"EarlyR2","EarlyR2")
sum(Count_MatF25[as.numeric(Count_MatF25[,'EarlyR2']) >= 0,]$N)/sum(Count_MatF25$N)
sum(Count_MatF25[as.numeric(Count_MatF25[,'EarlyR2']) < 0,]$N)/sum(Count_MatF25$N)

mean(Both_M[Both_M$Line != "S06",]$EarlyR2) # -0.264 hours
mean(Both_M[Both_M$Line != "S06" & Both_M$EarlyR2 <= 0.5,]$EarlyR2)

median(rep(F2_mating$RelTime, F2_mating$Number))
median(Both_M[Both_M$Line != "S06",]$EarlyR2)
median(Both_M[Both_M$Line != "S06" & Both_M$EarlyR2 <= 0.5,]$EarlyR2)
wilcox.test(rep(F2_mating$RelTime, F2_mating$Number),Both_M[Both_M$Line != "S06",]$EarlyR2)




wilcox.test(rep(F2_mating$RelTime, F2_mating$Number),Both_M[Both_M$Line != "S06" & Both_M$EarlyR2 <= 0.5,]$EarlyR2)


F25_35_mating$Simple_Time <- factor(F25_35_mating$Simple_Time, levels = c("Day-neo","Intermediate","Dusk-try"))
p5a <- ggplot(data=F25_35_mating[F25_35_mating$Generation==25 & F25_35_mating$Line != "S06",], aes(x=Early)) + geom_bar(aes(y = (..prop..),group=1)) + labs(y="Proportion",x="Relative mating time, hours") + ylim(0,0.3) + xlim(-4,1)# + facet_grid(rows=vars(Line))

fig4dat <- Both_M[Both_M$Line != "S06",]
fig4dat <- data.frame(Line = "F25",EarlyR2 = as.numeric(levels(as.factor(fig4dat$EarlyR2))),Proportion=prop.table(table(Both_M[Both_M$Line != "S06",]$EarlyR2)))
fig4dat$MateShade <- ifelse(fig4dat$EarlyR2 >= 0, "Dusk maters",ifelse(fig4dat$EarlyR2 <= -1, "Early maters", "Other") )
fig4dat$MateShade <- factor(fig4dat$MateShade, levels = c("Dusk maters", "Early maters", "Other"))

p5b <- ggplot(data=fig4dat, aes(x=EarlyR2, fill = MateShade)) + geom_bar(aes(y=Proportion.Freq), stat='identity') + labs(y="Proportion",x="Relative mating time, hours") + ylim(0,0.3) + xlim(-4,1.25) + ggtitle("F25 mating time, N = 294") + theme(plot.title = element_text(size=12)) + 
  geom_vline(xintercept=c(-3,0),linetype='dashed') + annotate("text", x = c(-3.5,-1.5,0.75), y=0.275, label = c("Day-neo", "Intermediate", "Dusk-try"), size = 5) + 
  scale_fill_manual("legend", values = c("Other" = "dimgray", "Early maters" = "Grey", "Dusk maters" = "Black")) # + scale_fill_grey()
sum(Both_M[Both_M$Line != "S06",]$EarlyR2 >= 0)/length(Both_M[Both_M$Line != "S06",]$EarlyR2) # 46.6% dusk mater, 53.4% intermediate




#ggplot(data=F25_35_mating[F25_35_mating$Generation==25 & F25_35_mating$Line == "R1",], aes(x=Early)) + geom_bar(aes(y = (..prop..),group=1)) + labs(y="Proportion",x="Relative mating time, hours") + ylim(0,0.8) + xlim(-4,2) #+ facet_grid(rows=vars(Line))
#ggplot(data=F25_35_mating[F25_35_mating$Generation==25 & F25_35_mating$Line == "R2",], aes(x=Early)) + geom_bar(aes(y = (..prop..),group=1)) + labs(y="Proportion",x="Relative mating time, hours") + ylim(0,0.8) + xlim(-4,2) #+ facet_grid(rows=vars(Line))

F25_callus$Callus3 <- factor(F25_callus$Callus3, levels = c("bbbb","bbbY","bbYY","bYYY","YYYY"))
F25_callus$Rep <- ifelse(F25_callus$Replicate %in% c('R1','1'),1,
                         ifelse(F25_callus$Replicate %in% c('R2','2'),2,3))
F25_calprop <- nest_prop2(subset(F25_callus,F25_callus$Sex %in% c("Male","Female")), c("Callus3"),"Callus3",'',"Rep")
F25_calpropDay <- nest_prop2(subset(F25_callus,F25_callus$Additional.details =="Day_selected"), c("Callus3"),"Callus3",'',"Rep")
F25_calpropDusk <- nest_prop2(subset(F25_callus,F25_callus$Additional.details =="Dusk_selected"), c("Callus3"),"Callus3",'',"Rep")
F25_cal <- ggplot(F25_calprop, aes(x = Callus3, y = Proportion)) + geom_bar(stat='identity', position='dodge') + geom_text(aes(label=N, y = Proportion+ifelse(is.na(CI_upper)==T,0,CI_upper)),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1))+ ylim(0,0.6) + ggtitle("F25 callus, n = 138") + labs(x="Callus class")
F25_calDay <- ggplot(F25_calpropDay, aes(x = Callus3, y = Proportion)) + geom_bar(stat='identity', position='dodge', fill = "Grey") + geom_text(aes(label=N, y = Proportion+ifelse(is.na(CI_upper)==T,0,CI_upper)),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1))+ ylim(0,0.8) + ggtitle("Callus of early maters, N = 26") + theme(plot.title = element_text(size=12)) + labs(x="Callus class") + theme(legend.position = "none")
F25_calDusk <- ggplot(F25_calpropDusk, aes(x = Callus3, y = Proportion)) + geom_bar(stat='identity', position='dodge', fill = "Black") + geom_text(aes(label=N, y = Proportion+ifelse(is.na(CI_upper)==T,0,CI_upper)),stat='identity', vjust = -0.15,position=position_dodge(width = 1)) + geom_errorbar(aes(ymin=Proportion-CI_lower,ymax = Proportion+CI_upper), width=0.2, position=position_dodge(width = 1))+ ylim(0,0.8) + scale_x_discrete(drop=FALSE) + ggtitle("Callus of dusk maters, N = 85") + theme(plot.title = element_text(size=12)) + labs(x="Callus class") + theme(legend.position = "none")


F25_calpropDay$Callus <- c(0,1,2,3,4)
F25_calpropDusk$Callus <- c(0,1,2,3,4)
F25_calpropDusk$N <- ifelse(is.na(F25_calpropDusk$N)==T, 0, F25_calpropDusk$N)
wilcox.test(rep(F25_calpropDay$Callus,F25_calpropDay$N),rep(F25_calpropDusk$Callus,F25_calpropDusk$N) )
F25_DayDuskCal <- data.frame(Callus = c(rep(F25_calpropDay$Callus,F25_calpropDay$N),rep(F25_calpropDusk$Callus,F25_calpropDusk$N)),
                             Time = c(rep('Day',length(rep(F25_calpropDay$Callus,F25_calpropDay$N))), rep('Dusk',length(rep(F25_calpropDusk$Callus,F25_calpropDusk$N))) ))

chisq.test(F25_DayDuskCal$Callus,F25_DayDuskCal$Time)
fisher.test(F25_DayDuskCal$Callus,F25_DayDuskCal$Time)

# Figure 5
jpeg(".../Figure5.jpg",width = 5000, height = 4000, res = 600)
ggdraw() +
  #draw_plot(F2_mat, x = 0.05, y = 0.62, width = .475, height = 0.28) +
  #draw_plot(F2_cal, x = .525, y = 0.62, width = .475, height = 0.28) +
  draw_plot(p5b, x = 0.15, y = 0.52, width = 0.8, height = 0.48) +
  # draw_plot(F25_cal, x = 0.05, y = 0.5, width = 0.95, height = 0.25) +
  draw_plot(F25_calDay, x = 0.02, y = 0, width = 0.48, height = 0.5) +
  draw_plot(F25_calDusk, x = 0.52, y = 0, width = 0.48, height = 0.5) +
  draw_plot_label(label = c('A','B','C'), size = 15,
                  x = c(0.13,0.01,0.51), y = c(0.98,0.48,0.48), fontface = "bold")
dev.off()

Male26iso <- aggregate(Male26Day1$EarlyR1, by=list(Male26Day1$Line), FUN=mean)
names(Male26iso) <- c("Isofemale","Mean26")

F25LateDay <- subset(Male25DayDusk, Male25DayDusk$Line %in% c("R1_Backup","R2_Backup","R3_Backup") & Male25DayDusk$RelTime2 < 0)
F25Dusk <- subset(Male25DayDusk, Male25DayDusk$Select == "Dusk" & Male25DayDusk$Remarks == "Selected")

F25_mean <- rbind(data.frame(Group.1=Fem25[order(Fem25$Isofemale),]$Isofemale, x=Fem25[order(Fem25$Isofemale),]$RelTime2), aggregate(F25LateDay$RelTime2, by=list(F25LateDay$Line), FUN=mean), aggregate(F25Dusk$RelTime2, by=list(F25Dusk$Line), FUN=mean))

F26_mean <- rbind(aggregate(Male26Day$RelTime2, by=list(Male26Day$Line), FUN=mean)[1:18,],aggregate(Male26Late$RelTime2, by=list(interaction(Male26Late$Replicate, Male26Late$Line)), FUN=mean))
F26_length <- rbind(aggregate(Male26Day$RelTime2, by=list(Male26Day$Line), FUN=length)[1:18,], aggregate(Male26Late$RelTime2, by=list(interaction(Male26Late$Replicate, Male26Late$Line)), FUN=length))

F25_26 <- data.frame(Line = F25_mean$Group.1, Overall = mean(Both_M[Both_M$Line != "S06",]$RelTime2), F25 = F25_mean$x, F26 = F26_mean$x, N = F26_length$x)
F25_26$Heritability <- (F25_26$F26 - F25_26$Overall)/(F25_26$F25 - F25_26$Overall)

write.csv(F25_26,"C:/Users/yea017/Dropbox/CSIRO/Publication/Mating time choice/Table1.csv",row.names = F)



F25_35_selected <- F25_35_mating[F25_35_mating$Line2 %in% c('1','1_1','1_5','2'),]
F25_35_selected$Line2 <- ifelse(F25_35_selected$Line2 == '1', '1_1',ifelse(F25_35_selected$Line2 == 'R2', '2',as.character(F25_35_selected$Line2))) 
F25_1_5 <- subset(F25_35_selected, F25_35_selected$Line2 == "1_1" & F25_35_selected$Generation == 25)
F25_1_5$Line2 <- "1_5"

F25_35_selected <- rbind(F25_1_5,F25_35_selected)
F25_35_selected$Selection <- ifelse(F25_35_selected$Generation %in% c(25,26), "Early selected","No selection")
F25_35_selected$Generation <- factor(F25_35_selected$Generation) # levels = c('25','26','27','28','29','30','31','32','33','34','35'))



F25_35_selected <- subset(F25_35_selected,is.na(F25_35_selected$Early) == F)

# Figure 6

F252635 <- F25_35_selected[is.na(F25_35_selected$Early) == F & F25_35_selected$Generation %in% c(25,26,35),]

stat.test1 <- compare_means(Early~Generation, data = F252635, group.by = "Line2")

stat.test1$p[c(1,3,4,6,7,9)]
stat.test1$p.adj2 <- NA
stat.test1$p.adj2[c(4,1,6,3,7,9)] <- sort(stat.test1$p[c(1,3,4,6,7,9)])*c(6:1)

stat.test1$sigadj <- ifelse(stat.test1$p.adj2 > 0.05, "NS", 
                            ifelse(stat.test1$p.adj2 > 0.01, "*",
                                   ifelse(stat.test1$p.adj2 > 0.001, "**",
                                          ifelse(stat.test1$p.adj2 > 1e-4, "***",
                                                 ifelse(stat.test1$p.adj2 > 1e-5, "****",
                                                        ifelse(stat.test1$p.adj2 > 1e-6, "*****","******"))))))


Mate35p <- ggplot(data=F25_35_selected[is.na(F25_35_selected$Early) == F & F25_35_selected$Generation %in% c(25,26,35),], aes(x = Generation, y = Early)) + geom_violin() + facet_grid(rows=vars(Line2)) + 
  stat_summary(aes(y = Early,group=1), fun.y=mean, colour="red", geom="line",group=1) + geom_boxplot(width=0.1) + scale_color_manual(values=c("black","darkgrey")) + labs(y = "Relative mating time, hours") + ylim(-4.5,2.5) +
  stat_pvalue_manual(stat.test1, label = "sigadj", y.position = c(1.25,2.5,1.7)) + theme(legend.position = "none")
#stat_compare_means(aes(label=..p..),comparisons= list(c('25','26'), c('26','35'),c('25','35')),method="wilcox.test")
# geom_signif(comparisons= list(c('25','26'), c('26','35'),c('25','35')))

wilcox.test(F25_35_selected[F25_35_selected$Generation == 25,]$Early,F25_35_selected[F25_35_selected$Generation == 26,]$Early)
wilcox.test(F25_35_selected[F25_35_selected$Generation == 26,]$Early,F25_35_selected[F25_35_selected$Generation == 34,]$Early)
wilcox.test(F25_35_selected[F25_35_selected$Generation == 34,]$Early,F25_35_selected[F25_35_selected$Generation == 35,]$Early)

wilcox.test(F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 25,]$Early,F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 26,]$Early, alternative = "greater")
wilcox.test(F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 26,]$Early,F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 34,]$Early)
wilcox.test(F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 34,]$Early,F25_35_selected[F25_35_selected$Line == "1_1" & F25_35_selected$Generation == 35,]$Early)

wilcox.test(F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 25,]$Early,F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 26,]$Early, alternative = "greater")
wilcox.test(F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 26,]$Early,F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 34,]$Early)
wilcox.test(F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 34,]$Early,F25_35_selected[F25_35_selected$Line == "1_5" & F25_35_selected$Generation == 35,]$Early)

wilcox.test(F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 25,]$Early,F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 26,]$Early, alternative = "greater")
wilcox.test(F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 26,]$Early,F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 34,]$Early)
wilcox.test(F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 34,]$Early,F25_35_selected[F25_35_selected$Line == "2" & F25_35_selected$Generation == 35,]$Early)

#write.csv(as.data.frame(compare_means(Early~Generation, data =F25_35_selected[is.na(F25_35_selected$Early) == F & F25_35_selected$Generation != 28,], group.by="Line")), "C:/Users/yea017/Dropbox/CSIRO/Publication/Mating time choice/Table2.csv",row.names = F)

as.data.frame(compare_means(Early~Generation, data =F25_35_selected[is.na(F25_35_selected$Early) == F & F25_35_selected$Generation %in% c(25,26,35),], group.by="Line"))

F25_calprop$Callus <- as.numeric(str_count(F25_calprop$Callus3,"Y"))
Iso26_callus <- subset(General26,General26$Line == "Day" & !(General26$Replicate %in% c(1,2,3)))

Mean26isocal <- aggregate(Iso26_callus$Frequency, by =list(Iso26_callus$Replicate), FUN=sum)
Mean26isocal$sumCal <- aggregate(Iso26_callus$Frequency*(Iso26_callus$Callus-1), by =list(Iso26_callus$Replicate), FUN=sum)$x
Mean26isocal$Average <- Mean26isocal$sumCal/Mean26isocal$x

F26_isocal <- c(rep(subset(General26, General26$Replicate == '1_1')$Callus-1, subset(General26, General26$Replicate == '1_1')$Frequency),
                rep(subset(General26, General26$Replicate == '1_5')$Callus-1, subset(General26, General26$Replicate == '1_5')$Frequency), 
                rep(subset(General26, General26$Replicate %in% c('2_1','2_2','2_3','2_4','2_5','2_6'))$Callus-1, subset(General26, General26$Replicate %in% c('2_1','2_2','2_3','2_4','2_5','2_6') )$Frequency))


## Parent offspring estimate
Fem25$Callus2 <- as.numeric(ifelse(Fem25$Callus3 == "try", 4, as.character(Fem25$Callus3)))
Mean26isocal$F25 <- Fem25$Callus2[c(1,3,5,6,7,8,9,10,11,12,13,14,16,17,18)]
summary(lm(Average~F25, data=Mean26isocal[1:6,]))
summary(lm(Average~F25, data=Mean26isocal[7:12,]))
summary(lm(Average~F25, data=Mean26isocal[13:15,]))


Final35$Callus4 <- as.numeric(str_count(Final35$Callus3,"Y"))
Final35only <- subset(Final35, Final35$Line %in% c("R1_1","R1_5","R2"))
Final35only$Line <- substring(as.character(Final35only$Line),2)
Final35only$Line <- factor(Final35only$Line, levels = c("1_1","1_5","2"))
levels(Final35only$Line) <- c("1_1","1_5","2")


F25_26_35cal <- data.frame(Callus = c(rep(rep(F25_calprop$Callus,F25_calprop$N),3), F26_isocal, Final35only$Callus4),
                           Generation = c(rep(25,3*(9+2+7+49+71)),rep(26,length(F26_isocal)), rep(35,length(Final35only$Callus4))),
                           Replicate = c(rep('1_1',9+2+7+49+71),rep('1_5',9+2+7+49+71),rep('2',9+2+7+49+71), rep('1_1', length(rep(subset(General26, General26$Replicate == '1_1')$Callus-1, subset(General26, General26$Replicate == '1_1')$Frequency))),
                                         rep('1_5',length(rep(subset(General26, General26$Replicate == '1_5')$Callus-1, subset(General26, General26$Replicate == '1_5')$Frequency))),
                                         rep('2',length(rep(subset(General26, General26$Replicate %in% c('2_1','2_2','2_3','2_4','2_5','2_6'))$Callus-1, subset(General26, General26$Replicate %in% c('2_1','2_2','2_3','2_4','2_5','2_6') )$Frequency))),
                                         as.character(Final35only$Line)))
F25_26_35cal$NewCallus <- ifelse(F25_26_35cal$Callus %in% c(3,4), 3, F25_26_35cal$Callus)

F25_26_35cal$Selected <- ifelse(F25_26_35cal$Generation > 26, "No selected","Early selected")
stat.test2 <- compare_means(NewCallus~Generation, data = F25_26_35cal, group.by = "Replicate")

stat.test2$p[c(1,3,4,6,7,9)]
stat.test2$p.adj2 <- NA
stat.test2$p.adj2[c(4,1,6,3,7,9)] <- sort(stat.test2$p[c(1,3,4,6,7,9)])*c(6:1)


stat.test2$sigadj2 <- ifelse(stat.test2$p.adj2 > 0.05, "NS", 
                             ifelse(stat.test2$p.adj2 > 0.01, "*",
                                    ifelse(stat.test2$p.adj2 > 0.001, "**",
                                           ifelse(stat.test2$p.adj2 > 1e-4, "***",
                                                  ifelse(stat.test2$p.adj2 > 1e-5, "****",
                                                         ifelse(stat.test2$p.adj2 > 1e-6, "*****","******"))))))

abbrev_y <- c('bbbb','bbbY','bbYY','bYYY &\nYYYY')

AggF252634C <- aggregate(F25_26_35cal$NewCallus, by=list(F25_26_35cal$NewCallus,F25_26_35cal$Generation,F25_26_35cal$Replicate,F25_26_35cal$Selected), FUN=length)
str(AggF252634C)
names(AggF252634C) <- c("Callus","Generation","Replicate","Selected","N")


Cal35p <- ggplot(data=F25_26_35cal, aes(x = as.factor(Generation), y = Callus)) +geom_point(aes(size=N), data = AggF252634C) + scale_size(range = c(1,10))  + facet_grid(rows=vars(Replicate)) + 
  stat_summary(aes(y = Callus,group=1), fun.y=mean, colour="red", geom="line",group=1)  + scale_color_manual(values=c("black","darkgrey")) + labs(x = "Generation", y = "Callus class") +
  theme(legend.position = "none") + stat_pvalue_manual(stat.test2, label = "sigadj2", y.position = c(3.65,4.5,4)) +
  scale_y_continuous(breaks = 0:3, labels = abbrev_y, limits = c(0,4.2))# + geom_dotplot(aes(color=Selected),binaxis='y',stackdir='center', stackratio=1, dotsize=0.5, binwidth = 1/50) + geom_boxplot(aes(color=Selected),width=0.1) 
#+ geom_violin(aes(color=Selected)) + geom_signif(comparisons= list(c('25','26'), c('26','35')))


jpeg(".../Figure6.jpg",width = 3800, height = 4500, res = 600)
ggdraw() +
  draw_plot(Mate35p, x = 0.01, y = 0, width = 0.49, height = 1) +
  # draw_plot(F25_cal, x = 0.05, y = 0.5, width = 0.95, height = 0.25) +
  draw_plot(Cal35p, x = 0.51, y = 0, width = 0.49, height = 1) +
  draw_plot_label(label = c('A','B'), size = 15,
                  x = c(0.005,0.505), y = c(1,1), fontface = "bold")

dev.off()


### Heterogeneity in lines after selection (F26)
kruskal.test(Male26Early2$EarlyR1,Male26Early2$Line) ## Use this