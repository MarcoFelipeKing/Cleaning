
#RUN THIS PIECE FIRST BEFORE RUNNING FUNCTION CODE

# using Marco's Markov chains to simulate different contact behaviors
if("readr" %in% rownames(installed.packages())==FALSE){install.packages("readr"); require(readr)}else{require(readr)}
if("plyr" %in% rownames(installed.packages())==FALSE){install.packages("plyr"); require(plyr)}else{require(plyr)}
if("markovchain" %in% rownames(installed.packages())==FALSE){install.packages("markovchain"); require(markovchain)}else{require(markovchain)}

#library(readr)
movsdf.rbind <- read.csv("movsdf.rbind_orientationcorrected.csv")

for(i in 1:length(table(movsdf.rbind$ActivityID))){
  if (i==1){
    movsdf.rbind.new<-movsdf.rbind[movsdf.rbind$ActivityID==i & movsdf.rbind$Surface!="AlcOutside" &
                                     movsdf.rbind$Surface!="ApronOn" & movsdf.rbind$Surface!="ApronOff" &
                                     movsdf.rbind$Surface!="GlovesOn" & movsdf.rbind$Surface!="GlovesOff"&
                                     movsdf.rbind$Surface!="Alc" ,]
  }else{
    movsdf.rbindtemp<-movsdf.rbind[movsdf.rbind$ActivityID==i & movsdf.rbind$Surface!="AlcOutside" &
                                     movsdf.rbind$Surface!="ApronOn" & movsdf.rbind$Surface!="ApronOff" &
                                     movsdf.rbind$Surface!="GlovesOn" & movsdf.rbind$Surface!="GlovesOff"&
                                     movsdf.rbind$Surface!="Alc",]
    movsdf.rbind.new<-rbind(movsdf.rbind.new,movsdf.rbindtemp)
  }
}


#####
# 2.3 Aggregating surfaces into categories for Transition Matrices
#detach("package:dplyr", unload = TRUE)
#library(plyr)
movsdf.rbind.new$SurfaceCategories<-revalue(movsdf.rbind.new$Surface,c("In"="In",
                                                                       "Door"="FarPatient",
                                                                       "Other"="FarPatient",
                                                                       "Table"="NearPatient",
                                                                       "Tray"="Equipment" ,
                                                                       "Patient"="Patient",
                                                                       "Sharps"="Equipment"  ,
                                                                       "Waste"="HygieneInside" ,
                                                                       "Sink"="HygieneInside",
                                                                       "Soap"="HygieneInside" ,
                                                                       "PaperTowel"="HygieneInside" ,
                                                                       "IV"="Equipment" ,
                                                                       "Out"="Out" ,
                                                                       #"Alc"="Alcohol" ,
                                                                       "ObsTrolley"="Equipment",
                                                                       "Wipes"="HygieneInside",
                                                                       "Bed"="NearPatient",
                                                                       "Chair"="NearPatient",
                                                                       "Stethoscope"="Equipment"))

# Creating lists so we don't count transitions between "out" and first state of next observed episode
# It looksl ike Activity ID is unique per observation

TObs.list<-list()
TIV.list<-list()
TRounds.list<-list()

for (i in 1:max(movsdf.rbind.new$ActivityID)){
  TObs.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Obs" & movsdf.rbind.new$ActivityID==i]
  
  TIV.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="IV" & movsdf.rbind.new$ActivityID==i]
   
  TRounds.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Rounds" & movsdf.rbind.new$ActivityID==i]
  
}

require(markovchain)
TObs<-markovchainFit(TObs.list)

TIV<-markovchainFit(TIV.list)

TRounds<-markovchainFit(TRounds.list)
