# Marco-Felipe King, Amanda Wilson

# SET PARENT DIRECTORY & LOAD/INSTALL PACKAGES -----

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}
if("triangle" %in% rownames(installed.packages())==FALSE){install.packages("triangle"); require(triangle)}else{require(triangle)}
if("tidyverse" %in% rownames(installed.packages())==FALSE){install.packages("tidyverse"); require(tidyverse)}else{require(tidyverse)}
if("ggpubr" %in% rownames(installed.packages())==FALSE){install.packages("ggpubr"); require(ggpubr)}else{require(ggpubr)}
if("ramify" %in% rownames(installed.packages())==FALSE){install.packages("ramify"); require(ramify)}else{require(ramify)}
if("data.table" %in% rownames(installed.packages())==FALSE){install.packages("data.table"); require(data.table)}else{require(data.table)}
if("reshape2" %in% rownames(installed.packages())==FALSE){install.packages("reshape2"); require(reshape2)}else{require(reshape2)}
if("magrittr" %in% rownames(installed.packages())==FALSE){install.packages("magrittr"); require(magrittr)}else{require(magrittr)}
if("magrittr" %in% rownames(installed.packages())==FALSE){install.packages("magrittr"); require(magrittr)}else{require(magrittr)}
if("deSolve" %in% rownames(installed.packages())==FALSE){install.packages("deSolve"); require(deSolve)}else{require(deSolve)}

# IMPORT DATA NEEDED FOR ALL OPERATIONS DURING AUTOMATION -----

  #load in duration data
  durations<-read.csv('ContactDuration.csv')
  #read in bootstrapped values for dose-response
  exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')
 
NSIM=100  
# -------- Cleaning ODE ---------------------------------- 
  
  #Define functon to solve ODE
  ode_model <- function(time_space, initial_contamination, parameters){
    with(
      as.list(c(initial_contamination, parameters)),{
        dContamination<- Contamination*r*(1-Contamination/C)-d*exp(-g*time_space)*Contamination
        return(list(dContamination))
      }
    )
  }
  
  #sim<-ode(initial_contamination, time_space, ode_model, parameters, method = "radau",atol = 1e-4, rtol = 1e-4)
  
  deterministic_run<-function(precision,initial_contamination,parameters){
    tmax = 24
    time_space = seq(0, tmax, length.out = precision+1)
    
    #sim=odeint(ode_model,initial_contamination,time_space,args=(r,C,d,g,l))
    #parameters=c(r,C,d,g,l)
    sim<- ode(initial_contamination, time_space, ode_model, parameters, method = "radau",atol = 1e-4, rtol = 1e-4)
    num_at_0=sim[(precision*0.01/50.0),2]
    num_at_1=sim[(precision*1/50.0),2]
    num_at_2=sim[(precision*2/50.0),2]
    num_at_4=sim[(precision*4/50.0),2]
    num_at_8=sim[(precision*8/50.0),2]
    num_at_24=sim[(precision*24/50.0),2]
    
    
    return(cbind(num_at_0,num_at_1,num_at_2,num_at_4,num_at_8,num_at_24))
  }

  initial_contamination=c(Contamination=59)
  #This is a palce-holder for when the ABC results are used
  parameters<-data.frame(r=runif(NSIM,0,5),C=runif(NSIM,8,16),d=runif(NSIM,0,5),g=runif(NSIM,0,4)) #This is a palce-holder for when the ABC results are used
  precision=5000 #For the ODE solver
  CleaningType<-c("Alcohol","Control","Detergent","DistilledWater")
  
  #This saves the results from ODE solver over time and for 4 cleaning types.
  surfDecay=list()
  surfDecay <- vector(mode = "list", length=NROW(CleaningType))
  temp=matrix(data=NA,nrow=NSIM,ncol=6)

  for(j in 1:NROW(CleaningType)){
    for(i in 1:NSIM){
      temp[i,]=deterministic_run(precision,initial_contamination,parameters[i,]) 
    }
    surfDecay[[j]]=temp #   Gives the total CFU count on the experimental surface 
  }
names(surfDecay)<-CleaningType



#Currently assume room is homogeneously contaminated. No gradient between near and far patient surfaces
# Set up surface concentrations based on Stephanie's % distribution: NearPatient= Bed rails, FarPatient=Table and Notes
#surfconc=list()
SurfaceTypes<-c("NearPatient","Patient","Equipment","FarPatient","Hygiene")


# --------- Behaviour -------- 

#Extracts transition matrices
suppressMessages(source("adjust_behaviors_covid.R"))

caretype<-rep(c("IV","Obs","Rounds"),each=NSIM)
behaviourList=list()
#behaviour model section
for(ii in 1:NROW(caretype)){
  if (caretype[ii]=="IV"){
    prob.mat<-TIV$estimate #IV
  }else if (caretype[ii]=="Obs"){
    prob.mat<-TObs$estimate #Obs
  }else{
    prob.mat<-TRounds$estimate #Rounds
  }
  
  sample.space<-c("Equipment","FarPatient","HygieneInside","In","NearPatient","Out","Patient")
  
  behaviour<-"In" #first behaviour is "In"
  k<-2 #setting up k for while loop
  
  #while the previous behaviour is not "Out," make new behaviours
  while (behaviour[k-1]!="Out"){
    
    #If previous step was not out, create new behaviour
    #where the prob of selecting a behaviour from sample.space is determined by the row
    #of probabilities corresponding to the previous behaviour contact (prob.mat[behaviour[k-1],])
    
    behaviour[k]<-sample(sample.space,1,replace=TRUE,prob=prob.mat[behaviour[k-1],])
    
    #now we have to handle issues where someone is taking GlovesOff having never put them on
    
    #if the current selected behaviour is GlovesOff and they've never put GlovesOn or the maximum position of a previous GlovesOff is greater than
    # the maximum position of a previous GlovesOn
    behaviourcount<-1:length(behaviour)
    
    #advance contact number by 1
    k<-k+1
  }
  behaviourList[[ii]]=behaviour
  
}
names(behaviourList)=caretype

# --------- TRANSFER EFFICIENCY -------------------------------------------------------  
# Could have different lambdas and betas depending on surface type

#initialize transfer efficiencies
#lambda = hand--> surf
#beta = surf --> hand

x=unique(unlist(behaviourList, use.names = TRUE)) #Name all of the surface types
lambda<-matrix(rtrunc(NSIM*NROW(x),spec="norm",a=0,b=1,mean=0.51,sd=0.16), ncol=NROW(x))%>%as.data.frame()%>%set_colnames(x)
lambda$Patient<-rtrunc(NSIM,spec="norm",a=0,b=1,mean=0.21,sd=0.16)
beta<-matrix(rtrunc(NSIM*NROW(x),spec="norm",a=0,b=1,mean=0.51,sd=0.16), ncol=NROW(x))%>%as.data.frame()%>%set_colnames(x)
beta$Patient<-rtrunc(NSIM,spec="norm",a=0,b=1,mean=0.21,sd=0.16)

# --------- Surface Concentrations ---------- Missing caretypes. You need NSIM*(IV,Obs and DRs)=15

surfConc<-list()
surfConcPerClean<-list()#vector(mode="list", length=5)
tempsurfconcList<-list()
x=unique(unlist(behaviourList, use.names = TRUE)) #Name all of the surface types
tempsurfconc<-vector(mode="list", length=7)# temporary holder

for(t in 1:NROW(CleaningType)){     # Number of cleaning types

  for(j in 1:NSIM){   # For each simulation
    
    for(i in 1:NCOL(surfDecay[[1]])){ # Number of time-points to check
      
      tempsurfconc=surfDecay[[t]][j,i]/4.5^2*pi*runif(7) # CFU per cm^2. Random value
      
      tempsurfconcList[[i]]=as.data.frame(t(tempsurfconc))%>%set_colnames(x)# A single row of the decay is now in a list
    }
    
    surfConcPerClean[[j]]=tempsurfconcList #A single cleaning regime is entirely in a list
  }
  
  surfConc[[t]]=surfConcPerClean #This gives all the surfaces in the room at any given time, regardless of whether they are touched.
}
names(surfConc)=c(CleaningType) #Done
rm(tempsurfconc,tempsurfconcList,temp)
  


# -------------- FRACTIONAL HAND SURFACE AREA ------------------------------------------------------------------ 

SH<-list()
#Initialize fraction of hand used for contact
for(i in 1:(NSIM*NROW(unique(caretype)))){ 
  SHtemp<-rep(NA,length(behaviourList[[i]]))
  
  #fractional surface area for open hand grip; from AuYeung
  SHtemp[behaviourList[[i]]=="In" | behaviourList[[i]]=="Out"]<-runif(length(SHtemp[behaviourList[[i]]=="In" | behaviourList[[i]]=="Out"]),0.10,0.17) #min and max of left and right hands closed hand grip in AuYeung et al. (2008)
  
  #fractional surface area for patient contact (front partial fingers to full front palm with fingers)
  SHtemp[behaviourList[[i]]=="Patient"]<-runif(length(SHtemp[behaviourList[[i]]=="Patient"]),0.04,0.25) 
  
  #fractional surface area for variety of grip types (non "in"/"out" contacts)
  SHtemp[behaviourList[[i]]=="Equipment"|behaviourList[[i]]=="FarPatient"|behaviourList[[i]]=="NearPatient"|behaviourList[[i]]=="HygieneInside"]<-runif(length(SHtemp[behaviourList[[i]]=="Equipment"|behaviourList[[i]]=="FarPatient"|behaviourList[[i]]=="NearPatient"|behaviourList[[i]]=="HygieneInside"]),0.04/5,0.25)
  #min and max of left and right hands in AuYeung et al. (2008) for various hand grip and hand press contacts (hand immersion contacts not included)
  #from single fingertip up to full palm
  SH[[i]]=SHtemp
  print(i)
}

# --------- PAM ---------- 



C<-list()
CNSIM<-list()
CTime<-list()
CMaster<-list()
#NSIM=5
for(clean in 1:NROW(CleaningType)){
  for(t in 1:6){
    for(NSIM in 1:NSIM){ 
      for(i in 1:NROW(behaviourList)){ #This includes NSIM because NSIM*CareTypes=5*3=15
        Cmatrix=matrix(nrow=1,ncol=NROW(behaviourList[[i]])) #Temporary matrix to store PAM
        for(n in 1:(NROW(behaviourList[[i]]))){
          s=1 #Counter to run through the surfaces touched, for the SH list.
          if(n==1){ 
            Cmatrix[1,n]=lambda[NSIM,1]*0.04*(surfConc[[clean]][[NSIM]][[t]]$In)
          } else {  
            
            if(behaviourList[[i]][n]=="FarPatient"){
              Cmatrix[1,n]= Cmatrix[1,n-1]-SH[[i]][s]*(lambda$FarPatient[NSIM]*Cmatrix[1,n-1]-beta$FarPatient[NSIM]*surfConc[[clean]][[NSIM]][[t]]$FarPatient)
            } else if(behaviourList[[i]][n]=="Out"){
              Cmatrix[1,n]=Cmatrix[1,n-1]-SH[[i]][s]*(lambda$Out[NSIM]*Cmatrix[1,n-1] -beta$Out[NSIM]*surfConc[[clean]][[NSIM]][[t]]$Out)
            } else if(behaviourList[[i]][n]=="NearPatient"){
              Cmatrix[1,n]=Cmatrix[1,n-1]-SH[[i]][s]*(lambda$NearPatient[NSIM]*Cmatrix[1,n-1]-beta$NearPatient[NSIM]*surfConc[[clean]][[NSIM]][[t]]$NearPatient)
            } else if(behaviourList[[i]][n]=="Equipment"){
              Cmatrix[1,n]=Cmatrix[1,n-1]-SH[[i]][s]*(lambda$Equipment[NSIM]*Cmatrix[1,n-1]-beta$Equipment[NSIM]*surfConc[[clean]][[NSIM]][[t]]$Equipment)
            } else if(behaviourList[[i]][n]=="Patient"){
              Cmatrix[1,n]=Cmatrix[1,n-1]-SH[[i]][s]*(lambda$Patient[NSIM]*Cmatrix[1,n-1]-beta$Patient[NSIM]*surfConc[[clean]][[NSIM]][[t]]$Patient)
            } else if(behaviourList[[i]][n]=="HygieneInside"){
              Cmatrix[1,n]=Cmatrix[1,n-1]-SH[[i]][s]*(lambda$HygieneInside[NSIM]*Cmatrix[1,n-1]-beta$HygieneInside[NSIM]*surfConc[[clean]][[NSIM]][[t]]$HygieneInside)
            }
            
          }
          s=s+1 #Counter for the surface fraction to avoid another for loop
        }
        C[[i]]=Cmatrix #For a single episode of care
      }
      CNSIM[[NSIM]]=C  #For all the simulated lambda values #There's something missing with this!
    }
    CTime[[t]]=CNSIM #For all the simulated times 
  } 
  CMaster[[clean]]=CTime #For all simulated cleaningtype
}
#Hasta aqui
# a=match(behaviourList[[1]],colnames(temp))
# sum(temp(a))
# test=CMaster
# names(test)<-CleaningType
# 
# names(test$Alcohol[[1]])=c(1:NSIM)

#------------- Plotting---------------------------------------------------------------------------
# D<-list()
# DNSIM<-list()
# DTime<-list()
# DMaster<-list()
# 
#   
# for(CC in 1:4){
#   for(t in 1:6){
#     for(N in 1:1){
#     for(i in 1:15){
#       
#       Dmatrix=runif(15)
#       D[[i]]=Dmatrix
#     }
#     DTime[[t]]=D
#     }
#     DNSIM[[N]]=DTime
#   }
#   DMaster[[CC]]=DTime
# }


df <- tibble(lists = CMaster) %>% 
  mutate(CleaningType = row_number()) %>% 
  unnest_longer(lists, indices_to = "TimePoint") %>% 
  unnest_longer(lists, indices_to = "Replicate") %>%
  unnest_longer(lists, indices_to = "BehaviourObservation") %>%
  mutate(lists = map(lists, as.vector)) %>% unnest_longer(lists, indices_to = "sub_sub_observation")
  
df$CleaningType<-as.factor(df$CleaningType)

df %>% 
  ggplot(aes(TimePoint, lists, group = TimePoint,fill=CleaningType)) + 
  #stat_summary(fun.data = "mean_cl_boot", size = 0.5,position = position_dodge(width=0.75))+
  geom_boxplot() +
  scale_y_continuous(name = "CFUs on hands")+
  scale_x_discrete(name="Time after cleaning (h)")+
  theme_pubclean()+
  facet_wrap(~ CleaningType)


