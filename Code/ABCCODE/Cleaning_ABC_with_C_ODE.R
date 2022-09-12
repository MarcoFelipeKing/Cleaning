# Cleaning ODE in C
# dC/dt=C*r*(1-C/Cap)-d*exp(-g*time_space)*C

pacman::p_load(dplyr,deSolve,ggpubr,hrbrthemes)
setwd("Code/ABCCODE/")

# Check the system architecture before loading the compiled function
if(Sys.info()["machine"]=="arm64"){
  system("R CMD SHLIB cleaning_ODE_arm64.c")
  dyn.load("cleaning_ODE_arm64.so")  
}else{
  system("R CMD SHLIB cleaning_ODE_X86.c")
  dyn.load("cleaning_ODE_X86.c")  
}
dynLoad <- function(dynlib){
  dynlib <- paste(dynlib, .Platform$dynlib.ext, sep = "")
  dyn.load(dynlib)
}
#ARM
dynLoad("cleaning_ODE_arm64")
#OR 
dynLoad("cleaning_ODE_X86")

# ODE Compiled C version ----

deterministic_run_C<-function(precision,initial_contamination,parameters){
  tmax = 24
  times = seq(0, tmax, length.out = precision) #precision+1
  Y <-c(Contamination=initial_contamination) 
  sim <- ode(initial_contamination, times, func = "derivs", dllname = "cleaning_ODE_arm64", initfunc = "initmod",
             parms = parameters, method = "radau", atol = 1e-4, rtol = 1e-4)
  num_at_0=initial_contamination
  num_at_1=sim[(precision*1/tmax),2]
  num_at_2=sim[(precision*2/tmax),2]
  num_at_4=sim[(precision*4/tmax),2]
  num_at_8=sim[(precision*8/tmax),2]
  num_at_24=sim[(precision*24/tmax),2]
  
  
  return(cbind(num_at_1,num_at_2,num_at_4,num_at_8,num_at_24))
}


##################################################################################################################
## Applying the ABC rejection algorithm
df <- vroom::vroom("../../LabData/cleaningExptBeth_summary.csv") %>% janitor::clean_names()
initial_contamination=c(Contamination=df$m[df$cleaning_type=="Detergent" & df$time_after_cleaning=="0"] ) #median 34
experimental_data = df$m[df$cleaning_type=="Detergent"]# c(19,5.4,5.4,2.4,8.6)#c(11.5,5,2,0,6)#
s=df$s[df$cleaning_type=="Detergent"]#c(26,2.3,4.7,4.3,4.3)
sample_size = 1000

parameter_sample<- matrix(data=NA,nrow=1,ncol=4)
total_trials=0  #Starts the counter of run trials
accepted_trials=0 

#Create empty vector for the distances
distances<-c()
delta=80#0.006136007
precision=5000 #this is delta t

while (NROW(parameter_sample) <= sample_size){
  #foreach(total_trials=1:2, .packages=c("deSolve")) %dopar% {
  
  # Define parameters ----------
  # The prior distributions we use are d ~ U(0.001,10.0), C ~ U(200,1200), r ~ U(0.001,1.0), g ~ U(0.001,1.0). 
  # We begin by sampling from these distributions and simulating the process
  trial_r = runif(1,1e-5,5)
  # trial_rU = runif(1,1E-1,1000)
  trial_C = runif(1,1E1,100)
  # trial_CU = runif(1,1e2,1e6)
  trial_d = runif(1,1e-5,10)
  trial_g = runif(1,1e-5,10)
  # trial_l = runif(1,1E-5,1e2)
  total_trials=total_trials+1.0
  #print(total_trials)
  # parameters=c(r=trial_r,C=trial_C,d=trial_d,g=trial_g) #For R code - but it doesn't matter what order they are in.
  # one_run = deterministic_run(precision,initial_contamination,parameters)# R code
  parameters <- c(C = trial_C, d = trial_d, g=trial_g,r=trial_r); #For C code
  #Treated
  one_run = deterministic_run_C(precision,initial_contamination,parameters)# For C code
  #one_runU = deterministic_run_C(precision,initial_contaminationU,c(C = trial_C, d = 0, g=0,r=trial_rU))# For C code
  # Now we find the Euclidean distance between the simulated output and the
  # experimental results. delta is the threshold that the Euclidean distance
  # must be less than for us to accept the trial parameters into our sample.
  euclidean_distance = abs(dist(rbind(one_run,experimental_data)))#+abs(dist(rbind(one_runU,experimental_dataU)))#Distance(one_run, experimental_data,s)+Distance(one_runU, experimental_dataU,sU) # abs(dist(rbind(one_run,experimental_data))) #
  #print(euclidean_distance)# print(parameter_sample,euclidean_distance)
  #print(total_trials)
  if (euclidean_distance < delta){
    #parameter_sample = parameter_sample[!is.na(parameter_sample)];
    # print(euclidean_distance)
    parameter_sample=rbind(parameter_sample, c(trial_C,trial_d,trial_g,trial_r))
    distances=rbind(distances,euclidean_distance)
    accepted_trials=accepted_trials+1.0
    print(paste0("Trial number accepted: ",accepted_trials))
  }
  #else{
  #  print(euclidean_distance)
  #  }
}
print(paste0("Percentage of trials accepted: ",(100*accepted_trials/total_trials)))#Order the distances and then the parameter sample
parameter_sample<-parameter_sample[order(distances), ]
distances <- distances[order(distances)]
parameter_sample<-parameter_sample%>%as.data.frame()#%>%set_colnames(c("r", "C", "d","g","l"))
colnames(parameter_sample)<-c( "Cmax", "d","g","r")
summary(distances)



#  Pairs plot -------------------------------------------------------------


GGally::ggpairs(parameter_sample)

#  Individual histograms plot -------------------------------------------------------------

parameter_sample %>% 
drop_na() %>% 
  pivot_longer(c(r,Cmax,d,g)) %>% 
  ggplot()+
  # geom_density(aes(x=value),fill="#008080",alpha=0.7)+
  geom_histogram(aes(x=value),fill="purple",bins = 15,alpha=0.7)+
  facet_wrap(~name,scales="free")+
  hrbrthemes::theme_ipsum() ->a

# Experimental data
tmax = 24
time_space = seq(0, tmax, length.out = precision+1)

ode_model <- function(time_space, initial_contamination, parameters){
  with(
    as.list(c(initial_contamination, parameters)),{
      dContamination <- Contamination*r*(1-Contamination/C)-d*exp(-g*time_space)*Contamination
      return(list(dContamination))
    }
  )
}
sim1<-ode(initial_contamination, time_space, ode_model, parameter_sample[1,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim2<-ode(initial_contamination, time_space, ode_model, parameter_sample[2,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim3<-ode(initial_contamination, time_space, ode_model, parameter_sample[3,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim4<-ode(initial_contamination, time_space, ode_model, parameter_sample[4,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim<-rbind(sim1[,c(1,2)],sim2[,c(1,2)],sim3[,c(1,2)],sim4[,c(1,2)])%>%as.data.frame()
sim$CurveNum<-rep(1:4,each=nrow(sim1))
#lines(sim,col="blue")
sim<-sim%>%as.data.frame()%>%set_colnames(c("Time","Contamination","CurveNum"))
df<- read.csv("../../LabData/cleaningExptBeth.csv")#data.frame(Time=c(0,1,2,4,8,24),Contamination=c(initial_contamination,experimental_data))
df<-df[,-2]%>%set_colnames(c("Contamination","CleaningType","Time","Replica"))
df$CleaningType<-factor(df$CleaningType, levels = c("Control","Alcohol","Detergent", "Distilled Water"))

ggplot() + 
  # geom_point(data=subset(df,CleaningType=="Alcohol"), aes(TimeAfterCleaning, Count),color="#008080") +
  geom_errorbar(data=df %>% 
               group_by(TimeAfterCleaning,CleaningType) %>% 
               summarise(Mean=mean(Count,na.rm=TRUE),SE=sd(Count,na.rm =TRUE)/sqrt(5)) %>% ungroup() %>% filter(CleaningType=="Detergent"),
               aes(x=TimeAfterCleaning,y=Mean,ymin=Mean-SE,ymax=Mean+SE),color="#008080",size=1.1)+
  geom_point(data=df %>% 
                  group_by(TimeAfterCleaning,CleaningType) %>% 
                  summarise(Mean=mean(Count,na.rm=TRUE),SE=sd(Count,na.rm =TRUE)/sqrt(5)) %>% ungroup() %>% filter(CleaningType=="Detergent"),
                aes(x=TimeAfterCleaning,y=Mean),color="#008080",size=1.1)+
  # stat_summary(data=subset(df,CleaningType=="Alcohol"), aes(TimeAfterCleaning, Count),color="#008080",
  #              fun.data = "mean_cl_boot", size = 0.5,position = position_dodge(width=0.75))+
  geom_line(data = sim,aes(color=as.factor(CurveNum),x=time,y=Contamination),size=0.3)+ #"#008080"
  scale_y_continuous(name  = "S. aureus colonies",limits = c(0,110))+ #trans = "log10",
  scale_colour_brewer(palette = "Set2")+
  labs(color='Trial curve')+# +scale_color_discrete(name = "Trial curve")+  
  theme_pubclean()->b

ggarrange(b,a,ncol = 2,common.legend = FALSE,labels = "AUTO")

ggsave(filename = "Outputs/prediction_alcohol.png",dpi = 600,bg = "white",width=12)

write.csv("Outputs/parameter_sample_Alcohol_1K.csv",row.names = FALSE)
