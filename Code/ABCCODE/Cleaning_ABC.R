

#author: MFK

#if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
#if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}
if("tidyverse" %in% rownames(installed.packages())==FALSE){install.packages("tidyverse"); require(tidyverse)}else{require(tidyverse)}
if("ggpubr" %in% rownames(installed.packages())==FALSE){install.packages("ggpubr"); require(ggpubr)}else{require(ggpubr)}
#if("ramify" %in% rownames(installed.packages())==FALSE){install.packages("ramify"); require(ramify)}else{require(ramify)}
#if("data.table" %in% rownames(installed.packages())==FALSE){install.packages("data.table"); require(data.table)}else{require(data.table)}
#if("reshape2" %in% rownames(installed.packages())==FALSE){install.packages("reshape2"); require(reshape2)}else{require(reshape2)}
if("deSolve" %in% rownames(installed.packages())==FALSE){install.packages("deSolve"); require(deSolve)}else{require(deSolve)}
if("magrittr" %in% rownames(installed.packages())==FALSE){install.packages("magrittr"); require(magrittr)}else{require(magrittr)}
if("mnormt" %in% rownames(installed.packages())==FALSE){install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/mnormt_1.5-7.tgz", repos=NULL, type="source"); require(mnormt)}else{require(mnormt)}
if("EasyABC" %in% rownames(installed.packages())==FALSE){install.packages("EasyABC"); require(EasyABC)}else{require(EasyABC)}
if("smfsb" %in% rownames(installed.packages())==FALSE){install.packages("smfsb"); require(smfb)}else{require(smfsb)}


Distance<-function(x,y,s){  # computes the Euclidean distance between two lists of the same length
  if (length(x) == length(y)){
    sqrt(sum(((x)-(y))/(s))^2)
              }
  else{
    print( 'lists not the same length')}
}




#parameters=c(r=runif(1),C=runif(1),d=runif(1),g=runif(1),l=runif(1))
#time_space=seq(0, tmax, length.out = precision+1)
#initial_contamination <- c(Contamination = 23)

# l= recontamination, r= , g, C
ode_model <- function(time_space, initial_contamination, parameters){
  with(
    as.list(c(initial_contamination, parameters)),{
      dContamination <- Contamination*r*(1-Contamination/C)-d*exp(-g*time_space)*Contamination
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
#num_at_0=sim[int(precision*0.1/50.0)]
num_at_1=sim[(precision*1/50.0),2]
num_at_2=sim[(precision*2/50.0),2]
num_at_4=sim[(precision*4/50.0),2]
num_at_8=sim[(precision*8/50.0),2]
num_at_24=sim[(precision*24/50.0),2]


return(cbind(num_at_1,num_at_2,num_at_4,num_at_8,num_at_24))
}

##################################################################################################################
## Applying the ABC rejection algorithm
initial_contamination=c(Contamination=59) #median 34
experimental_data = c(19,5.4,5.4,2.4,8.6)#c(11.5,5,2,0,6)#
s=c(26,2.3,4.7,4.3,4.3)
sample_size = 5
#parameter_sample <- c()
parameter_sample<- matrix(data=NA,nrow=1,ncol=5)
total_trials=0
accepted_trials=0

# Files with posteriors
# Posterior_r = open("Posterior_r.txt","w")
# Posterior_C = open("Posterior_C.txt","w")
# Posterior_d = open("Posterior_d.txt","w")
# Posterior_g = open("Posterior_g.txt","w")
# Posterior_l = open("Posterior_l.txt","w")
distances<-c()

precision=5000 #this is delta t
# library(foreach)
# library(doParallel)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
while (NROW(parameter_sample) < sample_size){
#foreach(total_trials=1:2, .packages=c("deSolve")) %dopar% {
  
  # The prior distributions we use are d ~ U(0.001,10.0), C ~ U(200,1200), r ~ U(0.001,1.0), g ~ U(0.001,1.0). 
  # We begin by sampling from these distributions and simulating the process
trial_r = runif(1,1E-1,4)
trial_C = runif(1,6,15)
trial_d = runif(1,1,4)
trial_g = runif(1,1E-2,1)
trial_l = runif(1,1E-5,2)
total_trials=total_trials+1.0
#print(total_trials)
parameters=c(r=trial_r,C=trial_C,d=trial_d,g=trial_g,l=trial_l)


one_run = deterministic_run(precision,initial_contamination,parameters)#

#experimental_data_noP2 = [652.0, 556.0, 424.5, 467.5, 428.0, 550.0, 672.0]

# Now we find the Euclidean distance between the simulated output and the
# experimental results. delta is the threshold that the Euclidean distance
# must be less than for us to accept the trial parameters into our sample.
delta = 10 #dictates distance
euclidean_distance = abs(dist(rbind(one_run,experimental_data))) #Distance(one_run, experimental_data,s)#_noP2)
#print(len(parameter_sample),euclidean_distance)
#print(total_trials)

if (euclidean_distance < delta){
  #parameter_sample = parameter_sample[!is.na(parameter_sample)];
  parameter_sample=rbind(parameter_sample, c(trial_r,trial_C,trial_d,trial_g,trial_l))
  distances=rbind(distances,euclidean_distance)
  accepted_trials=accepted_trials+1.0
  print(paste0("Trial number accepted: ",accepted_trials))
}
  #else{
  #  print(euclidean_distance)
  #  }

}
print(paste0("Percentage of trials accepted: ",(100*accepted_trials/total_trials)))
#Order the distances and then the parameter sample
parameter_sample<-parameter_sample[order(distances), ]
parameter_sample<-parameter_sample%>%as.data.frame()%>%set_colnames(c("r", "C", "d","g","l"))

#stopImplicitCluster()

# Posterior_r.write(str(trial_r))
# Posterior_r.write("\n")
# Posterior_C.write(str(trial_C))
# Posterior_C.write("\n")
# Posterior_d.write(str(trial_d))
# Posterior_d.write("\n")
# Posterior_g.write(str(trial_g))
# Posterior_g.write("\n")
# Posterior_l.write(str(trial_l))
# Posterior_l.write("\n")
#print(parameter_sample) #This prints all the values accepted


#write.csv(parameter_sample,"parameter_sample_Alcohol.csv",row.names = FALSE)

## Plotting a single best curve


tmax = 24
time_space = seq(0, tmax, length.out = precision+1)

sim1<-ode(initial_contamination, time_space, ode_model, parameter_sample[1,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim2<-ode(initial_contamination, time_space, ode_model, parameter_sample[2,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim3<-ode(initial_contamination, time_space, ode_model, parameter_sample[3,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim4<-ode(initial_contamination, time_space, ode_model, parameter_sample[4,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim<-rbind(sim1[,c(1,2)],sim2[,c(1,2)],sim3[,c(1,2)],sim4[,c(1,2)])%>%as.data.frame()
sim$CurveNum<-rep(1:4,each=nrow(sim1))
#lines(sim,col="blue")
sim<-sim%>%as.data.frame()%>%set_colnames(c("Time","Contamination","CurveNum"))
df<- read_csv("~/Downloads/Surface Cleaning/Data/cleaningExptBeth.csv")#data.frame(Time=c(0,1,2,4,8,24),Contamination=c(initial_contamination,experimental_data))
df<-df[,-2]%>%set_colnames(c("Contamination","CleaningType","Time","Replica"))
df$CleaningType<-factor(df$CleaningType, levels = c("Control","Alcohol","Detergent", "Distilled Water"))

ggplot(subset(df,CleaningType=="Alcohol"), aes(Time, Contamination,)) + 
  geom_point(color="#008080") +
  #stat_summary(fun.data = "mean_cl_boot", size = 0.5,position = position_dodge(width=0.75))+
  geom_line(data = sim,aes(color=as.factor(CurveNum)),size=1.1)+ #"#008080"
  scale_y_continuous(name  = "S. aureus colonies")+ #trans = "log10",
  labs(color='Trial curve')+# +scale_color_discrete(name = "Trial curve")+  
  theme_pubclean()


#Just the experiemntal curves
ggplot(subset(df), aes(x=Time, y=Contamination/(4.5^2*pi), colour=CleaningType)) +
  geom_point(alpha=0.2,color="black")+
  stat_summary(fun.data = "mean_cl_boot", size = 0.5,position = position_dodge(width=0.75))+
  scale_y_continuous(name  = "CFU/cm^2")+ #trans = "log10",
  scale_x_continuous(name  = "Time(h) after cleaning")+ #trans = "log10",
  scale_color_brewer(palette = "Set1",name="Care Type")+
  labs(color='Cleaning Method')+# +scale_color_discrete(name = "Trial curve")+
  facet_wrap(~CleaningType,nrow = 3)+
  theme_pubr()
  

 #ggsave(filename = 'Alcohol.pdf')
 #write.csv(parameter_sample,"parameter_sample_Alcohol_500.csv",row.names = FALSE)
 # This will save a 400x400 file at 100 ppi
 #ggsave("plot.png", width=4, height=4, dpi=100)
# 
# ggplot(melt(parameter_sample),aes(x=value,))+
#   geom_histogram(color="black", fill="white")+
#   facet_wrap(~variable,ncol=2,scales = "free")+
#   theme_pubclean()

#Matrix Plot
#ppi <- 600
#png("hist_Alcohol.png", width=6*ppi, height=6*ppi, res=ppi)
#hist(parameter_sample[,-5], pch=20 , cex=1.5 , col="#69b3a220")
df_tidy<-gather(parameter_sample[,-5], Parameter, value)
ggplot(df_tidy,aes(x=value))+
  geom_histogram(fill="#69b3a2",bins = 15)+
  theme_pubclean()+
  facet_wrap(~Parameter,nrow = 2,scales = "free")
#dev.off()
png("matrix_Alcohol.png", width=6*ppi, height=6*ppi, res=ppi)
plot(parameter_sample[,-5] , pch=20 , cex=1.5 , col="#69b3a220")
dev.off()


###################################
# ABC SMC - Setting it up - Doesn't work


priors<-function(){
  # we define here the priors. In this case, 3 uniform distributions
  trial_r = runif(1,1E-1,4)
  trial_C = runif(1,6,15)
  trial_d = runif(1,1,4)
  trial_g = runif(1,1E-2,1)
  trial_l = runif(1,1E-5,2)
return(parameters=c(r=trial_r,C=trial_C,d=trial_d,g=trial_g,l=trial_l))
}

perturbation_kernel_r<-function(r,posterior,fraction){
  perturbed_r=runif(1,r-fraction*range(posterior),r+fraction*range(posterior))
return(perturbed_r)
}

perturbation_kernel_C<-function(C,posterior,fraction){
  perturbed_C=runif(1,C-fraction*range(posterior),C+fraction*range(posterior))
  return(perturbed_C)
}

perturbation_kernel_d<-function(d,posterior,fraction){
  perturbed_d=runif(1,d-fraction*range(posterior),d+fraction*range(posterior))
  return(perturbed_d)
}

perturbation_kernel_g<-function(g,posterior,fraction){
  perturbed_g=runif(1,g-fraction*range(posterior),g+fraction*range(posterior))
  return(perturbed_g)
}

perturbation_kernel_l<-function(l,posterior,fraction){
  perturbed_l=runif(1,l-fraction*range(posterior),l+fraction*range(posterior))
  return(perturbed_l)
}



prior_density<-function(r,C,d,g,l){
  #density_r=0
if (r>=1E-1 & r<=4){ ###### need to catch NANs
  density_r=1.0/(4-1E-1)
  }
  #density_C=0
if (C>=6 & C<=15){
  density_C=1.0/(15-6)
  }
  #density_d=0
if (d>=1 & d<=4){
  density_d=1.0/(4-1)
  }
if (g>=1E-2 & g<=1){
    density_g=1.0/(1-1E-2)
  }
if (l>=1E-5 & l<=2){
    density_l=1.0/(2-1E-5)
}
  else{density_r=0
  density_C=0
  density_d=0
  density_g=0
  density_l=0}
return(density_r*density_C*density_d*density_g*density_l)
}

perturbation_kernel_density<-function(r1,C1,d1,g1,l1,r2,C2,d2,g2,l2,posterior_r,posterior_C,posterior_d,posterior_g,posterior_l,fraction){
  density_r=0
if (abs(r1-r2)<fraction*range(posterior_r)){
  density_r=1.0/(2.0*fraction*range(posterior_r))
  }
density_C=0
if (abs(C1-C2)<fraction*range(posterior_C)){
  density_C=1.0/(2.0*fraction*range(posterior_C))
}
density_d=0
if (abs(d1-d2)<fraction*range(posterior_d)){
  density_d=1.0/(2.0*fraction*range(posterior_d))
}
density_g=0
if (abs(g1-g2)<fraction*range(posterior_g)){
  density_g=1.0/(2.0*fraction*range(posterior_g))
}
density_l=0
if (abs(l1-l2)<fraction*range(posterior_l)){
  density_l=1.0/(2.0*fraction*range(posterior_l))
}
return(density_r*density_C*density_d*density_g*density_l)
}
#################################################################################################################
## Applying the ABC-SMC algorithm
initial_contamination=c(Contamination=59) #median 34
experimental_data = c(19,5.4,5.4,2.4,8.6)#c(11.5,5,2,0,6)#
s=c(26,2.3,4.7,4.3,4.3)


initial_distance_threshold=100

T=10 #The number of times you will run the ABC rejection algorithm
number_of_particles=10 #These are the number of parameter estimates you make
# We will store here the posteriors for each iteration t=0,1,...,T-1. Each posterior is a vector with number_of_particles values of the parameters
# posteriors_ss = [[] for i in range(number_of_particles)]
# posteriors_sf = [[] for i in range(number_of_particles)]
# posterior_lam = [[] for i in range(number_of_particles)]

posterior_r=c()
posterior_C=c()
posterior_d=c()
posterior_g=c()
posterior_l=c()

# For each (joint) posterior parameter values; that is, each particle, we store its weight and the corresponding distance
weights = matrix(data=1,nrow=T,ncol=number_of_particles) #c()#[[] for i in range(number_of_particles)]
distances =matrix(data=NA,nrow=T,ncol=number_of_particles)#c()# [[] for i in range(number_of_particles)]

# We input the data
# cf0_data = np.asarray([8.25E+05,1.25E+05,7.30E+05,5.20E+04,2.34E+05,8.25E+04,7.70E+05,3.05E+05,1.83E+02,2.70E+02,3.70E+07,5.20E+05,2.52E+07,2.02E+08,3.35E+05,4.80E+05])
# cs0_data = np.asarray([2.40E+05,5.35E+06,2.72E+04,3.10E+06,3.50E+06,2.50E+06,5.30E+07,1.41E+08,3.80E+06,2.30E+06,3.35E+06,7.35E+02,1.75E+06,1.32E+06,4.50E+02,5.00E+00])
# 
# cf_data=np.asarray([6.67E+05,3.32E+05,6.51E+05,4.93E+04,2.39E+05,5.92E+04,9.28E+05,1.67E+06,4.92E+04,1.12E+04,1.59E+07,5.75E+05,2.67E+07,2.07E+08,4.45E+05,6.40E+05])
# cs_data=np.asarray([2.12E+05,4.30E+06,2.71E+04,3.05E+06,9.16E+05,2.90E+06,5.35E+07,1.71E+08,2.15E+06,1.46E+06,1.98E+06,4.52E+03,2.14E+06,2.54E+06,3.81E+03,1.15E+03])

#print_range=[100,200,300,400,500,600,700,800,900]
# We initialize the first iteration of the ABC-SMC
t=1
while (t<T){
  print(t)

# We consider the first particle
i=1
while (i<number_of_particles){
  if (t==1){
    # If it is the first iteration, the particle is sampled from the prior distributions
    parameters=priors()
    r=parameters[1]
    C=parameters[2]
    d=parameters[3]
    g=parameters[4]
    l=parameters[5]
  } else {
    # If it is not the first iteration, we choose the particle from the last posterior
    # I choose the particle from the previous posterior according to the vector weights (of probabilities)
    #Just choose 1 value from 1-number particles with a weighting
    k=sample(x=number_of_particles,size=1,replace=TRUE,prob=weights[t-1,])#weight[t-1,]))
    #int(np.random.choice(np.arange(number_of_particles),size=1,replace=True,p=weights[t-1]))
  
  r=posterior_r[t-1,k]
  C=posterior_C[t-1,k]
  d=posterior_d[t-1,k]
  g=posterior_g[t-1,k]
  l=posterior_l[t-1,k]
  
  # We perturb the values obtained. We perturb according to perturbation kernel K, which is defined as a Uniform in a range defined in terms of 
  # a fraction f=0.1 of the previous marginal posterior for this parameter
  r=perturbation_kernel_r(r,posterior_r[t-1],0.1)
  C=perturbation_kernel_C(C,posterior_C[t-1],0.1)
  d=perturbation_kernel_d(d,posterior_d[t-1],0.1)
  g=perturbation_kernel_g(g,posterior_g[t-1],0.1)
  l=perturbation_kernel_g(l,posterior_l[t-1],0.1)
  
  }
  # We check if the proposed parameter values have positive prior probability
  if (prior_density(r,C,d,g,l)==0)
    next
  
  # We simulate the concetrations after the contact, from our model and the chosen parameter values. cf0 and cs0 should be two vectors of values, so that cf and cs are also vectors of predicted concentrations
  # cf=cf0_data/sf-lam*(cf0_data/sf-cs0_data/ss)
  # cs=cs0_data/ss-lam*(cs0_data/ss-cf0_data/sf)
  parameters=c(r,C,d,g,l)
  
  
  one_run = deterministic_run(5000,initial_contamination,parameters)#
  
  
  # We compute the Euclidean distance between predicted and observed values
  
  #delta = 5 #dictates distance
  euclidean_distance = abs(dist(rbind(one_run,experimental_data)))
  
  # distance=0
  # for (k in range(len(cf))){
  #   distance+=(math.log10(cf[k]*sf)-math.log10(cf_data[k]))**2+(math.log10(cs[k]*ss)-math.log10(cs_data[k]))**2
  # distance=math.sqrt(distance)
  
  if(t==1){
    # In the first iteration, we use the initial distance threshold
    if (euclidean_distance<initial_distance_threshold){
      # If the distance for the proposed parameter values is smaller than the threshold, we keep the parameter values in the new posterior vector
      posterior_r<-rbind(posterior_r, r)
      posterior_C<-rbind(posterior_C, C)
      posterior_d<-rbind(posterior_d, d)
      posterior_g<-rbind(posterior_g, g)
      posterior_l<-rbind(posterior_l, l)
      
      # And we store the weights (1.0 for each particle) for this iteration, and the distance obtained
      weights<-rbind(weights,1)
      distances=rbind(distances,euclidean_distance)
    }
    
    # If not, we just try again with a new particle
    else {next} 
       
    
} else {
    # If it is not the first iteration, we use a new distance_threshold, which has been obtained as the median for the previous vector of distances
    if (distance<distance_threshold){
      # And we do the same here
      posterior_r<-rbind(posterior_r,r)
      posterior_C<-rbind(posterior_C,C)
      posterior_d<-rbind(posterior_d,d)
      posterior_g<-rbind(posterior_g,g)
      posterior_l<-rbind(posterior_l,l)
      
      # The weight has a more complicated expression now
      numerator=prior_density(r,C,d,g,l)
      denominator=0
      for (j in number_of_particles){
        denominator=denominator+weights[t-1,j]*perturbation_kernel_density(r,C,d,g,l,posterior_r[t-1,j],posterior_C[t-1,j],posterior_d[t-1,j],posterior_g[t-1,j],posterior_l[t-1,j],posterior_r[t-1],posterior_C[t-1],posterior_d[t-1],posterior_g[t-1],posterior_l[t-1],0.1)
        
        weights<-rbind(weights,numerator/denominator)
        distances=rbind(distances,euclidean_distance)
      }
      
      } else {
        next}
}
    i=i+1
    
    distance_threshold=median(distances)
    suma_weights=sum(weights)	
    
    weights=weights/suma_weights
    # for (j in (number_of_particles)){
    #   weights[t,j]=weights[t,j]/suma_weights
    #   weights/suma_weights
    # }
  
    t=t+1     
  } #Close first while statement

t=1 
}

