pacman::p_load(dplyr,deSolve,ggplot2,ggpubr,tidyr)
# Read in cleaningExptBeth.csv from LabData

df<-read.csv("LabData/cleaningExptBeth.csv")


#Plot the data with ggplot2


df$TimeAfterCleaning <- as.numeric(df$TimeAfterCleaning)
df %>% 
  mutate(MeanCountPerMl=Count*20)->df
# subset(df,CleaningType=="Control" )

df%>%
    filter(Count > 0) %>%
    ggplot(aes(x=TimeAfterCleaning,y=Count,color=CleaningType))+
    geom_point(alpha=0.3)+
    #Add errorbars to the points using stat_summary
    stat_summary(fun = mean, geom = "point",size=2,colour="black") + 
    stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=0.2),colour="black")+
    scale_color_brewer(palette="Set1")+
    scale_y_continuous(trans="log10")+
    facet_grid(CleaningType~.,scales="free")+
    hrbrthemes::theme_ipsum()+
    labs(x="Time after cleaning (hours)",y="CFU",title="Cleaning Experiment - Beth")+
    theme(legend.position="",legend.title=element_text(face="bold"))->a

ggsave(a,filename="outputs/images/cleaningExptBeth.png",width=8,height=8,dpi=300,units="in",device="png",bg="white")
a

# Plots all values together
df%>%
    filter(Count > 0) %>%
    ggplot(aes(x=TimeAfterCleaning,y=Count))+
    geom_point(alpha=0.3)+
    #Add errorbars to the points using stat_summary
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_se, geom = "errorbar")+
    scale_color_brewer(palette="Set1")+
    scale_y_continuous(trans="log10")+
    hrbrthemes::theme_ipsum()+
    labs(x="Time after cleaning (hours)",y="CFU",title="Cleaning Experiment")+
    theme(legend.position="",legend.title=element_text(face="bold"))

####### Waseem's data
dW <- read.table("LabData/cleaningExptWaseem.tsv",header = TRUE)
dW%>%
  filter(Time<24) %>% 
  ggplot(aes(x=Time,y=cfu.9cm2,color=Technique))+
  geom_point(alpha=0.3)+
  #Add errorbars to the points using stat_summary
  stat_summary(fun = mean, geom = "point",size=2,colour="black") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=0.2),colour="black")+
  facet_grid(Technique~.,scales="free")+
  scale_color_brewer(palette="Set1")+
  scale_y_continuous(trans="log10")+
  hrbrthemes::theme_ipsum()+
  labs(x="Time after cleaning (hours)",y="CFU",title="Cleaning Experiment Waseem")+
  theme(legend.position="",legend.title=element_text(face="bold"))->b
ggsave(b,filename="outputs/images/cleaningExptWaseem.png",width=8,height=8,dpi=300,units="in",device="png",bg="white")
b

# Total orginal amount on coupons is about 4.9 x 10^7
dW %>% filter(Time==0 & Technique=="C") %>% mutate(x=cfu_by_swabb_eff/0.0488) %>% summarise(mean(x))



# Plot posteriors from BethExpt_ABC.py

parameter_sample<-vroom::vroom("outputs/Posterior_Beth_Expt_with_dw.txt",col_names = c("r","C","m_de","g_de","m_di","g_di","m_dw","g_dw","l"))


ode_model<- function(t, y, parameters){
  with(
    as.list(c(y, parameters)),{
      dContamination <- Contamination*r -d * exp(-g * t) * Contamination
      return(list(dContamination))
    }
  )
}

parameter_sample %>% 
  drop_na() %>% 
  summarise(across(.cols=everything(),quantile)) ->p_summary

sim1<-ode(initial_contamination, times, ode_model, p_summary[2,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim2<-ode(initial_contamination, times, ode_model, p_summary[3,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim3<-ode(initial_contamination, times, ode_model, p_summary[4,], method = "radau",atol = 1e-4, rtol = 1e-4)
# sim4<-ode(initial_contamination, times, ode_model, p_summary[1,], method = "radau",atol = 1e-4, rtol = 1e-4)
sim<-rbind(sim1[,c(1,2)],sim2[,c(1,2)],sim3[,c(1,2)])%>%as.data.frame()
sim$CurveNum<-rep(1:3,each=nrow(sim1))
#lines(sim,col="blue")
sim<-sim%>%as.data.frame()#%>%tidyr::set_colnames(c("Time","Contamination","CurveNum"))
colnames(sim)<-c("Time","Contamination","CurveNum")

ribbon_data <- data.frame(ymin=sim %>% as_tibble()%>% filter(CurveNum==1) %>% select(Time,Contamination),
                          ymax=sim %>% as_tibble()%>% filter(CurveNum==3) %>% select(Contamination)) %>% rename(Time=ymin.Time,ymin=ymin.Contamination,ymax=Contamination)



ggplot()+
  geom_point(data=df %>% filter(treatment=="Treated" & time <30),aes(x=time,y=conc))+
  geom_line(data=sim %>% as_tibble()%>% filter(CurveNum==2) ,aes(x=Time,y=Contamination))+
  geom_ribbon(data=ribbon_data ,aes(x=Time,ymin=ymin, ymax=ymax), alpha=0.2)+
  # stat_summary(data=df %>% filter(treatment=="Treated"),aes(x=time,y=conc),fun.data = "mean_cl_boot", colour = "red", size = 2)+
  # scale_y_continuous(trans="log")+
  scale_y_continuous(trans="log")+
  xlab("Time (h)")+
  ylab("Concentration")+
  hrbrthemes::theme_ipsum()+
  theme(legend.position = "none")








