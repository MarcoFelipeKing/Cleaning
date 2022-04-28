
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