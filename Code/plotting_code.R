
# Read in cleaningExptBeth.csv from LabData

df<-read.csv("LabData/cleaningExptBeth.csv")

#Plot the data with ggplot2


df$TimeAfterCleaning <- as.numeric(df$TimeAfterCleaning)

df%>%
    filter(Count > 0) %>%
    ggplot(aes(x=TimeAfterCleaning,y=Count,color=CleaningType))+
    #geom_point(alpha=0.3)+
    #Add errorbars to the points using stat_summary
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_se, geom = "errorbar")+
    scale_color_brewer(palette="Set1")+
    scale_y_continuous(trans="log10")+
    facet_grid(CleaningType~.,scales="free")+
    hrbrthemes::theme_ipsum()+
    labs(x="Time after cleaning (hours)",y="CFU",title="Cleaning Experiment")+
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
