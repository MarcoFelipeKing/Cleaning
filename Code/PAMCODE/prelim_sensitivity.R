# sensitivity analysis


this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    


SIM <- c("DPLPNS", "DPLPS","DPHPNS","DPHPS", "DPallNS", "DPallS",
         "DOLPNS","DOLPS","DOHPNS","DOHPS", "DOallNS", "DOallS",
         "DGLPNS","DGLPS","DGHPNS","DGHPS", "DGallNS", "DGallS",
         "control")   

NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

for(j in 1:NUM.SIM){
  
  sim.num <- j; sim.name <- SIM[j]
  
  if(sim.name=="DPLPNS"){if(dir.exists("DPLPNS")==FALSE){dir.create("DPLPNS"); setwd("DPLPNS")}else{setwd("DPLPNS")}}
  if(sim.name=="DPLPS"){if(dir.exists("DPLPS")==FALSE){dir.create("DPLPS"); setwd("DPLPS")}else{setwd("DPLPS")}}
  if(sim.name=="DPHPNS"){if(dir.exists("DPHPNS")==FALSE){dir.create("DPHPNS"); setwd("DPHPNS")}else{setwd("DPHPNS")}}
  if(sim.name=="DPHPS"){if(dir.exists("DPHPS")==FALSE){dir.create("DPHPS"); setwd("DPHPS")}else{setwd("DPHPS")}}
  if(sim.name=="DPallNS"){if(dir.exists("DPallNS")==FALSE){dir.create("DPallNS"); setwd("DPallNS")}else{setwd("DPallNS")}}
  if(sim.name=="DPallS"){if(dir.exists("DPallS")==FALSE){dir.create("DPallS"); setwd("DPallS")}else{setwd("DPallS")}}
  
  if(sim.name=="DOLPNS"){if(dir.exists("DOLPNS")==FALSE){dir.create("DOLPNS"); setwd("DOLPNS")}else{setwd("DOLPNS")}}
  if(sim.name=="DOLPS"){if(dir.exists("DOLPS")==FALSE){dir.create("DOLPS"); setwd("DOLPS")}else{setwd("DOLPS")}}
  if(sim.name=="DOHPNS"){if(dir.exists("DOHPNS")==FALSE){dir.create("DOHPNS"); setwd("DOHPNS")}else{setwd("DOHPNS")}}
  if(sim.name=="DOHPS"){if(dir.exists("DOHPS")==FALSE){dir.create("DOHPS"); setwd("DOHPS")}else{setwd("DOHPS")}}
  if(sim.name=="DOallNS"){if(dir.exists("DOallNS")==FALSE){dir.create("DOallNS"); setwd("DOallNS")}else{setwd("DOallNS")}}
  if(sim.name=="DOallS"){if(dir.exists("DOallS")==FALSE){dir.create("DOallS"); setwd("DOallS")}else{setwd("DOallS")}}
  
  if(sim.name=="DGLPNS"){if(dir.exists("DGLPNS")==FALSE){dir.create("DGLPNS"); setwd("DGLPNS")}else{setwd("DGLPNS")}}
  if(sim.name=="DGLPS"){if(dir.exists("DGLPS")==FALSE){dir.create("DGLPS"); setwd("DGLPS")}else{setwd("DGLPS")}}
  if(sim.name=="DGHPNS"){if(dir.exists("DGHPNS")==FALSE){dir.create("DGHPNS"); setwd("DGHPNS")}else{setwd("DGHPNS")}}
  if(sim.name=="DGHPS"){if(dir.exists("DGHPS")==FALSE){dir.create("DGHPS"); setwd("DGHPS")}else{setwd("DGHPS")}}
  if(sim.name=="DGallNS"){if(dir.exists("DGallNS")==FALSE){dir.create("DGallNS"); setwd("DGallNS")}else{setwd("DGallNS")}}
  if(sim.name=="DGallS"){if(dir.exists("DGallS")==FALSE){dir.create("DGallS"); setwd("DGallS")}else{setwd("DGallS")}}
  
  if(sim.name=="control"){if(dir.exists("control")==FALSE){dir.create("control"); setwd("control")}else{setwd("control")}}
  
  
  numiterations<-500 
  
  IV<-readRDS(file=sprintf("%s.IV.exposure.frame.rds",sim.name))
  Rounds<-readRDS(file=sprintf("%s.Rounds.exposure.frame.rds",sim.name))
  Obs<-readRDS(file=sprintf("%s.Obs.exposure.frame.rds",sim.name))
  

  for(i in 1:numiterations){
    IVframe<-IV[[i]]
    Roundsframe<-Rounds[[i]]
    Obsframe<-Obs[[i]]
   
    if (i ==1){
      
      RNA<-c(mean(IVframe$RNAinfectiousall[!is.na(IVframe$RNAinfectiousall)]),mean(Roundsframe$RNAinfectiousall[!is.na(Roundsframe$RNAinfectiousall)]),mean(Obsframe$RNAinfectiousall[!is.na(Obsframe$RNAinfectiousall)]))
      SM<-c(mean(IVframe$SMall[IVframe$SMall>0]),mean(Roundsframe$SMall[Roundsframe$SMall>0]),mean(Obsframe$SMall[Obsframe$SMall>0]))
      TEmouth<-c(mean(IVframe$TE.tempall[IVframe$TE.tempall>0]),mean(Roundsframe$TE.tempall[Roundsframe$TE.tempall>0]),mean(Obsframe$TE.tempall[Obsframe$TE.tempall>0]))
      Ah.dose<-c(mean(IVframe$AH.doseall[IVframe$AH.doseall>0]),mean(Roundsframe$AH.doseall[Roundsframe$AH.doseall>0]),mean(Obsframe$AH.doseall[Obsframe$AH.doseall>0]))
      beta.dose<-c(mean(IVframe$beta.dose.all),mean(Roundsframe$beta.dose.all),mean(Obsframe$beta.dose.all))
      alpha<-c(mean(IVframe$alphaall),mean(Roundsframe$alphaall),mean(Obsframe$alphaall))
      lambda<-c(mean(IVframe$lambda),mean(Roundsframe$lambda),mean(Obsframe$lambda))
      beta<-c(mean(IVframe$beta),mean(Roundsframe$beta),mean(Obsframe$beta))
      rnainfect<-c(mean(IVframe$RNAinfectiousall),mean(Roundsframe$RNAinfectiousall),mean(Obsframe$RNAinfectiousall))
      duration<-c(mean(IVframe$duration),mean(Roundsframe$duration),mean(Obsframe$duration))
      SH<-c(mean(IVframe$SH),mean(Roundsframe$SH),mean(Obsframe$SH))
      surfconc<-c(mean(IVframe$surfconc),mean(Roundsframe$surfconc),mean(Obsframe$surfconc))
      k.sall<-c(mean(IVframe$k.sall),mean(Roundsframe$k.sall),mean(Obsframe$k.sall))
      k.hall<-c(mean(IVframe$k.hall),mean(Roundsframe$k.hall),mean(Obsframe$k.hall))
      infect<-c(max(IVframe$infect),max(Roundsframe$infect),max(Obsframe$infect))
      care<-c("IV","Rounds","Obs")
      
    }else{
      RNAtemp<-c(mean(IVframe$RNAinfectiousall[!is.na(IVframe$RNAinfectiousall)]),mean(Roundsframe$RNAinfectiousall[!is.na(Roundsframe$RNAinfectiousall)]),mean(Obsframe$RNAinfectiousall[!is.na(Obsframe$RNAinfectiousall)]))
      SMtemp<-c(mean(IVframe$SMall[IVframe$SMall>0]),mean(Roundsframe$SMall[Roundsframe$SMall>0]),mean(Obsframe$SMall[Obsframe$SMall>0]))
      TEmouthtemp<-c(mean(IVframe$TE.tempall[IVframe$TE.tempall>0]),mean(Roundsframe$TE.tempall[Roundsframe$TE.tempall>0]),mean(Obsframe$TE.tempall[Obsframe$TE.tempall>0]))
      Ah.dosetemp<-c(mean(IVframe$AH.doseall[IVframe$AH.doseall>0]),mean(Roundsframe$AH.doseall[Roundsframe$AH.doseall>0]),mean(Obsframe$AH.doseall[Obsframe$AH.doseall>0]))
      beta.dosetemp<-c(mean(IVframe$beta.dose.all),mean(Roundsframe$beta.dose.all),mean(Obsframe$beta.dose.all))
      alphatemp<-c(mean(IVframe$alphaall),mean(Roundsframe$alphaall),mean(Obsframe$alphaall))
      lambdatemp<-c(mean(IVframe$lambda),mean(Roundsframe$lambda),mean(Obsframe$lambda))
      betatemp<-c(mean(IVframe$beta),mean(Roundsframe$beta),mean(Obsframe$beta))
      rnainfecttemp<-c(mean(IVframe$RNAinfectiousall),mean(Roundsframe$RNAinfectiousall),mean(Obsframe$RNAinfectiousall))
      durationtemp<-c(mean(IVframe$duration),mean(Roundsframe$duration),mean(Obsframe$duration))
      SHtemp<-c(mean(IVframe$SH),mean(Roundsframe$SH),mean(Obsframe$SH))
      surfconctemp<-c(mean(IVframe$surfconc),mean(Roundsframe$surfconc),mean(Obsframe$surfconc))
      k.salltemp<-c(mean(IVframe$k.sall),mean(Roundsframe$k.sall),mean(Obsframe$k.sall))
      k.halltemp<-c(mean(IVframe$k.hall),mean(Roundsframe$k.hall),mean(Obsframe$k.hall))
      infecttemp<-c(max(IVframe$infect),max(Roundsframe$infect),max(Obsframe$infect))
      caretemp<-c("IV","Rounds","Obs")
      
      RNA<-c(RNA,RNAtemp)
      SM<-c(SM,SMtemp)
      TEmouth<-c(TEmouth,TEmouthtemp)
      Ah.dose<-c(Ah.dose,Ah.dosetemp)
      beta.dose<-c(beta.dose,beta.dosetemp)
      alpha<-c(alpha,alphatemp)
      rnainfect<-c(rnainfect,rnainfecttemp)
      lambda<-c(lambda,lambdatemp)
      beta<-c(beta,betatemp)
      duration<-c(duration,durationtemp)
      SH<-c(SH,SHtemp)
      surfconc<-c(surfconc,surfconctemp)
      k.sall<-c(k.sall,k.salltemp)
      k.hall<-c(k.hall,k.halltemp)
      infect<-c(infect,infecttemp)
      care<-c(care,caretemp)
    }
    
 
  }
  
  if (j==1){
    scenario<-rep(sprintf("%s",sim.name),length(lambda))
    frame<-data.frame(lambda=lambda,beta=beta,duration=duration,SH=SH,surfconc=surfconc,k.sall=k.sall,k.hall=k.hall,
                      infect=infect,care=care,scenario=scenario,beta.dose=beta.dose,alpha=alpha,
                      SM=SM,TEmouth=TEmouth,Ah.dose=Ah.dose,RNA=RNA)
    
  }else{
    scenariotemp<-rep(sprintf("%s",sim.name),length(lambda))
    frametemp<-data.frame(lambda=lambda,beta=beta,duration=duration,SH=SH,surfconc=surfconc,k.sall=k.sall,k.hall=k.hall,
                          infect=infect,care=care,scenario=scenariotemp,beta.dose=beta.dose,alpha=alpha,
                          SM=SM,TEmouth=TEmouth,Ah.dose=Ah.dose,RNA=RNA)
    frame<-rbind(frame,frametemp)
    print(length(frame$lambda)/(j*1500))
          print(j)
  }
  
  #reset directory to parent folder so we can go to correct subfolder within parent folder for next sim run
  setwd(this.dir)
}

require(reshape2)
require(ggplot2)
require(ggpubr)

frametemp<-frame[frame$infect>1e-15,]
framecor = subset(frametemp, select = -c(care,scenario) )

cormat<-round(cor(framecor,method=c("spearman")),2)
melted_cormat<-melt(cormat)
ggplot(data=melted_cormat,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
  geom_text(aes(label = signif(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue") 

cor(frametemp$infect,frametemp$beta.dose,method="spearman")
#cor(frametemp$infect,frametemp$alpha,method="spearman")
#cor(frametemp$infect[!is.na(frametemp$RNA)],frametemp$RNA[!is.na(frametemp$RNA)],method="spearman")
#cor(frametemp$infect[!is.na(frametemp$TEmouth)],frametemp$TEmouth[!is.na(frametemp$TEmouth)],method="spearman")
#cor(frametemp$infect[!is.na(frametemp$SM)],frametemp$SM[!is.na(frametemp$SM)],method="spearman")
#cor(frametemp$infect[!is.na(frametemp$Ah.dose)],frametemp$Ah.dose[!is.na(frametemp$Ah.dose)],method="spearman")

setwd(this.dir)

write.csv(frame,'frame_sensitivity.csv')

#--------------working on tile plots to investigate relationships between------------------------------------------
#---------------input parameters that result in highest infection risks--------------------------------------------------------
#A<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(beta,3),fill=log10(infect)))+
#  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
#  scale_y_continuous(name="Surface-to-Hand Transfer Efficiency")+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

#B<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(SH,3),fill=log10(infect)))+
#  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
#  scale_y_continuous(name="Fraction of Hand Contact Area")+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

#C<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(log10(surfconc),1),fill=log10(infect)))+
#  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
#  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

#D<-ggplot(frame)+geom_tile(aes(x=round(beta,3),y=round(log10(surfconc),1),fill=log10(infect)))+
#  scale_x_continuous(name="Surface-to-Hand Transfer Efficiency")+
#  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

#E<-ggplot(frame)+geom_tile(aes(x=round(SH,3),y=round(log10(surfconc),1),fill=log10(infect)))+
#  scale_x_continuous(name="Fraction of Hand Contact Area")+
#  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

#G<-ggplot(frame)+geom_tile(aes(x=round(SH,3),y=round(beta,3),fill=log10(infect)))+
#  scale_x_continuous(name="Fraction of Hand Contact Area")+
#  scale_y_continuous(name="Surface-to-Hand Transfer Efficiency")+
#  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
#  theme_pubr()

windows()
ggarrange(A,B,D,C,E,G,common.legend = TRUE)



A<-ggplot(frametemp[!is.na(frametemp$RNA),])+geom_point(aes(x=RNA,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Fraction of Infectious RNA")
B<-ggplot(frametemp)+geom_point(aes(x=surfconc,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name=expression("Surface Concentration (viral particles/cm"^2*")"),trans="log10")
C<-ggplot(frametemp)+geom_point(aes(x=TEmouth,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Hand-to-Face Transfer Efficiency")
D<-ggplot(frametemp)+geom_point(aes(x=SM,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="FSA for Hand-to-Face Contact")
E<-ggplot(frametemp)+geom_point(aes(x=Ah.dose,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name=expression("Total Hand Surface Area (cm"^2*")"))
G<-ggplot(frametemp)+geom_point(aes(x=beta,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Surface-to-Hand Transfer Efficiency")
H<-ggplot(frametemp)+geom_point(aes(x=lambda,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")
I<-ggplot(frametemp)+geom_point(aes(x=duration,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Duration (seconds)")
J<-ggplot(frametemp)+geom_point(aes(x=SH,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="FSA of Hand-to-Surface")
K<-ggplot(frametemp)+geom_point(aes(x=k.sall,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Surface Inactivation Constant")
L<-ggplot(frametemp)+geom_point(aes(x=k.hall,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Hand Inactivation Constant")
M<-ggplot(frametemp)+geom_point(aes(x=alpha,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Alpha (Dose-response)",trans="log10")
N<-ggplot(frametemp)+geom_point(aes(x=beta.dose,y=infect))+
  scale_y_continuous(trans="log10",name="Infection Risk")+theme_pubr()+
  scale_x_continuous(name="Beta (Dose-response)",trans="log10")


windows()
ggarrange(A,B,C,D,E,G,H,I,
          J,K,L,M,N,common.legend = TRUE)


