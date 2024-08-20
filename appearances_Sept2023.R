library(tidyverse)

### Real data ###

data <- read.csv("~/Google Drive/My Drive/MPIIB/first_appearances/Simulations_Manuela/real_data_clones_all.csv")
count<-data %>%
  group_by(Timepoint,first_appearance) %>%
  summarise(count=n())

count<-subset(count,first_appearance==1)
count$Timepoint<-gsub("13","May13",count$Timepoint)
count$Timepoint <- factor(count$Timepoint,levels=c("0","1","2","3","4","5","6","7","8","9","10","11"))
ggplot(count,aes(factor(Timepoint),count,group=first_appearance))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle=0,size=10),
        axis.text.y = element_text(size=10),
        legend.position="bottom",
        panel.background = element_rect(fill=NA,size=1,color="black"))+
  labs(x="Timepoint",y="first_appearances_per_timepoint")
  

## Contingency tables

myTab <- table("first_appearance" = data$first_appearance, "peak" = data$is_peak)
FA_real<-data %>% group_by(first_appearance) %>% summarise(count=n()) ## Number is 907

first_appearance_per_100_peaks_observed <- round(myTab[2,2]/(myTab_realDat[1,2]+myTab_realDat[2,2])*100)

first_appearance_per_100_peaks_observed <- round(myTab[2,2]/(myTab[1,2]+myTab[2,2])*100)

#fisher <- fisher.test(myTab2)

FA_real1<-real1 %>% group_by(first_appearance) %>% summarise(count=n()) ## Number is 660
FA_real<-real %>% group_by(first_appearance) %>% summarise(count=n()) ## Number is 717

first_appearance_per_100_peaks_observed <-round(myTab2[2,2]/(myTab2[2,2]+myTab2[2,1])*100)


### Simulated data first attempt ###

sim1<-read.csv("~/Google Drive/My Drive/MPIIB/first_appearances/Simulations_Manuela/simulations_weighted.csv.gz",h=T)
sim2<-read.csv("~/Google Drive/My Drive/MPIIB/first_appearances/Simulations_Manuela/simulations_unweighted.csv.gz",h=T)

FA_sim1<-sim1 %>% group_by(simulation,first_appearance) %>% summarise(count=n())
FA_sim2<-sim2 %>% group_by(simulation,first_appearance) %>% summarise(count=n())


simulations<-sim1 %>% group_split(simulation)
first_appearance_per_100_peaks <- simulation_number <- c()
for(k in 1:length(simulations)){
  mytab <- table("first_appearance" = simulations[[k]]$first_appearance, "peak" = simulations[[k]]$is_peak)
  simulation_number[k] <- k
  first_appearance_per_100_peaks[k] <- round(mytab[2,2]/(mytab[1,2]+mytab[2,2])*100)
}

weighted <- data.frame(simulation = simulation_number, number_of_first_appearances_100_peaks = first_appearance_per_100_peaks)


simulations<-sim2 %>% group_split(simulation)
first_appearance_per_100_peaks <- simulation_number <- c()
for(k in 1:length(simulations)){
  mytab <- table("first_appearance" = simulations[[k]]$first_appearance, "peak" = simulations[[k]]$is_peak)
  simulation_number[k] <- k
  first_appearance_per_100_peaks[k] <- round(mytab[2,2]/(mytab[1,2]+mytab[2,2])*100)
}

unweighted <- data.frame(simulation = simulation_number, number_of_first_appearances_100_peaks = first_appearance_per_100_peaks)

weighted$type<-"weighted"
unweighted$type<-"unweighted"

all<-rbind(weighted,unweighted)


ggplot(all, aes(x = number_of_first_appearances_100_peaks)) +
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="black", fill="gray", alpha = 0.4) +
  geom_density(alpha=.4) + 
  #scale_x_continuous(breaks = c(61:69)) +
  #geom_vline(xintercept = first_appearance_per_100_peaks_observed, linetype = "dashed", color = "red") +
  theme_bw()+
  facet_wrap(~type)
  
  
  
  
  

sim3<-read.csv("~/Desktop/MPIIB/simulations_weighted_v2.csv.gz",h=T)
sim4<-read.csv("~/Desktop/MPIIB/simulations_unweighted_v2.csv.gz",h=T)

FA_sim3<-sim3 %>% group_by(simulation,first_appearance) %>% summarise(count=n())
FA_sim4<-sim4 %>% group_by(simulation,first_appearance) %>% summarise(count=n())


simulations<-sim4 %>% group_split(simulation)
first_appearances_in_peaks <- simulation_number <- c()
for(k in 1:length(simulations)){
  mytab <- table("first_appearance" = simulations[[k]]$first_appearance, "peak" = simulations[[k]]$is_peak)
  simulation_number[k] <- k
  first_appearances_in_peaks[k] <- round(mytab[2,1]/(mytab[2,2]+mytab[2,1])*100)
}

weighted <- data.frame(simulation = simulation_number, proportion_of_first_appearances_peaks = first_appearances_in_peaks)
#weighted <- data.frame(simulation = simulation_number, number_of_first_appearances_100_peaks = first_appearance_per_100_peaks)

#weighted$type<-"weighted"
#unweighted$type<-"unweighted"

#all<-rbind(weighted,unweighted)


ggplot(weighted, aes(x = proportion_of_first_appearances_peaks)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="black", fill="gray", alpha = 0.4) +
  #scale_x_continuous(breaks = c(61:69)) +
  geom_vline(xintercept = first_appearance_per_100_peaks_observed, linetype = "dashed", color = "red") +
  theme_bw()+
  labs(x="First appearances in 100 peaks")



