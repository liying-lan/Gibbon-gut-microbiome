####This code is to make bar plot of contribution of each group and bin####

####history clean####
rm(list=ls())
####please setting the working directory manually####

####load package####
library(tidyverse)
library(dplyr)


####load data#####
bin_data <- read.csv("Test.BinContributeToProcess_EachGroup_3spe.csv") %>% 
  dplyr::select(-leaf_group)


####1.data clean#####

df_plot1 <- bin_data %>% 
  #dplyr::filter(Species=="Hoolock_tianxing") %>% 
  dplyr::filter(Process=="HeS"|Process=="HoS"|Process=="DL") %>% 
  dplyr::group_by(Species,Process) %>% 
  dplyr::summarise(across(Bin1:Bin11, ~ sum(., na.rm = TRUE))) %>% 
  tidyr::pivot_longer(cols =Bin1:Bin11,
                      names_to = "bin_name",
                      values_to = "value") %>% 
  ddply(.,'Process',transform,percent_value=value/sum(value)*100) %>% 
  dplyr::mutate(Process=factor(Process,
                               level=c("HeS","HoS","DL")),
                bin_name=factor(bin_name,
                                level=c("Bin1","Bin2","Bin3","Bin4","Bin5","Bin6",
                                        "Bin7","Bin8","Bin9","Bin10","Bin11")))


####2.plot making#####
my_pal <- c("#A6CEE3FF","#FDBF6FFF",	"#2986CC",	"#b2df8a","#E06666",
            "#FEC211","#CAB2D6FF","#6A3D9AFF", "#FF7F00FF","#B45F06", "#7E6148B1")


P1 <- ggplot(df_plot1,aes(x=Process,y=percent_value,fill=bin_name))+
  geom_bar(position = "fill",width=0.5,stat="identity",colour='black')+
  facet_wrap(vars(Species),ncol=1) +
  theme_bw()+
  scale_fill_manual(values =my_pal)+
  labs(y="Percentage (%)",fill="Bin_name")+
  scale_y_continuous(position = "right")

