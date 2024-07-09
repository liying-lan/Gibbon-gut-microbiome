####This code is for showing the dynamic of community assembly process

####history clean####
rm(list=ls())
####please setting the working directory manually####

####load packages#####
library(dplyr)
library(doParallel)
library(ggpubr)


####load data#####

raw_group <- read.csv("metadata.csv")

dataset <- read.csv("icamp_out.csv")

####1.data clean######
spe_month <- unique(raw_group$spe_month)

sel_data <- foreach(k = spe_month, 
                      .combine = rbind,
                      .packages = c("tidyverse",
                                    "tidymodels",
                                    "rstatix",
                                    "doParallel")) %do% {
                                      #k <- spe_month[2]
                                      
                                      colnames(raw_group)
                                      
                                      group <- subset(raw_group, 
                                                           spe_month==k) %>% 
                                        dplyr::rename(samp1=ID)
                                      
                                      ID <- group$samp1
                                      
                                      sub_data <- subset(dataset, 
                                                         samp1 %in% ID & samp1 %in% ID) %>%
                                        left_join(.,group) 
                                      
                                      return(sub_data)
                                      
                                    }


point_plot_data <-  sel_data %>% 
  dplyr::group_by(spe_month) %>% 
  dplyr::mutate(Leaf_pro=mean(Leaf_pro),
                HeS=mean(HeS),
                HoS=mean(HoS),
                DL=mean(DL),
                HD=mean(HD),
                DR=mean(DR),
                Deterministic.process=(HeS+HoS),
                Stochastic.process=(1-Deterministic.process)) %>% 
  dplyr::select(-samp1,-samp2) %>% 
  distinct() %>% 
  tidyr::separate(col=spe_month,
                  into=c("s1","s2","Group","month")) %>% 
  tidyr::unite("Species",s1:s2)%>% 
  dplyr::mutate(Species=factor(Species,
                               levels = c("Hoolock_tianxing",
                                          "Nomascus_concolor",
                                          "Nomascus_hainanus")))
####Fig.3A_cicle plot######
#####1. plot data####
venn_plot_data <- point_plot_data %>%  
  dplyr::select(Leaf_level, HeS,HoS,
                DL,HD,DR) %>% 
  dplyr::group_by(Leaf_level) %>% 
  dplyr::mutate(HeS=mean(HeS)*100,
                HoS=mean(HoS)*100,
                DL=mean(DL)*100,
                HD=mean(HD)*100,
                DR=mean(DR)*100) %>% distinct() %>% 
  tidyr::pivot_longer(cols = HeS:DR,
                      names_to = "process",
                      values_to = "percent_value") %>% 
  dplyr::mutate(process=factor(process,
                               levels=c("HeS","HoS","DL","HD","DR"))) 

#####2. figure making####
my_pal <- c("#C08F55","#E5DA48","#819FB5","#C3DFF2","#D4D4D5")


total_venn <- ggplot(venn_plot_data, aes(x = Leaf_level, y = percent_value, fill = process)) +
  geom_col(width = 0.7, color = "black", show.legend = TRUE) +  
  scale_x_discrete(limits = c(" ", "A", "B", "C", "D", "E", "F", "G", "H", "I")) +
  coord_polar("y") +
  scale_fill_manual(values = my_pal) +
  theme_minimal(base_size = 14) 


####Fig.3B_bar plot######
#####1. plot data####
bar_data <- sel_data %>% 
  dplyr::select(HeS:Species) %>% 
  group_by(Species) %>%
  dplyr::mutate(HeS=sum(HeS),HoS=sum(HoS),DL=sum(DL), HD=sum(HD),DR=sum(DR)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(process_total=HeS+HoS+DL+HD+DR,
                HeS=HeS/process_total,
                HoS=HoS/process_total,
                DL=DL/process_total,
                HD=HD/process_total,
                DR=DR/process_total) %>% 
  dplyr::select(Species,HoS,HeS,DL,HD,DR) %>% 
  tidyr::pivot_longer(.,
                      cols=HoS:DR,
                      names_to = "process",
                      values_to = "persentage") %>% 
  dplyr::mutate(process=factor(process,
                               levels=c("HeS","HoS","DL","HD","DR")),
                Species=factor(Species,
                               levels = c("Hoolock_tianxing",
                                          "Nomascus_concolor",
                                          "Nomascus_hainanus")))


#####2. figure making####

my_pal <- c("#C08F55","#E5DA48","#819FB5","#C3DFF2","#D4D4D5")

p2_bar_plot <- ggplot(bar_data,aes(x=Species,y=persentage,fill=process))+
  geom_bar(position = "stack",width=0.5,stat="identity",colour='black')+
  theme_bw()+
  scale_fill_manual(values = my_pal)+
  labs(y="Percentage",fill="Assemble process")

####Fig.3G-F_point plot######
#####1. plot function####
plot.function <- function(input_tab) {
    plot <- ggplot(input_tab,
                 aes(Leaf_pro,y_value,
                     color = Species,fill = Species))+
      geom_point(size=3)+
      geom_smooth(method = 'lm', formula = y ~ x, se = T)+ 
      stat_cor(size=4)+
      theme_bw()+
      scale_color_manual(values = c("#9ECAE1","#CAB2D6FF", "#7E6148B1"))+
      scale_fill_manual(values = c("#9ECAE1","#CAB2D6FF", "#7E6148B1"))+
      geom_hline(yintercept = 0.5,linetype=2,linewidth=1,alpha=0.7)+
      theme(axis.text.x=element_text(colour="black",family="Times",size = 20,face="plain"),
            axis.text.y=element_text(colour="black",family="Times",size=20, hjust = 1), 
            axis.title.y=element_text(colour="black",family="Times",size = 20,face="plain"), 
            axis.title.x=element_blank(), 
            panel.background = element_rect(color = '#636363', 
                                        fill = 'transparent'),
            legend.title = element_text(size = 15,face = "bold",
                                    vjust = 1, hjust = 0), 
            legend.text = element_text(size = 12),
            legend.key.size=unit(.5,'cm'),
            legend.position = "none",
            strip.text = element_text(size = 15),
            strip.text.x = element_text(size = 15, colour = "black",
                                        face = "bold.italic"),
            panel.grid.minor = element_blank())+ylim(0,1)+xlim(0,1)

}

#####2. plot_begin#### 
y_values <- colnames(point_plot_data)[c(12,13,2,3)]

colnames(point_plot_data)

for (i in y_values) {
  #i = y_values[1]
  plot_data <- point_plot_data %>% 
    dplyr::select(Species,month,Leaf_pro,i) %>% 
    dplyr::mutate(month=factor(month,level=c("1","2","3","4","5","6",
                                             "7","8","9","10","11","12"))) %>% 
    dplyr::rename(y_value=colnames(.)[4])
  y_value=paste0(i)
  sel_plot <- plot.function(input_tab = plot_data)
  assign(paste0(i,"_plot"),sel_plot)
}

print(Deterministic.process_plot)
print(Stochastic.process_plot)
print(HoS_plot)
print(DL_plot)




####The end######
