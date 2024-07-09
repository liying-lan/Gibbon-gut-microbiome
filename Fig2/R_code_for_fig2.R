####This code is for the random forest analysis and the dynamic of gut bacterial diversity and composition over time

####history clean####
rm(list=ls())
####please setting the working directory manually####

####load package###
require(ggplot2)
require(tidyverse)
library(doParallel)
library(patchwork)
library(dplyr)
library(randomForest)
library(rfPermute)
library(data.table)
library(ggpubr)

####load data#####
raw_group <- read.csv("ENV.csv")

total_div <- read.csv("total_div.csv")

raw_ASV <- read.csv("asv_tab.csv",
                    header=T, check.names=F ,row.names=1, fileEncoding = 'GBK')

taxonomy <- read.csv("taxonomy.csv")



####Fig.2A--random forest model######
#####1. explain variation calculate####

ENV <- raw_group %>% 
  na.omit(Leaf_pro:shannon_entropy)

set.seed(1234)#设置种子，确保结果一致
RF<- randomForest(shannon_entropy~., ENV[,-1:-2], importance = T)
RFs<- rfPermute(shannon_entropy~., ENV[,-1:-2], nperm=99, ntree=501)
RF;importance(RFs)[,1:2]

df1<- subset(ENV, Species == "Hoolock_tianxing")[,-1:-2]
set.seed(1234)
RF1<- randomForest(shannon_entropy~., df1, importance = T)
RF1s<- rfPermute(shannon_entropy~., df1, nperm=99, ntree=501)
RF1;importance(RF1s)[,1:2]

df2<- subset(ENV, Species == "Nomascus_concolor")[,-1:-2]
set.seed(1234)
RF2<- randomForest(shannon_entropy~., df2, importance = T)
RF2s<- rfPermute(shannon_entropy~., df2, nperm=99, ntree=501)
RF2;importance(RF2s)[,1:2]


df3<- subset(ENV, Species == "Nomascus_hainanus")[,-1:-2]
set.seed(1234)
RF3<- randomForest(shannon_entropy~., df3, importance = T)
RF3s<- rfPermute(shannon_entropy~., df3, nperm=99, ntree=501)
RF3;importance(RF3s)[,1:2]

bar_data <- data.frame(variable =c("Whole", 
                                   "H.tianxing",
                                   "N.concolor",
                                   "N.hainanus"),
                 Exp = c(14.13,11.27,17.68,8.89)) %>%  #from %Var explained above
  dplyr::mutate(variable=factor(variable,
                                level=c("Whole", 
                                        "H.tianxing",
                                        "N.concolor",
                                        "N.hainanus")))

#####2. bar_plot figure####

p1A_bar_plot <- ggplot(bar_data, aes(variable, Exp))+
  geom_bar(stat = "identity", fill = "steelblue")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25))+
  theme_classic(base_line_size = 0.75)+
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, color = "black",
                                   size = 12, vjust = 0.6), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 12))+
  labs(title = "Explained variation (%)", size = 12)


#####3. cor-relation calculate####

r<- data.frame(ENV = colnames(ENV)[3:8],
               Whole = (cor(ENV[,-1:-2],method = "spearman"))[7,1:6],
               H.tianxing = (cor(df1))[7,1:6],
               N.concolor = (cor(df2))[7,1:6],
               N.hainanus = (cor(df3))[7,1:6])%>%
  melt(id = "ENV", value.name = "Correlation") %>% 
  dplyr::mutate(ENV=factor(ENV))

circle<- data.frame(ENV = colnames(ENV)[3:8]) %>%
  left_join(data.frame(ENV = row.names(importance(RFs)),
                       Whole = ifelse(importance(RFs)[,2]<0.05,
                                      importance(RFs)[,1], NA))) %>%
  left_join(data.frame(ENV = row.names(importance(RF1s)),
                       H.tianxing = ifelse(importance(RF1s)[,2]<0.05,
                                     importance(RF1s)[,1], NA))) %>%
  left_join(data.frame(ENV = row.names(importance(RF2s)),
                       N.concolor = ifelse(importance(RF2s)[,2]<0.05,
                                     importance(RF2s)[,1], NA))) %>%
  left_join(data.frame(ENV = row.names(importance(RF3s)),
                       N.hainanus = ifelse(importance(RF3s)[,2]<0.05,
                                     importance(RF3s)[,1], NA))) %>%
  tidyr::pivot_longer(cols = Whole:last_col(),
                      names_to = "variable",
                      values_to = "Importance")
  

#####4. heatmap plot making####
p1A_heatmap_plot <- ggplot()+
  geom_tile(data = r, aes(x = variable,  y = ENV, fill = Correlation))+ #热图
  scale_fill_gradientn(colors = c('#2D6DB1', 'white', '#DC1623'),
                       limit = c(-1, 1))+
  geom_point(data = circle, aes(x = variable, y = ENV, 
                                size = Importance), shape = 21)+ #圆圈大小
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, color = "black",
                                   size = 12, vjust = 0.6), 
        axis.text.y = element_text(color = 'black', size = 12),
        legend.title = element_text(size = 10),
        legend.position = "bottom")+
  labs(y = '', x = '')


####Fig.2B--fruit consumption over year######
#####1. data clean####

plot_data <- ENV %>% 
  dplyr::select(Species,month,Leaf_pro,Fruit_pro) %>% 
  dplyr::group_by(Species,month) %>% 
  dplyr::mutate(Leaf_pro=mean(Leaf_pro),
                Fruit_pro=mean(Fruit_pro)) %>% 
  distinct() %>% 
  dplyr::mutate(month=factor(month,levels=c("1","2","3","4","5","6","7",
                                            "8","9","10","11","12"))) 

#####2. point figure plot ####
p2B_point_plot<-ggplot(plot_data, aes( x = month, 
                                       y = Fruit_pro, 
                                       color=Leaf_pro,
                                       group=Leaf_pro)) +   
  geom_point(aes(color=Leaf_pro), 
             size=5)+  
  geom_line(aes(x=month, y=Fruit_pro,group = Species),
            stat="identity",
            color= "black",
            size=1) +  
  facet_wrap(vars(Species),  ncol = 3) +  
  scale_colour_gradientn(colours = c("orange", "#edf8e9", "#006d2c"))+
  theme_bw()+
  labs(x="Month",y= "Fruit eating proportion (%)")


####Fig.2C_2E--######
#####1. function build####

plot.function <- function(input_tab,
                          select_index) {
  plot_data <- input_tab %>% 
    dplyr::mutate(month=factor(month,
                               levels=c("1","2","3","4","5","6",
                                        "7","8","9","10","11","12")),
                  Species=factor(Species,
                                 levels=c("Hoolock_tianxing",
                                          "Nomascus_concolor",
                                          "Nomascus_hainanus"))) %>%
    dplyr::group_by(Species,month,test_index) %>% 
    dplyr::summarise(mean_value=mean(value),
                     sd_value=sd(value),
                     Leaf_pro=mean(Leaf_pro)) %>%
    dplyr::filter(test_index==select_index)
  
  
  plot <- ggplot(plot_data, aes( x = month, y = mean_value, color=Leaf_pro, fill=Leaf_pro)) +   
    geom_point(aes( fill = Leaf_pro), 
               size=5)+  
    geom_line(aes(x=month, y=mean_value,group = Species),
              stat="identity",
              color= "black",
              size=1) +  
    facet_wrap(vars(Species),  ncol = 3) +  
    scale_fill_gradientn(colours = c("orange", "#edf8e9", "#006d2c"))+
    scale_color_gradientn(colours = c("orange", "#edf8e9", "#006d2c"))+
    labs(x="Month",y= select_index)+  
    theme_bw()
  
  return(plot)
}

#####2. plot####

index_list <- unique(total_div$test_index)


for (select_index in index_list) {
  sel_plot <- plot.function(input_tab = total_div,
                            select_index = select_index)
  assign(paste0(select_index,"_plot"),sel_plot)
}

Faith_pd_plot/Nestedness_plot/Turnover_plot

####Fig.2F-2I--######
#####1. data clean####

design <- raw_group %>% 
  dplyr::select(Species,month,Leaf_pro) %>%
  na.omit() %>% 
  dplyr::group_by(Species,month) %>% 
  dplyr::mutate(Leaf_pro=mean(Leaf_pro))%>% 
  tidyr::unite("Species_month",
               Species:month,
               na.rm = T,
               remove=T)

group_list <- unique(design$Species_month)

data <- foreach(k = group_list, 
                     .combine = rbind,
                     .packages = c("tidyverse",
                                   "tidymodels",
                                   "rstatix",
                                   "doParallel")) %do% {
                                     #k <-group_list[1]
                                     sub_group <- raw_group %>% 
                                       dplyr::select(Species,month,everything()) %>% 
                                       tidyr::unite("Species_month",
                                                    Species:month,
                                                    na.rm = T,
                                                    remove=T) %>% 
                                       dplyr::filter(Species_month==k)
                                     
                                     sample_id <- sub_group$ID
                                     
                                     test_tab <- raw_ASV[,sample_id] %>% 
                                         data.frame() %>% 
                                         dplyr::mutate(row_all=rowSums(.)) %>% 
                                         dplyr::filter(row_all!=0) %>% 
                                         dplyr::mutate(ASV_ID=rownames(.)) %>% 
                                         dplyr::left_join(.,taxonomy)
                                     
                                     asv_num <- sum(test_tab$row_all)
                                     
                                     data <- data.frame(phylum_num=length(unique(test_tab$Phylum)),
                                                        family_num=length(unique(test_tab$Family)),
                                                        ASV_num=length(unique(test_tab$ASV_ID)),
                                                        Total_ASV=asv_num,
                                                        k) %>% 
                                         dplyr::rename(Species_month=k) %>% 
                                         dplyr::left_join(.,design) %>%    
                                         tidyr::separate(col = Species_month,
                                                         into = c("s1","s2","month")) %>% 
                                         tidyr::unite("Species",
                                                      s1:s2,
                                                      remove=T) %>%    
                                         dplyr::distinct(Species,month,.keep_all = T)
                                    
                                     return(data)
                                     
                                   }

#####2. plot_function####

plot.function <- function(input_tab) {
  plot <- ggplot(input_tab, aes( x = month, y = y_value, color=Leaf_pro, fill=Leaf_pro)) +  
    geom_point(aes( fill = Leaf_pro), 
               size=5)+  
    geom_line(aes(x=month, y=y_value,group = Species),
              stat="identity",
              color= "black",
              size=1) +  
    facet_wrap(vars(Species),  ncol = 3) +  
    scale_fill_gradientn(colours = c("orange", "#edf8e9", "#006d2c"))+
    scale_color_gradientn(colours = c("orange", "#edf8e9", "#006d2c"))+
    theme_bw()+
    xlab("Month")
  
  return(plot)
  
}

#####3. plot_begin####
colnames(data)

y_values <- colnames(data)[1:4]

colnames(data)

for (i in y_values) {
  #i = y_values[1]
  plot_data <- data %>% 
    dplyr::select(Species,month,Leaf_pro,i) %>% 
    dplyr::mutate(month=factor(month,level=c("1","2","3","4","5","6",
                                             "7","8","9","10","11","12"))) %>% 
    dplyr::rename(y_value=colnames(.)[4])
  y_value=paste0(i)
  sel_plot <- plot.function(input_tab = plot_data)
  assign(paste0(i,"_plot"),sel_plot)
}

print(phylum_num_plot)
print(family_num_plot)
print(ASV_num_plot)
print(Total_ASV_plot)


