####This code is for the analysis of the composition and diversity of gibbon gut bacterial community

####history clean####
rm(list = ls())

####please setting the working directory manually####

####load packages####
library(tidyverse)###data process
library(ggvenn)###venn figure
library(vegan)###beta div
library(ggsci)


####load data####

raw_ASV <- read.csv("asv_tab.csv",row.names=1) 

group <- read.csv("metadata.csv") 

alpha_div <- read.csv("alpha_div.csv")

raw_phylum <- read.csv("phylum_tab.csv") 

raw_family <- read.csv( "family_tab.csv")


####Fig.1A--venn figure######
#####1. data clean####

new_tab <- raw_ASV %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::mutate(ID=rownames(.)) %>% 
  dplyr::select(ID,everything()) %>% 
  dplyr::left_join(.,group) %>% 
  dplyr::select(ID,Species,everything()) %>% 
  tidyr::pivot_longer(cols = 3:last_col(),
                      names_to = "microbe_name",
                      values_to = "relative_abun")%>% 
  dplyr::group_by(Species,microbe_name) %>% 
  dplyr::summarise(all_ASV=sum(relative_abun))


sep_tab <- within(new_tab, {
  new_count <- NA
  new_count[all_ASV > 0] <- "1"
  new_count[all_ASV == 0] <- "0"
  new_count[is.na(all_ASV)] <- "0"
})%>% 
  dplyr::select(-all_ASV) %>% 
  tidyr::pivot_wider(names_from = microbe_name ,
                     values_from = new_count) 


#####2. list exited ASV  ####
#spe1

sp_name =unique(sep_tab$Species)[1]


spe1_list <- sep_tab %>% 
  dplyr::filter(Species == sp_name) %>% 
  dplyr::select(-Species) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::rename(value=colnames(.)[1]) %>% 
  dplyr::filter(value != 0) %>% 
  dplyr::mutate(asv_name=row.names(.)) %>% 
  dplyr::select(-value)

 
#spe2
sp_name =unique(sep_tab$Species)[2]

spe2_list <- sep_tab %>% 
  dplyr::filter(Species == sp_name) %>% 
  dplyr::select(-Species) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::rename(value=colnames(.)[1]) %>% 
  dplyr::filter(value != 0) %>% 
  dplyr::mutate(asv_name=row.names(.)) %>% 
  dplyr::select(-value)


#spe3
sp_name =unique(sep_tab$Species)[3]


spe3_list <- sep_tab %>% 
  dplyr::filter(Species == sp_name) %>% 
  dplyr::select(-Species) %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::rename(value=colnames(.)[1]) %>% 
  dplyr::filter(value != 0) %>% 
  dplyr::mutate(asv_name=row.names(.)) %>% 
  dplyr::select(-value)

#creat list
venn_data <- list(
  "H. tianxing" = c(spe1_list$asv_name),
  "N. concolor" = c(spe2_list$asv_name),
  "N. hainanus" = c(spe3_list$asv_name)
)

#####3. fig make####
p1_venn <- ggvenn(venn_data, fill_color = c("#9ECAE1", "#CAB2D6", "#7E6148")) +
  #labs(title = "Venn Diagram of Species Genes") +
  theme_classic() +  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(-1, -1, -1, -1), "lines")  
  )


####Fig.1B--number of unique ASVs (bar figure)######
#####1. data clean####

t_ASV <- t(raw_ASV) 

t_ASV[t_ASV > 0] <- "1"

summary_df <- t_ASV %>% 
  data.frame() %>% 
  dplyr::mutate(across(everything(), as.numeric)) %>% 
  dplyr::mutate(ASV_num=rowSums(.)) %>% 
  rownames_to_column(var = "ID") %>% 
  dplyr::select(ID,ASV_num) %>% 
  left_join(.,group) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    mean_ASV_Count = mean(ASV_num),
    se_ASV_Count = sd(ASV_num) / sqrt(n())
  ) 

#####2. fig make####
p2_bar <- ggplot(summary_df, aes(Species, mean_ASV_Count)) +
  geom_col(aes(fill = Species), width = 0.6) +  
  scale_fill_manual(values = c("#9ECAE1", "#CAB2D6FF", "#7E6148B1")) +
  geom_errorbar(data=summary_df,mapping=aes(x = Species,
                                            ymin = mean_ASV_Count-se_ASV_Count, 
                                            ymax = mean_ASV_Count+se_ASV_Count),
                width = 0.1, color = 'black',size=0.8)+ 
  geom_text(aes(y = mean_ASV_Count + 10, label = mean_ASV_Count), 
            vjust = 0, size = 4, color = "black") +  
  theme_classic() +  
  ylab("Number of unique ASVs")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   face = "italic", 
                                   color = "black"),  
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 14, color = "black"),  
        text = element_text(size = 12, color = "black"))

####Fig.1C--Faith's phylogenetic diversity######
#####1. data clean####

alpha_df <- alpha_div %>% 
  group_by(Species) %>% 
  dplyr::summarise(sd_value=sd(faith_pd),
                   faith_pd=mean(faith_pd)) %>% 
  dplyr::mutate(faith_pd=round(faith_pd,3)) %>% 
  distinct()

#####2. fig make####
p3_bar <- ggplot(alpha_df, aes(Species, faith_pd)) +
  geom_col(aes(fill = Species), width = 0.6) +  
  scale_fill_manual(values = c("#9ECAE1", "#CAB2D6FF", "#7E6148B1")) +
  geom_errorbar(data=alpha_df,mapping=aes(x = Species,
                                          ymin = faith_pd-sd_value, 
                                          ymax = faith_pd+sd_value), # error line
                width = 0.1, 
                color = 'black', 
                size=0.8)+ 
  geom_text(aes(y = faith_pd + 5, label = faith_pd), vjust = 0, 
            size = 4, color = "black") + 
  theme_classic() + 
  ylab("Faith Phylogenetic Diversity")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   face = "italic", color = "black"),  
        axis.text.y = element_text(size = 12, color = "black"),  
        axis.title = element_text(size = 14, color = "black"),  
        text = element_text(size = 12, color = "black")) 


####Fig.1D--disimilarity######
#####1. data clean####
otu <- t(raw_ASV)

otu.distance <- vegan::vegdist(otu,method="bray")

dis <-as.matrix(otu.distance)

env1 <- subset(group, Species == "H. tianxing")$ID
dis_env1 <- dis[env1,env1] %>% 
  as.dist() %>% as.vector()

env2 <- subset(group, Species == "N. concolor" )$ID
dis_env2 <- dis[env2,env2] %>% 
  as.dist() %>% as.vector()

env3 <- subset(group, Species == "N. hainanus" )$ID
dis_env3 <- dis[env3,env3] %>% 
  as.dist() %>% as.vector()

#####2. data for figure making####

dis_df <- data.frame(
  dis = c(dis_env1, dis_env2, dis_env3),  
  Species = factor(c(                 
    rep('H. tianxing', length(dis_env1)), 
    rep('N. concolor', length(dis_env2)), 
    rep('N. hainanus', length(dis_env3))
  ))
) %>% 
  dplyr::group_by(Species) %>% 
  dplyr::summarise(dis_sd=sd(dis,na.rm = TRUE),
                   dis=mean(dis)
  ) %>% 
  dplyr::mutate(dis=round(dis,3)) %>% 
  distinct()

#####3. fig make####
p4_bar <- ggplot(dis_df, aes(Species, dis)) +
  geom_col(aes(fill = Species), width = 0.6) +  
  scale_fill_manual(values = c("#9ECAE1", "#CAB2D6FF", "#7E6148B1")) +
  geom_errorbar(data=dis_df,mapping=aes(x = Species,
                                        ymin = dis-dis_sd, 
                                        ymax = dis+dis_sd), 
                width = 0.1, 
                color = 'black', 
                size=0.8)+ 
  geom_text(aes(y = dis + 0.3, label = dis), vjust = 0, size = 4, color = "black") + 
  theme_classic() +  
  ylab("Dissimilarity")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, 
                                   face = "italic", color = "black"),  
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 14, color = "black"), 
        text = element_text(size = 12, color = "black"),
        legend.position = 'none') 
  


####Fig.1E--NMDS######
#####1. data clean####

df_nmds <- metaMDS(otu.distance, k = 2)

#extract stress value（should be <=0.2）
label <- paste('Stress=',round(df_nmds$stress, 3))

#extract data for making figure
df_nmds <- as.data.frame(df_nmds$points) %>% 
  dplyr::mutate(ID=row.names(.)) %>% 
  dplyr::rename(NMDS1=colnames(.)[1],
                NMDS2=colnames(.)[2]) %>% 
  left_join(.,group)

#####2. fig make####
unique(df_nmds$Species)

p5_nmds <-ggplot(data=df_nmds,aes(x=NMDS1,y=NMDS2))+
  geom_point(aes(color = Species), shape = 19, size=3)+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  stat_ellipse(data=df_nmds,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Species),
               alpha=0.2)+
  annotate("text", x = -0.5, y = -1, 
           label = label,
           col = "black", size = 5)+
  theme_bw()+
  scale_color_manual(values =c("#9ECAE1", "#CAB2D6FF", "#7E6148B1")) +
  scale_fill_manual(values =c("#9ECAE1", "#CAB2D6FF", "#7E6148B1"))+
  theme(axis.text.x=element_text(colour="black",family="Times",size = 20,face="plain"), 
        axis.text.y=element_text(colour="black",family="Times",size=20, hjust = 1),   
        axis.title.y=element_text(colour="black",family="Times",size = 22,face="plain"),   
        axis.title.x=element_text(colour="black",family="Times",size = 22,face="plain"), 
        panel.background = element_blank(),
        legend.title = element_text(size = 15,face = "bold",
                                    vjust = 1, hjust = 0), 
        legend.text = element_text(size = 12),
        legend.key.size=unit(.5,'cm'),
        strip.text = element_text(size = 15),
        strip.text.x = element_text(size = 15, colour = "black",
                                    face = "bold.italic"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  ylim(-1,1)



####Fig.1F--community composition at phylum level######
#####1. data clean####

df_phylum <- raw_phylum %>%
  mutate(row_all = rowSums(select(., -phylum))) %>%  
  dplyr::arrange(desc(row_all)) %>% 
  dplyr::slice(1:10) %>% 
  select(-row_all) %>% # 
  tidyr::pivot_longer(
    cols = -phylum,  
    names_to = "ID", 
    values_to = "value") %>%
  left_join(.,group,by="ID") %>%
  na.exclude(phylum) %>%
  group_by(Species,month, phylum)%>%
  dplyr::summarise(cal_value=mean(value)) %>% 
  tidyr::unite("Species_month",Species:month)%>%
  plyr::ddply(.,'Species_month',transform,percent_value=cal_value/sum(cal_value)*100) %>%
  tidyr::separate(col=Species_month,
                  into = c("s1","s2","month"))%>%
  tidyr::unite("Species",s1:s2)%>%
  dplyr::mutate(phylum=factor(phylum,levels = c( "Firmicutes", "Bacteroidota", 
                                                 "Fibrobacterota","Actinobacteriota", 
                                                 "Campilobacterota", "Cyanobacteria",
                                                 "Parabasalia", "Proteobacteria",  
                                                  "Spirochaetota","Thermoplasmatota")),
                month=factor(month,levels=c("1","2","3","4","5","6",
                                            "7","8","9","10","11","12")))

#####2. fig make####
p6_phylum <- ggplot(df_phylum,aes(x=month,y=percent_value,fill=phylum))+
  geom_bar(position = "stack",width=0.8,stat="identity",colour='black')+
  facet_grid(~Species,scales= "free" ,space= "free") +
  theme_bw()+
  scale_fill_jco()+
  labs(y="Percentage (%)",fill="Phylum")+
  theme(axis.title.x=element_blank (),
        axis.text.x=element_text(colour="black",family="Times",size=14, hjust = 1), 
        axis.text.y=element_text(colour="black",family="Times",size=14,face="plain", hjust = 1), 
        axis.title.y=element_text(colour="black",family="Times",size = 20,face="plain"), 
        panel.background = element_rect(color = '#636363', 
                                        fill = 'transparent'),
        legend.title = element_text(size = 15,face = "bold",
                                    vjust = 1, hjust = 0), 
        legend.text = element_text(size = 12),
        legend.key.size=unit(.5,'cm'),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        strip.text.x = element_text(size = 12, colour = "black",
                                    face = "bold.italic"),
        strip.background = element_rect(color="black", 
                                        fill="grey95",
                                        size=0.5, 
                                        linetype="solid"))


####Fig.1G--community composition at family level######
#####1. data clean####
df_family <- raw_family %>%
  mutate(row_all = rowSums(select(., -family))) %>%  
  dplyr::arrange(desc(row_all)) %>% 
  dplyr::slice(1:10) %>% 
  select(-row_all) %>% 
  tidyr::pivot_longer(
    cols = -family, 
    names_to = "ID",  
    values_to = "value" ) %>%
  left_join(.,group,by="ID") %>%
  na.exclude(family) %>%
  group_by(Species,month, family)%>%
  dplyr::summarise(cal_value=mean(value)) %>% 
  tidyr::unite("Species_month",Species:month)%>%
  plyr::ddply(.,'Species_month',transform,percent_value=cal_value/sum(cal_value)*100) %>%
  tidyr::separate(col=Species_month,
                  into = c("s1","s2","month"))%>%
  tidyr::unite("Species",s1:s2)%>%
  dplyr::mutate(family=factor(family,levels = c( "Prevotellaceae", "Acholeplasmataceae", 
                                                 "Lachnospiraceae","Fibrobacteraceae",
                                                 "Erysipelatoclostridiaceae", "Ruminococcaceae",
                                                 "Chloroplast",  "Oscillospiraceae",  
                                                 "Atopobiaceae","Spirochaetaceae")),
                month=factor(month,levels=c("1","2","3","4","5","6",
                                            "7","8","9","10","11","12")))

#####2. fig make####

colors <- c("#A6CEE3FF","#1F78B4FF", "#B2DF8AFF", "#33A02CFF",
            "#FB9A99FF", "#E31A1CFF", "#FDBF6FFF", "#FF7F00FF",
            "#CAB2D6FF", "#6A3D9AFF", "#7E6148B1")

p7_family <- ggplot(df_family,aes(x=month,y=percent_value,fill=family))+
  geom_bar(position = "stack",width=0.8,stat="identity",colour='black')+
  facet_grid(~Species,scales= "free" ,space= "free") +
  theme_bw()+
  scale_fill_manual(values=colors)+
  labs(y="Percentage (%)",fill="Family")+
  theme(axis.title.x=element_blank (),
        axis.text.x=element_text(colour="black",family="Times",size=14, hjust = 1), 
        axis.text.y=element_text(colour="black",family="Times",size=14,face="plain", hjust = 1), 
        axis.title.y=element_text(colour="black",family="Times",size = 20,face="plain"), 
        panel.background = element_rect(color = '#636363', 
                                        fill = 'transparent'),
        legend.title = element_text(size = 15,face = "bold",
                                    vjust = 1, hjust = 0), 
        legend.text = element_text(size = 12),
        legend.key.size=unit(.5,'cm'),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        strip.text.x = element_text(size = 12, colour = "black",
                                    face = "bold.italic"),
        strip.background = element_rect(color="black", 
                                        fill="grey95",
                                        size=0.5, 
                                        linetype="solid"))

####The end######
