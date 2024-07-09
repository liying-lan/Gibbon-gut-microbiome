####This code is to do LEfSe analysis####

####history clean####
rm(list=ls())
####please setting the working directory manually####

####load packages
library(tidyverse)
library(dplyr)
library(microeco)
library(ggtree)
library(magrittr)

####load data#####
raw_env <- read.csv("metadata.csv",
                    header = TRUE,row.names = 1, fileEncoding = 'GBK')

raw_otu <- read.csv("asv_tab.csv",
                    header = T, check.names = F,row.names=1,fill=T)

tax <- read.csv("taxonomy.csv",
                header = T, check.names = F,row.names=1,fill=T)


####1.data clean#####
###please change the species name here#####
sample_table <-  raw_env %>%
  dplyr::filter(Species=="Hoolock_tianxing") %>%
  dplyr::mutate(leaf_group=factor(leaf_group,levels=c("low","middle","high"))) %>%
  select(leaf_group,everything()) %>%
  #dplyr::rename(group=leaf_group) %>%
  na.omit(group)

feature_table <- raw_otu[,row.names(sample_table)] %>%
  dplyr::mutate(row_all=rowSums(.,na.rm =T)) %>%
  dplyr::arrange(row_all=desc(row_all)) %>%
  dplyr::slice(1:1000) %>%
  dplyr::select(-last_col()) %>%
  data.frame() 

tax %<>% as.data.frame()

####2.build dataset and lefse object#####
dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax)

lefse <- trans_diff$new(dataset = dataset, 
                         method = "lefse", 
                         group = "leaf_group", 
                         alpha = 0.01, 
                         lefse_subgroup = NULL)


####3.plot_diff_cladogram####
cols <- c("#7570B3","#D95F02","#1b9e77")

lefse$plot_diff_cladogram(use_taxa_num = 100, 
                          use_feature_num = 50, 
                          clade_label_level = 5, 
                          group_order = c("low", "middle", "high")) +
  ggplot2::scale_colour_manual(values=cols) +
  ggplot2::scale_fill_manual(values=cols) 

####4.plot_diff_bar#####

lefse$plot_diff_bar(threshold = 4,
                    use_number = 1:30, 
                    width = 0.6, 
                    group_order = c("low", "middle", "high")) +
  ggplot2::scale_colour_manual(values=cols) +
  ggplot2::scale_fill_manual(values=cols) 


