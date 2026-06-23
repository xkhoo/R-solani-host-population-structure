#install.packages("forecast")
require(ggplot2)
require(ggtree)
library(treeio)
library(tidyverse)
library(countrycode)
library(rBt)
library(dplyr)
library(igraph)
library(ggraph)
library(ggforce)
library(cowplot)
library(phytools)
library(forecast)

# Host-transition network analysis for Rhizoctonia solani
#
# This script is the original working script used to summarize host-transition
# counts from BEAST annotated history files and generate host-transition network
# visualizations.
#
# Large input files are not stored in this GitHub repository. To run this script,
# download/generate the processed BEAST output files


################### Empirical model (observed host) ############################
solani_beast1 <- read.annot.beast("R.solani.subsample1.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

is.rooted(solani_beast1)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast1$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast1$tip.label, metadata.filtered$accessions),]

#new_solani_beast1$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
 # add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast1$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree1 <- groupOTU(solani_beast1, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data1 <- fortify(tree1$metadata)

solani_beast2 <- read.annot.beast("R.solani.subsample2.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast2 <- root(solani_beast2, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast2)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast2$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast2$tip.label, metadata.filtered$accessions),]

#new_solani_beast2$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
 # add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast2$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree2 <- groupOTU(solani_beast2, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data2 <- fortify(tree2$metadata)

solani_beast3 <- read.annot.beast("R.solani.subsample3.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast3 <- root(solani_beast3, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast3)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast3$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast3$tip.label, metadata.filtered$accessions),]

#new_solani_beast3$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
 # add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast3$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree3 <- groupOTU(solani_beast3, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data3 <- fortify(tree3$metadata)

solani_beast4 <- read.annot.beast("R.solani.subsample4.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast4 <- root(solani_beast4, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast4)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast4$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast4$tip.label, metadata.filtered$accessions),]

#new_solani_beast4$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast4$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree4 <- groupOTU(solani_beast4, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data4 <- fortify(tree4$metadata)

solani_beast5 <- read.annot.beast("R.solani.subsample5.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast5 <- root(solani_beast5, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast5)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast5$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast5$tip.label, metadata.filtered$accessions),]

#new_solani_beast5$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
 # add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast5$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree5 <- groupOTU(solani_beast5, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data5 <- fortify(tree5$metadata)

solani_beast6 <- read.annot.beast("R.solani.subsample6.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast6 <- root(solani_beast6, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast6)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast6$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast6$tip.label, metadata.filtered$accessions),]

#new_solani_beast6$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
#  add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast6$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree6 <- groupOTU(solani_beast6, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data6 <- fortify(tree6$metadata)

solani_beast7 <- read.annot.beast("R.solani.subsample7.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast7 <- root(solani_beast7, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast7)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast7$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast7$tip.label, metadata.filtered$accessions),]

#new_solani_beast7$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
#  add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast7$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree7 <- groupOTU(solani_beast7, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data7 <- fortify(tree7$metadata)

solani_beast8 <- read.annot.beast("R.solani.subsample8.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast8 <- root(solani_beast8, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast8)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast8$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast8$tip.label, metadata.filtered$accessions),]

#new_solani_beast8$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
#  add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast8$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree8 <- groupOTU(solani_beast8, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data8 <- fortify(tree8$metadata)

solani_beast9 <- read.annot.beast("R.solani.subsample9.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast9 <- root(solani_beast9, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast9)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast9$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast9$tip.label, metadata.filtered$accessions),]

#new_solani_beast9$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
#  add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast9$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree9 <- groupOTU(solani_beast9, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data9 <- fortify(tree9$metadata)

solani_beast10 <- read.annot.beast("R.solani.subsample10.AG.host.200a.seqs.mafft.gins.chain2.host_groups.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_solani_beast10 <- root(solani_beast10, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(solani_beast10)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_beast10$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(solani_beast10$tip.label, metadata.filtered$accessions),]

#new_solani_beast10$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
#  add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(solani_beast10$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


tree10 <- groupOTU(solani_beast10, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
tree_data10 <- fortify(tree10$metadata)

#tree_data1
tree_data_expanded1 <- tree_data1 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches1 <- tree_data_expanded1 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded2 <- tree_data2 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches2 <- tree_data_expanded2 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded3 <- tree_data3 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches3 <- tree_data_expanded3 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded4 <- tree_data4 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches4 <- tree_data_expanded4 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded5 <- tree_data5 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches5 <- tree_data_expanded5 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded6 <- tree_data6 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches6 <- tree_data_expanded6 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded7 <- tree_data7 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches7 <- tree_data_expanded7 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded8 <- tree_data8 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches8 <- tree_data_expanded8 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded9 <- tree_data9 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches9 <- tree_data_expanded9 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

tree_data_expanded10 <- tree_data10 %>%
  separate_rows(host_groups.states.set, sep = ",") %>%
  mutate(host_groups.states.set = trimws(host_groups.states.set)) 
host_switches10 <- tree_data_expanded10 %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')


all_host_switches <- list(host_switches1, host_switches2, host_switches3, host_switches4, host_switches5, host_switches6, host_switches7, host_switches8, host_switches9, host_switches10) # Add more as needed

#Combine all the data frames into one
combined_switches <- bind_rows(all_host_switches) %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::summarize(average_jumps = mean(count), .groups = 'drop')

edges_avg <- combined_switches %>%
  filter(host_groups.states != host_groups.states.set) %>%
  dplyr::select(from = host_groups.states, to = host_groups.states.set, empirical_weight = average_jumps)

#For jumps per tree
combined_switches_sf <- bind_rows(all_host_switches) %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::select(from = host_groups.states, to = host_groups.states.set, empirical_weight = count)

empirical_summary <- combined_switches_sf %>%
  filter(from != to) %>%
  group_by(from, to) %>%
  dplyr::summarise(
    median_jumps = median(empirical_weight),
    lower_whisker = quantile(empirical_weight, 0.25),
    upper_whisker = quantile(empirical_weight, 0.75),
    mean_jumps = mean(empirical_weight)
  )

host_jumps_per_tree <- ggplot(empirical_summary, aes(x = from, fill = to)) +
                        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.9, aes(y = mean_jumps)) +
                        #geom_point(size = 3, position = position_jitter(width = 0.2)) +
                        geom_errorbar(aes(ymin=lower_whisker, ymax=upper_whisker), width=0.2,position=position_dodge(width=0.9)) +
                        #aes(ymin = lower_whisker, ymax = upper_whisker), 
                         #           width = 0.2, position = position_dodge(width = 0.7)) +
                        geom_point(aes(y = median_jumps),  
                                    size = 2, shape = 21,
                                    color = "black",
                                    position = position_dodge(width = 0.9)) +
                        scale_color_manual(values = rainbow(5)) +  # Use different colors for "to" species
                        scale_shape_manual(values = c(21, 22, 16, 24, 25)) +  # Different shapes
                        labs(title = "Number of Jumps Per Tree", 
                              y = "Number of Jumps Per Tree", 
                              x = "From Species", 
                              fill = "To Species") +
                        theme_minimal() +
                        theme(
                        #panel.border = element_rect(color = "black", fill = NA, size = 1), 
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(), 
                        axis.line.y = element_line(color = "black"),
                        axis.ticks.length.y = unit(0.2, "cm"),
                        axis.ticks.y.left = element_line(color = "black"),
                        #y.axis.line = element_line(color = "black"),
                        legend.position = "right")
host_jumps_per_tree

#edges_avg <- edges_avg %>%
#  mutate(Z = ifelse(is.na(Z), 0, Z),  # Fill NA Z-scores with 0
#        P.adj = ifelse(is.na(P.adj), 1, P.adj),  # Fill NA adjusted P-values with 1 (not significant)
#       significant = ifelse(is.na(significant), FALSE, significant))  # Mark as non-significant if NA


################### Null model (randomized host) ################
random_beast1 <- read.annot.beast("../random/R.solani.subsample1.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast1 <- root(random_beast1, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast1)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast1$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast1$tip.label, metadata.filtered$accessions),]

#new_random_beast1$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast1$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree1 <- groupOTU(random_beast1, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data1 <- fortify(random_tree1$metadata)

random_beast2 <- read.annot.beast("../random/R.solani.subsample2.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast2 <- root(random_beast2, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast2)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast2$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast2$tip.label, metadata.filtered$accessions),]

#new_random_beast2$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast2$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree2 <- groupOTU(random_beast2, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data2 <- fortify(random_tree2$metadata)

random_beast3 <- read.annot.beast("../random/R.solani.subsample3.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast3 <- root(random_beast3, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast3)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast3$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast3$tip.label, metadata.filtered$accessions),]

#new_random_beast3$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast3$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree3 <- groupOTU(random_beast3, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data3 <- fortify(random_tree3$metadata)

random_beast4 <- read.annot.beast("../random/R.solani.subsample4.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast4 <- root(random_beast4, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast4)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast4$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast4$tip.label, metadata.filtered$accessions),]

#new_random_beast4$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast4$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree4 <- groupOTU(random_beast4, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data4 <- fortify(random_tree4$metadata)

random_beast5 <- read.annot.beast("../random/R.solani.subsample5.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast5 <- root(random_beast5, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast5)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast5$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast5$tip.label, metadata.filtered$accessions),]

#new_random_beast5$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast5$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree5 <- groupOTU(random_beast5, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data5 <- fortify(random_tree5$metadata)

random_beast6 <- read.annot.beast("../random/R.solani.subsample6.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast6 <- root(random_beast6, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast6)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast6$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast6$tip.label, metadata.filtered$accessions),]

#new_random_beast6$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast6$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree6 <- groupOTU(random_beast6, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data6 <- fortify(random_tree6$metadata)

random_beast7 <- read.annot.beast("../random/R.solani.subsample7.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast7 <- root(random_beast7, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast7)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast7$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast7$tip.label, metadata.filtered$accessions),]

#new_random_beast7$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast7$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree7 <- groupOTU(random_beast7, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data7 <- fortify(random_tree7$metadata)

random_beast8 <- read.annot.beast("../random/R.solani.subsample8.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast8 <- root(random_beast8, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast8)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast8$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast8$tip.label, metadata.filtered$accessions),]

#new_random_beast8$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast8$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree8 <- groupOTU(random_beast8, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data8 <- fortify(random_tree8$metadata)

random_beast9 <- read.annot.beast("../random/R.solani.subsample9.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast9 <- root(random_beast9, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast9)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast9$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast9$tip.label, metadata.filtered$accessions),]

#new_random_beast9$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


#Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast9$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree9 <- groupOTU(random_beast9, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data9 <- fortify(random_tree9$metadata)

random_beast10 <- read.annot.beast("../random/R.solani.subsample10.AG.host.200a.seqs.mafft.gins.randomized.host.history.25burnin.combined.output.txt")

# Load metadata
metadata <- read.csv("table/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")
length(metadata$No)

#new_random_beast10 <- root(random_beast10, outgroup = "AY684917", resolve.root = TRUE) #Ceratobasidium:PP713314, athelia:GU080230;AY684917, waitea:DQ356417
#To check if the tree is rooted; should be TRUE
is.rooted(random_beast10)

#Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% random_beast10$tip.label,] #Extract 
length(metadata.filtered$accessions)

#Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(random_beast10$tip.label, metadata.filtered$accessions),]

#new_random_beast10$tip.metadata <- metadata.filtered.sorted$strain.beast
#Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

#Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

#Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

#Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

#Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

#outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

AG3  <- AG3  %>%
  add_column(AG3_focus = "AG3")

AG2  <- AG2  %>%
  add_column(AG3_focus = "AG2")

AG4  <- AG4  %>%
  add_column(AG3_focus = "AG4")

AG1 <- AG1 %>%
  add_column(AG3_focus = "AG1")

AG5  <- AG5  %>%
  add_column(AG3_focus = "AG5")

AG11  <- AG11  %>%
  add_column(AG3_focus = "AG11")

AG_BI  <- AG_BI  %>%
  add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

#outgroup <- outgroup  %>%
# add_column(AG3_focus = "outgroup")


# Bind datasets
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG10, AG11, AG_BI)
length(AG3.focus.metadata$No)

# Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(random_beast10$tip.label, AG3.focus.metadata$accessions),]

# Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus


random_tree10 <- groupOTU(random_beast10, AG3_focus, "AG3_focus")

#Extract metadata associated with the tree
random_tree_data10 <- fortify(random_tree10$metadata)

random_tree_data_expanded1 <- random_tree_data1 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches1 <- random_tree_data_expanded1 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded2 <- random_tree_data2 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches2 <- random_tree_data_expanded2 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded3 <- random_tree_data3 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches3 <- random_tree_data_expanded3 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded4 <- random_tree_data4 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches4 <- random_tree_data_expanded4 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded5 <- random_tree_data5 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches5 <- random_tree_data_expanded5 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded6 <- random_tree_data6 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches6 <- random_tree_data_expanded6 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded7 <- random_tree_data7 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches7 <- random_tree_data_expanded7 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded8 <- random_tree_data8 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches8 <- random_tree_data_expanded8 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded9 <- random_tree_data9 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches9 <- random_tree_data_expanded9 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

random_tree_data_expanded10 <- random_tree_data10 %>%
  separate_rows(host.states.set, sep = ",") %>%
  mutate(host.states.set = trimws(host.states.set)) 
random_host_switches10 <- random_tree_data_expanded10 %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(count = n(), .groups = 'drop')


all_random_host_switches <- list(random_host_switches1, random_host_switches2, random_host_switches3, random_host_switches4, random_host_switches5, random_host_switches6, random_host_switches7, random_host_switches8, random_host_switches9, random_host_switches10) # Add more as needed

# Combine all the data frames into one using bind_rows
random_combined_switches <- bind_rows(all_random_host_switches) %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::summarize(average_jumps = mean(count), .groups = 'drop')

random_edges_avg <- random_combined_switches %>%
  filter(host.states != host.states.set) %>%
  dplyr::select(from = host.states, to = host.states.set, null_weight = average_jumps)


#Combine empirical and null models
combined_edges_avg <- left_join(edges_avg, random_edges_avg, by = c("from", "to"))

sd_null_weight <- sd(combined_edges_avg$null_weight)

#Compute the z-scores
combined_edges_avg <- combined_edges_avg %>%
  mutate(z_score = (empirical_weight - null_weight) / sd_null_weight) %>%
  mutate(Significance = ifelse(abs(z_score) >= 1.96, TRUE, FALSE)) 

print(combined_edges_avg)

combined_edges_avg <- combined_edges_avg %>%
  filter(Significance ==TRUE)

edges_aligned_avg <- combined_edges_avg %>%
  dplyr::mutate(count_class = ifelse(empirical_weight >= 20, "high", "low"))
edges_aligned_avg <- edges_aligned_avg %>%
  dplyr::mutate(arrow_size = ifelse(count_class == "high", 5, 2))


#Create the graph object from the filtered aligned data
g_avg <- graph_from_data_frame(d = edges_aligned_avg, directed = TRUE)


#Plot the network with the filtered data
p_avg <- ggraph(g_avg, layout = "circle") +
  geom_node_tile(aes(fill = name), width = 0.5, height = 1) +
  geom_edge_parallel(aes( edge_width = arrow_size, color = z_score), 
                     arrow = arrow(length = unit(edges_aligned_avg$arrow_size, 'mm'), 
                                   type = "closed", ends = "last"),
                     show.legend = TRUE,
                     start_cap = circle(10, "mm"),
                     end_cap = circle(12, "mm"), 
                     lineend = "square",
                     linejoin = "round",
                     position = "identity",
                     n=100, angle_calc = "across", sep = unit(1.5, 'cm'), linemitre = 2) +
  #geom_node_point(aes(color = name), size = 8) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, color = "black") +
  #scale_fill_viridis_d(name = "Node Legend", option = "C") +
  theme_void() +
  scale_edge_color_gradient(high = "orange", low = "blue") +#, #limits = c(-7, 2)) +
  scale_edge_width(range = c(0.5, 5)) +
  guides(edge_width = guide_legend("Average Host Jumps"), 
         edge_color = guide_edge_colorbar("Significance (z-score)"),
         edge_alpha = guide_legend("Z-Score")) +
  labs(title = "Host-Switching Network", subtitle = "Frequency and Significance of Host Jumps")
p_avg


str(edges_aligned_avg$empirical_weight)
str(edges_aligned_avg$z_score)
str(edges_aligned_avg$arrow_size)

combined_switches_sf <- bind_rows(all_host_switches) %>%
  dplyr::group_by(host_groups.states, host_groups.states.set) %>%
  dplyr::select(from = host_groups.states, to = host_groups.states.set, empirical = count)

random_combined_switches_sf <- bind_rows(all_random_host_switches) %>%
  dplyr::group_by(host.states, host.states.set) %>%
  dplyr::select(from = host.states, to = host.states.set, null = count)


all_summary <- left_join(combined_switches_sf, random_combined_switches_sf, by = c("from", "to"))

all_summary <- all_summary %>%
  filter(from != to) %>%
  group_by(from, to) %>%
  dplyr::summarise(
    median_jumps = median(empirical),
    lower_whisker = quantile(empirical, 0.25),
    upper_whisker = quantile(empirical, 0.75),
    mean_jumps = mean(empirical)
  )

host_jumps_per_tree <- ggplot(empirical_summary, aes(x = from, fill = to)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.9, aes(y = mean_jumps)) +
  #geom_point(size = 3, position = position_jitter(width = 0.2)) +
  geom_errorbar(aes(ymin=lower_whisker, ymax=upper_whisker), width=0.2,position=position_dodge(width=0.9)) +
  #aes(ymin = lower_whisker, ymax = upper_whisker), 
  #           width = 0.2, position = position_dodge(width = 0.7)) +
  geom_point(aes(y = median_jumps),  
             size = 2, shape = 21,
             color = "black",
             position = position_dodge(width = 0.9)) +
  scale_color_manual(values = rainbow(5)) +  # Use different colors for "to" species
  scale_shape_manual(values = c(21, 22, 16, 24, 25)) +  # Different shapes
  labs(title = "Number of Jumps Per Tree (1950 - 2023)", 
       y = "Number of Jumps Per Tree", 
       x = "From Species", 
       point = "To Species") +
  theme_minimal() +
  theme(legend.position = "right")
host_jumps_per_tree
