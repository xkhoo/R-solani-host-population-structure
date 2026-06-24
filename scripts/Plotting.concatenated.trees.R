library(ggtree)
library(ggplot2)
library(ape)
library(countrycode)
library(phytools)
library(ggnewscale)
library(dplyr)
library(tidyverse)
library(treeio)
library(ggrepel)



#Load tree
tree_file <- "R.solani.concatenated.66.seqs.mafft.gins.tree.raxml.support"
solani_tree <- read.tree(tree_file)
length(solani_tree$tip.label)


length(metadata$No)
### tef1
# Supplementary Table S2
metadata <- read.csv("R.solani.tef1.229.filtered.ITS.complete.final.metadata.csv")

# Filter metadata
metadata.filtered <- metadata[metadata$accession %in% solani_tree$tip.label,] #Extract 
length(metadata.filtered$accession)

# Root tree
new_solani_tree <- root(solani_tree, outgroup = "MN509438_MN071107", resolve.root = TRUE)
#tef=#MN509438_MN071107

# To check if the tree is rooted; should be TRUE
is.rooted(new_solani_tree)

# Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(new_solani_tree$tip.label, metadata.filtered$accession),]

# Add continent information
metadata.filtered.sorted$continent <- countrycode(sourcevar = metadata.filtered.sorted[, "geographic_location"],
                                                  origin = "country.name",
                                                  destination = "continent")

# Add region
metadata.filtered.sorted$region <- countrycode(sourcevar = metadata.filtered.sorted[, "geographic_location"], origin = "country.name",destination = "region")

# Make all <NA> into real NAs
metadata.filtered.sorted[metadata.filtered.sorted=="<NA>"] = NA  

# Find unique AGs
unique(metadata.filtered.sorted$strain.beast)

# Extract AG3
AG3 <- subset(metadata.filtered.sorted, strain.beast== "AG3")

# Extract AG2
AG2 <- subset(metadata.filtered.sorted, strain.beast== "AG2")

# Extract AG4
AG4 <- subset(metadata.filtered.sorted, strain.beast== "AG4")

# Extract AG1
AG1 <- subset(metadata.filtered.sorted, strain.beast== "AG1")

AG5 <- subset(metadata.filtered.sorted, strain.beast== "AG5")

AG8 <- subset(metadata.filtered.sorted, strain.beast== "AG8")

#AG_BI <- subset(metadata.filtered.sorted, strain.summarized== "AG-BI")

#AG11 <- subset(metadata.filtered.sorted, strain.summarized== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

AG6 <- subset(metadata.filtered.sorted, strain.beast== "AG6")

#Extract FELCC
#FELCC <- subset(metadata.filtered.sorted, strain.summarized== "FELCC")

# Extract outgroup
outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

# Extract NAs
NAs <- metadata.filtered.sorted[is.na(metadata.filtered.sorted$strain.beast),]

#Extract all other AGs
#other.AGs <- metadata.filtered.sorted %>% filter(strain.summarized %in% c("AG5", "AG6","AG11","AG-BI", "AG10"))


# Check if all the values add up
length(metadata.filtered.sorted$No)
# Extracted values
length(AG3$No)
length(AG4$No)
length(outgroup$No)
length(FELCC$No)
length(AG2$No)
length(NAs$No)
length(other.AGs$No)
length(AG1$No)

# Add column with their respective category
# Note: the add_column function is from tiddyverse
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

#AG11  <- AG11  %>%
 # add_column(AG3_focus = "AG11")

#AG_BI  <- AG_BI  %>%
  #add_column(AG3_focus = "AG_BI")

AG10 <- AG10 %>%
  add_column(AG3_focus = "AG10")

AG8  <- AG8  %>%
  add_column(AG3_focus = "AG8")

AG6  <- AG6  %>%
  add_column(AG3_focus = "AG6")

#FELCC <- FELCC  %>%
 # add_column(AG3_focus = "FELCC")

outgroup <- outgroup  %>%
  add_column(AG3_focus = "outgroup")

NAs <- NAs  %>%
  add_column(AG3_focus = "NA")

#other.AGs <- other.AGs  %>%
 # add_column(AG3_focus = "other.AGs")

# Bind datasets
#AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, FELCC, outgroup, NAs, other.AGs)
AG3.focus.metadata <- rbind(AG3, AG2,AG1, outgroup)
length(AG3.focus.metadata$No)
head(AG3.focus.metadata)
# Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(new_solani_tree$tip.label, AG3.focus.metadata$accession),]

# Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accession", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accession, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus

#write.csv(AG3.focus.metadata.sorted, file = "test.csv")

tree1 <- groupOTU(new_solani_tree, AG3_focus, "AG3_focus")


# Plot
pdf("R.solani.concatenated.66.seqs.mafft.gins.tree.pdf", height = 6, width = 8)
tree <- ggtree(tree1, layout="rectangular", branch.length = 'none', size=0.5, aes(colour=AG3_focus), ladderize = T) +
  geom_treescale(x=0, y=45, width=1, color='black') +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  #scale_color_manual(name="Reported AGs (simplified)", values=c("AG1"="orange","AG2"="blue", "AG3"="red", "AG4"="green", "FELCC"= "purple", "other.AGs"="black", "outgroup"= "cyan", "NA"="grey"))+
  scale_color_manual(name="Reported AGs (simplified)", values=c("AG1"="orange","AG2"="blue", "AG3"="red", "AG4"="green" , "outgroup"= "cyan", "AG8" ="#65350F", "AG5"="#bcbd22", "AG10"="#e377c2", "FELCC" = "purple","NAs"="grey"))+ #A599CA
  theme(legend.text = element_text(size = 8), title = element_text(size = 8), legend.title = element_text(size=8, face="bold")) +
  geom_nodelab(color="black", geom='text', size=2, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label))) +
  geom_tiplab(size= 1, align=TRUE)
tree
dev.off()
