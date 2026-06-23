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


setwd("D:/UNL/Research/outgroup/Latest_4120/raxml/")


# Load tree
tree_file <- "R.solani.ITS.4120a.seqs.mafft.fftns2.gblocks.tree.raxml.support"
solani_tree <- read.tree(tree_file)
length(solani_tree$tip.label)

# Load metadata
metadata <- read.csv("../../../Table S1_Rsolani_4119_seqs_metadata.csv")
length(metadata$No)


# Filter metadata
metadata.filtered <- metadata[metadata$accessions %in% solani_tree$tip.label,] #Extract 
length(metadata.filtered$accessions)

# Root tree
#Root tree using the ape package; resolve.root: a logical specifying whether to resolve the new root as a bifurcating node.
new_solani_tree <- root(solani_tree, outgroup = "AY684917", resolve.root = TRUE) #athelia:AY684917

# To check if the tree is rooted; should be TRUE
is.rooted(new_solani_tree)

# Sort metadata
metadata.filtered.sorted <- metadata.filtered[match(new_solani_tree$tip.label, metadata.filtered$accessions),]

# Add continent information
metadata.filtered.sorted$continent <- countrycode(sourcevar = metadata.filtered.sorted[, "country"],
                                                  origin = "country.name",
                                                  destination = "continent")

# Add region
metadata.filtered.sorted$region <- countrycode(sourcevar = metadata.filtered.sorted[, "country"], origin = "country.name",destination = "region")

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

AG_BI <- subset(metadata.filtered.sorted, strain.beast== "AG-BI")

AG11 <- subset(metadata.filtered.sorted, strain.beast== "AG11")

AG10 <- subset(metadata.filtered.sorted, strain.beast== "AG10")

# Extract outgroup
outgroup <- subset(metadata.filtered.sorted, strain.beast== "Outgroup")

# Extract NAs
NAs <- metadata.filtered.sorted[is.na(metadata.filtered.sorted$strain.beast),]

# Add column with their respective category
#Note: the add_column function is from tiddyverse
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
AG3.focus.metadata <- rbind(AG3, AG2, AG4, AG1, AG5, AG11, AG10, AG_BI, outgroup)#, NAs)
length(AG3.focus.metadata$No)
head(AG3.focus.metadata)
#Match accessions to tree tips again because data has been moved around
AG3.focus.metadata.sorted <- AG3.focus.metadata[match(new_solani_tree$tip.label, AG3.focus.metadata$accessions),]

#Prepare data for tree
AG3.focus.metadata.sorted.tree  <- AG3.focus.metadata.sorted %>% dplyr::select(c("accessions", "AG3_focus"))
AG3.focus.metadata.sorted.tree <- aggregate(.~AG3_focus, AG3.focus.metadata.sorted.tree, FUN=paste, collapse=",")
AG3_focus <- lapply(AG3.focus.metadata.sorted.tree$accessions, function(x){unlist(strsplit(x,split=","))})
names(AG3_focus) <- AG3.focus.metadata.sorted.tree$AG3_focus

tree1 <- groupOTU(new_solani_tree, AG3_focus, "AG3_focus")

# Plot phylogeny
pdf("R.solani.ITS.4120a.seqs.mafft.fftns2.gblocks.raxml.tree.pdf", height = 5, width = 7)
tree <- ggtree(tree1, layout="fan", open.angle=15, branch.length = 'none', size=0.2, aes(colour=AG3_focus), ladderize = T) +
  geom_treescale(linesize = 0.5, fontsize = 3, color = "black") +
  #theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  #scale_color_manual(name="Reported AGs (simplified)", values=c("AG1"="orange","AG2"="blue", "AG3"="red", "AG4"="green", "FELCC"= "purple", "other.AGs"="black", "outgroup"= "cyan", "NA"="grey"))+
  scale_color_manual(name="Reported AGs (simplified)", values=c("AG1"="orange","AG2"="blue", "AG3"="red", "AG4"="green" , "outgroup"= "cyan", "AG11" ="#65350F", "AG5"="#bcbd22", "AG_BI"="#e377c2", "AG10"="purple"))+ #A599CA
  theme(legend.text = element_text(size = 12), title = element_text(size = 12), legend.title = element_text(size=12, face="bold")) +
  geom_nodelab(color="black", geom='text', size=2.5, aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50)) #+
  #geom_tiplab(size= 0.5, align=TRUE, color = "black")
tree
dev.off()

# Add the host family layer
unique(AG3.focus.metadata.sorted$host_group_summarized)
Host <- select(AG3.focus.metadata.sorted, host_group_summarized)
rownames(Host) <- new_solani_tree$tip.label

# Rename the column
colnames(Host)
names(Host)[names(Host) == "host_group_summarized"] <- "Host"

pdf("R.solani.subsample3.500a.seqs.mafft.gins.host.pdf", height = 5, width = 7)
tree2 <- gheatmap(tree, Host, color = FALSE, offset=1, width=.05,
                  colnames_angle=270, colnames_offset_x = .5, colnames_offset_y = 0.05)+
  scale_fill_viridis_d(option="A", name="Host", alpha=1) #+
  #geom_treescale(x=65, y=0, width=2, fontsize=5, linesize=1, offset=2)
tree2
dev.off()

#Plot Solanaceae
#Host_solanaceae <- select(AG3.focus.metadata.sorted, host_solanaceae)
#rownames(Host_solanaceae) <- new_solani_tree$tip.label

#tree2 <- gheatmap(tree, Host_solanaceae, color = FALSE, offset=25, width=.05,
#                  colnames_angle=270, colnames_offset_x = 9, colnames_offset_y = 0.05)+
#  scale_fill_viridis_d(option="C", name="Host_solanaceae", alpha=1)
#tree2

# Add the year layer
unique(AG3.focus.metadata.sorted$year)
Year <- select(AG3.focus.metadata.sorted, year)
rownames(Year) <- new_solani_tree$tip.label

year_ranges <- cut(Year$year,
                   breaks = c(-Inf, 1949, 1969, 1989, 1999, 2009, 2019, Inf),
                   labels = c("<= 1949", "1950-1969", "1970-1989", "1990-1999", "2000-2009", "2010-2019", "2020-Present"))

# Add year range to data frame
Year$year <- year_ranges
unique(Year$year)

tree3 <- tree2 + ggnewscale::new_scale_fill()

#tree3 <- gheatmap(tree3, Year, color = FALSE, offset=30, width=.05,
 #                 colnames_angle=270) +
  #scale_fill_viridis_c(option="G", name="Year", limits = c(1980, 2023), alpha=1, direction = 1)
#tree3

tree3 <- gheatmap(tree3, Year, color = FALSE, offset=5, width=.05,
                  colnames_angle=270) +
  scale_fill_manual(values = c("2020-Present" = "#003147", "2010-2019" = "#045275", "2000-2009" = "#3C93C2", "1990-1999" = "#6CB0D6", "1970-1989" = "#9EC9E2") , name = "Year Range")#, "1950-1969" = "#C5E1EF"), name = "Year Range")
tree3

# Add the region layer
Continent <- select(AG3.focus.metadata.sorted, continent)
rownames(Continent) <- new_solani_tree$tip.label
unique(AG3.focus.metadata.sorted$continent)
colnames(Continent)
names(Continent)[names(Continent) == "Continent"] <- "Continent"

pdf("R.solani.FELCC.subsample3.500a.seqs.mafft.gins.full.pink.pdf", height = 5, width = 7)
tree4 <- tree3 + ggnewscale::new_scale_fill()

#tree4 <- gheatmap(tree4, Region, color = FALSE, offset=60, width=.05,
 #                 colnames_angle=270) +
  #scale_fill_viridis_d(option="F", name="Region", alpha=1, direction = -1) #+
  #geom_treescale(x=65, y=0, width=2, fontsize=5, linesize=1, offset=2)
#tree4
tree4 <- gheatmap(tree4, Continent, color = FALSE, offset=9, width=.05,
                  colnames_angle=270) +
  scale_fill_manual(values = c("Africa" = "#8F003B", "Asia" = "#C40F5B" , "Europe" = "#E95694", "North America" = "#ED85B0", "Oceania" = "#F2ACCA", "South America" = "#F9D8E6"), name = "Continent")
tree4
dev.off()


tree4 <- gheatmap(tree4, Region, color = FALSE, offset=10, width=.05,
                  colnames_angle=270) +
  scale_fill_manual(values = c("East Asia & Pacific" = "#06592A", "Europe & Central Asia" = "#228B3B", "Latin America & Caribbean" = "#40AD5A", "Middle East & North Africa" = "#6FA253", "North America" = "#6CBA7D", "South Asia" = "#9CCEA7", "Sub-Saharan Africa" = "#CDE5D2"), name = "Region")
tree4

# Add the FELCC layer
#FELCC <- select(AG3.focus.metadata.sorted, FELCC)
#rownames(FELCC) <- new_solani_tree$tip.label

# Rename the column
#colnames(FELCC)
#names(FELCC)[names() == "host_group_summarized"] <- "Host_Family"

#tree5 <- tree4 + ggnewscale::new_scale_fill()
#tree5 <- gheatmap(tree, FELCC, color = FALSE, offset=0.1, width=.05,
 #                 colnames_angle=270, colnames_offset_x = 9, colnames_offset_y = 0.05)+
  #scale_fill_viridis_d(option="C", name="FELCC", alpha=1)
#tree5

####### Preparing the population datasets ######
# Extract specific clade (AG4)
# Find the MRCA node for the target clade 
# Replace 'tip_label1' and 'tip_label2' with the labels of the tips defining the clade
mrca_node <- MRCA(tree1, "MN160704", "KF686783")

#Extract the clade
AG4_clade_tree <- extract.clade(tree1, mrca_node)

AG4_clade_tips <- AG4_clade_tree$tip.label

#Get all tip labels (accessions) in the extracted clade
clade_tip_labels <- data.frame(AG4_accession = AG4_clade_tree$tip.label)

#Print the tip labels
print(clade_tip_labels)

#Convert the tip labels to a data frame with a column header "AG3_accession"
clade_tip_labels_df <- data.frame(AG1_accession = clade_tip_labels)

# Write the data frame to a CSV file
write.csv(clade_tip_labels, file = "AG4.1174.accessions.csv", row.names = FALSE)
