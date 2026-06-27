# Load Required Libraries
library(poppr)
library(igraph)
library(countrycode)
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(gganimate)
library(dartR)
library(directlabels)
library(tidyr)


#==========================================#
# Part A: Population Structure by Continent#
#==========================================#
#===========================
# Load Genlight Object
#===========================
genlight <- fasta2genlight(file = "AG3.914.seqs.mafft.fftns2.replaced.fasta")

# Load and Process Metadata
metadata <- read.csv("Rsolani.AG3.metadata.csv") # AG-3 population-specific metadata filtered from Supplementary Table S1
metadata.filtered <- metadata[metadata$accessions %in% genlight$ind.names,]
metadata.filtered <- replace(metadata.filtered, is.na(metadata.filtered), 0)
metadata.filtered.sorted <- metadata.filtered[match(genlight$ind.names, metadata.filtered$accessions),]
strata(genlight) <- metadata.filtered.sorted

# Check the order
head(genlight$ind.names)
head(metadata.filtered.sorted$accessions)

tail(genlight$ind.names)
tail(metadata.filtered.sorted$accessions)

#===========================
# Clustering and DAPC Setup
#===========================
setPop(genlight) <- ~continent
maxK <- 10
myMat <- matrix(nrow=4, ncol=maxK)
colnames(myMat) <- 1:maxK

for (i in 1:nrow(myMat)) {
  grp <- find.clusters(genlight, n.pca = 10, choose.n.clust = FALSE, max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

# BIC Boxplot for Optimal K
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)

pdf("ag3.1079.region.bic.maxk4.pdf", height = 10, width = 12)
ggplot(my_df, aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=15, face="bold")) +
  xlab("Number of groups (K)")
dev.off()

# Run DAPC for K = 4 to 10
my_k <- 4:10
grp_l <- vector("list", length(my_k))
dapc_l <- vector("list", length(my_k))

for (i in seq_along(my_k)) {
  set.seed(9)
  grp_l[[i]] <- find.clusters(genlight, n.pca = 10, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(genlight, pop = grp_l[[i]]$grp, n.pca = 10, n.da = my_k[i])
}

#===========================
# DAPC Scatter Plots
#===========================
my_df <- as.data.frame(dapc_l[[length(dapc_l)]]$ind.coord)
my_df$Group <- dapc_l[[length(dapc_l)]]$grp
my_pal <- c(brewer.pal(8, "Dark2"), brewer.pal(10, "Paired"))

pdf("ag4.1174.continent.LD1.LD2.pdf", height = 10, width = 12)
ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 6, shape = 21) +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size = 15, face="bold")) +
  scale_color_manual(values = my_pal) +
  scale_fill_manual(values = paste0(my_pal, "66"))
dev.off()

#===========================
# Principal Component Analysis
#===========================
solani.pca <- glPca(genlight, nf = 10)
var_explained <- solani.pca$eig / sum(solani.pca$eig) * 100
cumulative_var_explained <- cumsum(var_explained)
optimal_n_pca <- min(which(cumulative_var_explained >= 80))

barplot(100 * solani.pca$eig / sum(solani.pca$eig), col = heat.colors(50), 
        main="PCA Eigenvalues", ylab="Percent of variance\nexplained", xlab="Eigenvalues")

plot(cumulative_var_explained, type = "b", xlab = "Number of Principal Components",
     ylab = "Cumulative Variance Explained", main = "Cumulative Variance Explained")
abline(h = 80, col = "red", lty = 2)
abline(h = 90, col = "blue", lty = 2)
abline(v = optimal_n_pca, col = "green", lty = 2)

solani.pca.scores <- as.data.frame(solani.pca$scores)
solani.pca.scores$pop <- pop(genlight)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(metadata.filtered.sorted$continent)))

pdf("R.solani.ag3.1079.region.pc1xpc2.pdf", height = 10, width = 12)
ggplot(solani.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  stat_ellipse(level = 0.99, linewidth = 1) +
  scale_color_manual(values = mycolors) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()
dev.off()

#===========================
# DAPC Summary and Structure
#===========================
pnw.dapc <- dapc(genlight, n.pca = 6, n.da = 4)
scatter(pnw.dapc, col = mycolors, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft",
        scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch = 19, axesell = TRUE,
        label.inds = TRUE, main = "DAPC Scatter Plot", sub = "Colored by Group",
        xlab = "Discriminant Function 1", ylab = "Discriminant Function 2", cstar = 1)

# Admixture Plot
pdf("ag3.1079.continent.admixtureplot.pdf", height = 10, width = 12)
compoplot(pnw.dapc, col = mycolors, posi = 'bottom', space = 0.8,
          main = "Admixture Plot", xlab = "Individuals")
dev.off()

#===========================
# Population Structure by Group
#===========================
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(genlight)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
colnames(dapc.results) <- c("Original_Pop", "Sample", "Assigned_Pop", "Posterior_membership_probability")

pdf("ag3.1079.continent.population.structure.pdf", height = 10, width = 12)
ggplot(dapc.results, aes(x = Sample, y = Posterior_membership_probability, fill = Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = mycolors) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2),
        strip.text.x = element_text(size = 6))
dev.off()

#=============================================
# Part B: Population Structure by Host Species
#=============================================
setPop(genlight) <- ~host_group_summarized

# Test optimal number of clusters (K)
maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:maxK
for(i in 1:nrow(myMat)){
  grp <- find.clusters(genlight, n.pca = 10, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)

pdf("ag3.1079.host.all.bic.pdf", height = 10, width = 12)
ggplot(my_df, aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size = 15, face="bold")) +
  xlab("Number of groups (K)")
dev.off()

# Run DAPC for K = 2 to 10
my_k <- 2:10
grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] <- find.clusters(genlight, n.pca = 8, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(genlight, pop = grp_l[[i]]$grp, n.pca = 8, n.da = my_k[i])
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
my_pal <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 12))

# DAPC plots LD1-LD3
pdf("Rsolani.ag4.1079.all.host.LD1.LD2.pdf", height = 10, width = 12)
ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 5, shape = 21) +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size = 15, face="bold")) +
  scale_color_manual(values=c(my_pal)) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
dev.off()

pdf("Rsolani.ag3.1079.all.host.LD2.LD3.pdf", height = 10, width = 12)
ggplot(my_df, aes(x = LD2, y = LD3, color = Group, fill = Group)) +
  geom_point(size = 5, shape = 21) +
  theme_bw() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size = 15, face="bold")) +
  scale_color_manual(values=c(my_pal)) +
  scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
dev.off()

# PCA based on host.solanaceae
setPop(genlight) <- ~host.solanaceae
solani.pca <- glPca(genlight, nf = 10)
var_explained <- solani.pca$eig / sum(solani.pca$eig) * 100
cumulative_var_explained <- cumsum(var_explained)
optimal_n_pca <- min(which(cumulative_var_explained >= 80))

barplot(100*solani.pca$eig/sum(solani.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

plot(cumulative_var_explained, type = "b", xlab = "Number of Principal Components",
     ylab = "Cumulative Variance Explained", main = "Cumulative Variance Explained by Principal Components")
abline(h = 80, col = "red", lty = 2)
abline(h = 90, col = "blue", lty = 2)
abline(v = optimal_n_pca, col = "green", lty = 2)

solani.pca.scores <- as.data.frame(solani.pca$scores)
solani.pca.scores$pop <- pop(genlight)
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

pdf("ag3.solanaceae.only.1079.all.host.pc1xpc2.pdf", height = 10, width = 12)
ggplot(solani.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2) +
  stat_ellipse(level = 0.99, linewidth = 1) +
  scale_color_manual(values = mycolors) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()
dev.off()

# PCA eigen scatter
pdf("ag3.solanaceae.only.1079.all.host.PCA.pdf", height = 10, width = 12)
gl.pcoa.plot(solani.pca, genlight, ellipse=TRUE, plevel=0.95,hadjust=1.5, vadjust=1, xaxis=1, yaxis=3)
if (requireNamespace("plotly", quietly = TRUE)) {
  gl.pcoa.plot(solani.pca, genlight, interactive=T)
}
dev.off()

#===========================
# Cluster Stability Assessment (DAPC Randomization)
#===========================

# Repeat DAPC multiple times with different random seeds to assess cluster stability
repeat_dapc <- function(genlight_obj, n.pca = 8, n.da = 3, n.clust = 4, reps = 50) {
  cluster_assignments <- matrix(NA, nrow = nInd(genlight_obj), ncol = reps)
  rownames(cluster_assignments) <- indNames(genlight_obj)
  
  for (i in 1:reps) {
    set.seed(i)
    grp <- find.clusters(genlight_obj, n.pca = n.pca, n.clust = n.clust)
    dapc_res <- dapc(genlight_obj, pop = grp$grp, n.pca = n.pca, n.da = n.da)
    cluster_assignments[, i] <- as.character(dapc_res$assign)
  }
  
  return(cluster_assignments)
}

dapc_replicates <- repeat_dapc(genlight, n.pca = 8, n.da = 3, n.clust = 3, reps = 50)

assignment_summary <- t(apply(dapc_replicates, 1, function(x) table(factor(x, levels = unique(as.vector(dapc_replicates))))))
head(assignment_summary)
assignment_proportions <- assignment_summary / rowSums(assignment_summary)
head(assignment_proportions)
library(pheatmap)
pheatmap(assignment_proportions,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Cluster Assignment Stability",
         fontsize_row = 0.5)
dominant_cluster <- apply(assignment_proportions, 1, which.max)
assignment_proportions <- assignment_proportions[order(dominant_cluster), ]
majority_cluster <- apply(dapc_replicates, 1, function(x) names(sort(table(x), decreasing = TRUE))[1])
write.csv(data.frame(Sample = rownames(dapc_replicates), MajorityCluster = majority_cluster),
          "consensus_clusters.csv", row.names = FALSE)


# Reorder the vector to match genlight individual order
majority_cluster <- majority_cluster[match(indNames(genlight), names(majority_cluster))]

# Confirm the names match
stopifnot(all(names(majority_cluster) == indNames(genlight)))

# Set new population based on majority cluster
pop(genlight) <- as.factor(majority_cluster)
final.dapc <- dapc(genlight, n.pca = 8, n.da = 3)

# Scatter plot
scatter(final.dapc, col = mycolors, cex = 2, legend = TRUE, clabel = T,
        posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft",
        cleg = 0.75, pch = 19, axesell = TRUE, label.inds = TRUE,
        xlab = "Discriminant Function 1", ylab = "Discriminant Function 2",
        main = "Final DAPC using Majority Cluster", sub = "Based on 100 replicates", cstar = 1)

# DAPC (final)
pnw.dapc <- dapc(genlight, n.pca = 8, n.da = 3)
scatter(pnw.dapc, col = mycolors, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, 
        posi.pca = "topleft", cleg = 0.75, pch = 19, axesell = TRUE, label.inds = TRUE,
        xlab = "Discriminant Function 1", ylab = "Discriminant Function 2", main = "DAPC Scatter Plot", 
        sub = "Colored by Group", cstar = 1, xax = 1, yax = 3)

# DAPC Posterior Probabilities
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(genlight)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

pdf("ag3.solanaceae.only.1079.all.host.population.structure.final.pdf", height = 10, width = 12)
ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = mycolors) +
  facet_grid(~Original_Pop, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))
dev.off()



