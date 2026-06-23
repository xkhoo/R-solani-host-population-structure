# install.packages(c("ape", "pegas", "dplyr", "stringr", "purrr", "readr"))

library(ape)
library(pegas)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

setwd("F:/UNL/Research/")

metadata <- read.csv(
  "Table S1_Rsolani_4119_seqs_metadata.csv",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE
)

read_alignment_with_metadata <- function(fasta_file, metadata, strip_after_underscore = FALSE) {
  dna <- read.dna(fasta_file, format = "fasta")
  seq_names <- rownames(dna)
  
  if (strip_after_underscore) {
    seq_names_clean <- sub("_.*$", "", seq_names)
  } else {
    seq_names_clean <- seq_names
  }
  
  rownames(dna) <- seq_names_clean
  
  meta <- metadata %>%
    filter(accessions %in% seq_names_clean) %>%
    distinct(accessions, .keep_all = TRUE)
  
  meta <- meta[match(rownames(dna), meta$accessions), ]
  
  if (!all(meta$accessions == rownames(dna))) {
    stop("Metadata and FASTA sequence names are not matching correctly.")
  }
  
  return(list(dna = dna, meta = meta))
}

# Subset by host group
subset_dna_by_group <- function(dna, meta, group_col, groups_keep) {
  keep <- meta[[group_col]] %in% groups_keep & !is.na(meta[[group_col]])
  dna_sub <- dna[keep, ]
  meta_sub <- meta[keep, ]
  
  meta_sub[[group_col]] <- factor(meta_sub[[group_col]], levels = groups_keep)
  
  return(list(dna = dna_sub, meta = meta_sub))
}

# Pairwise distance-based Fst (Hudson-like pairwise Fst using ITS genetic distance)
mean_within_distance <- function(D, idx) {
  if (length(idx) < 2) return(NA_real_)
  vals <- D[idx, idx]
  vals <- vals[upper.tri(vals)]
  mean(vals, na.rm = TRUE)
}

mean_between_distance <- function(D, idx1, idx2) {
  mean(D[idx1, idx2], na.rm = TRUE)
}

pairwise_distance_fst <- function(dna, meta, group_col, ag_name, nperm = 999, seed = 123) {
  set.seed(seed)
  
  D <- as.matrix(dist.dna(
    dna,
    model = "raw",
    pairwise.deletion = TRUE
  ))
  
  groups <- as.character(meta[[group_col]])
  group_levels <- levels(factor(groups))
  
  pair_list <- combn(group_levels, 2, simplify = FALSE)
  
  results <- map_dfr(pair_list, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    idx1 <- which(groups == g1)
    idx2 <- which(groups == g2)
    
    dxy <- mean_between_distance(D, idx1, idx2)
    dxx <- mean_within_distance(D, idx1)
    dyy <- mean_within_distance(D, idx2)
    
    fst_obs <- (dxy - 0.5 * (dxx + dyy)) / dxy
    
    # permutation test by shuffling group labels within this pair
    idx_pair <- c(idx1, idx2)
    groups_pair <- groups[idx_pair]
    D_pair <- D[idx_pair, idx_pair]
    
    perm_fst <- replicate(nperm, {
      perm_groups <- sample(groups_pair)
      
      p1 <- which(perm_groups == g1)
      p2 <- which(perm_groups == g2)
      
      dxy_p <- mean_between_distance(D_pair, p1, p2)
      dxx_p <- mean_within_distance(D_pair, p1)
      dyy_p <- mean_within_distance(D_pair, p2)
      
      (dxy_p - 0.5 * (dxx_p + dyy_p)) / dxy_p
    })
    
    p_value <- (sum(perm_fst >= fst_obs, na.rm = TRUE) + 1) /
      (sum(!is.na(perm_fst)) + 1)
    
    tibble(
      AG = ag_name,
      population_1 = g1,
      population_2 = g2,
      n_1 = length(idx1),
      n_2 = length(idx2),
      mean_within_1 = dxx,
      mean_within_2 = dyy,
      mean_between = dxy,
      pairwise_Fst = fst_obs,
      p_value = p_value
    )
  })
  
  return(results)
}

# AMOVA
run_amova <- function(dna, meta, group_col, ag_name, nperm = 999) {
  D <- dist.dna(
    dna,
    model = "raw",
    pairwise.deletion = TRUE
  )
  
  dat <- data.frame(
    group = factor(meta[[group_col]])
  )
  
  amova_result <- amova(D ~ group, data = dat, nperm = nperm)
  
  sink(paste0(ag_name, "_AMOVA_output.txt"))
  cat("AMOVA for", ag_name, "\n\n")
  print(amova_result)
  sink()
  
  return(amova_result)
}


# AG-3 population: Potato vs Tobacco vs Non-Solanaceae
ag3_pop <- read_alignment_with_metadata(
  fasta_file = "haplotype/AG3.914.seqs.mafft.fftns2.replaced.renamed.ori.fasta",
  metadata = metadata,
  strip_after_underscore = TRUE
)

ag3_groups_keep <- c("Potato", "Tobacco", "Others", "Non_solanaceae")

ag3 <- subset_dna_by_group(
  dna = ag3_pop$dna,
  meta = ag3_pop$meta,
  group_col = "host_solanaceae",
  groups_keep = ag3_groups_keep
)

# Check sample sizes
table(ag3$meta$host_solanaceae)

ag3_pairwise_fst <- pairwise_distance_fst(
  dna = ag3$dna,
  meta = ag3$meta,
  group_col = "host_solanaceae",
  ag_name = "AG-3",
  nperm = 999
)

ag3_amova <- run_amova(
  dna = ag3$dna,
  meta = ag3$meta,
  group_col = "host_solanaceae",
  ag_name = "AG3",
  nperm = 999
)


# AG-4 population: major host families + Other
ag4_pop <- read_alignment_with_metadata(
  fasta_file = "haplotype/AG4/R.solani.AG4.1175.seqs.mafft.fftns2.gblocks.replaced.fasta",
  metadata = metadata,
  strip_after_underscore = FALSE
)

ag4_groups_keep <- c(
  "Amaranthaceae",
  "Brassicaceae",
  "Fabaceae",
  "Poaceae",
  "Solanaceae",
  "Other"
)

ag4 <- subset_dna_by_group(
  dna = ag4_pop$dna,
  meta = ag4_pop$meta,
  group_col = "host_group_summarized",
  groups_keep = ag4_groups_keep
)

# Check sample sizes
table(ag4$meta$host_group_summarized)

ag4_pairwise_fst <- pairwise_distance_fst(
  dna = ag4$dna,
  meta = ag4$meta,
  group_col = "host_group_summarized",
  ag_name = "AG-4",
  nperm = 999
)

ag4_amova <- run_amova(
  dna = ag4$dna,
  meta = ag4$meta,
  group_col = "host_group_summarized",
  ag_name = "AG4",
  nperm = 999
)


all_pairwise_fst <- bind_rows(ag3_pairwise_fst, ag4_pairwise_fst)

write.csv(
  all_pairwise_fst,
  "AG3_AG4_pairwise_distance_Fst_by_host_group.csv",
  row.names = FALSE
)

# View
all_pairwise_fst
