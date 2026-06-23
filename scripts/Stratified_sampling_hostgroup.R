library(dplyr)
library(Biostrings)
library(adegenet)

setwd("F:/UNL/Research/outgroup/athelia/4121")

metadata <- read.csv(file = "F:/UNL/Research/R.solani.published.ITS.host.country.year.4261.renamed.filtered.seqs.complete.metadata.csv")

genlight <- fasta2genlight(file = "R.solani.ITS.4121a.seqs.mafft.fftns2.fasta")

fastafile<- readDNAStringSet(file = "R.solani.ITS.4121a.seqs.mafft.fftns2.fasta")

# F metadata
metadata.filtered <- metadata[metadata$accessions %in% genlight$ind.names,] #Extract 
length(metadata.filtered$accessions)

metadata.filtered.sorted <- metadata.filtered[match(genlight$ind.names, metadata.filtered$accessions),]

# Create 10 stratified subsamples with complete AG data
multiple_stratified_subsamples <- function(data, host_col, group_col, n_samples, major_host_groups, n_subsamples) {
  all_subsamples <- vector("list", n_subsamples)
  
  for (i in 1:n_subsamples) {
    subsamples <- lapply(major_host_groups, function(host_group) {
      
      # Filter data for the current host group and ensure strain.beast is not NA
      host_group_data <- data %>%
        filter(!!sym(host_col) == host_group, !is.na(!!sym(group_col)))
      
      # Sample 40 rows from the host group (or fewer if less available)
      if (nrow(host_group_data) >= n_samples) {
        sample_n(host_group_data, n_samples)
      } else {
        host_group_data  # Include all available samples if fewer than 40
      }
    })
    
    # Combine all subsamples into one data frame and store in the list
    all_subsamples[[i]] <- bind_rows(subsamples)
  }
  
  return(all_subsamples)
}


set.seed(123)  #For reproducibility

# Define parameters
host_col <- "host_groups"
group_col <- "strain.beast"  # Column to ensure complete data
n_samples <- 40  # Number of samples per group
n_subsamples <- 10  # Number of subsamples to generate
major_host_groups <- c("Solanaceae", "Amaranthaceae", "Fabaceae", "Brassicaceae", "Poaceae")

# Run the function to generate 10 subsamples
subsamples_list <- multiple_stratified_subsamples(metadata.filtered, host_col, group_col, n_samples, major_host_groups, n_subsamples)

# Check the result
print(subsamples_list[[4]])

# Check the AG composition within each subsample
subsample_composition <- bind_rows(
  lapply(seq_along(subsamples_list), function(i) {
    subsamples_list[[i]] %>%
      count(host_groups, strain.beast, name = "n") %>%
      group_by(host_groups) %>%
      mutate(
        total_host = sum(n),
        prop = n / total_host,
        subsample = i
      ) %>%
      ungroup()
  })
)

subsample_composition %>%
  arrange(subsample, host_groups, desc(n))

dominant_AG_summary <- subsample_composition %>%
  group_by(subsample, host_groups) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(subsample, host_groups, dominant_AG = strain.beast, n, total_host, prop)

dominant_AG_summary

# Write each subsample to a CSV file
write_subsamples_to_csv <- function(subsamples_list, prefix = "subsample.AG.host") {
  for (i in seq_along(subsamples_list)) {
    filename <- paste0(prefix, "_", i, ".csv")
    write.csv(subsamples_list[[i]], filename, row.names = FALSE)
  }
}

write_subsamples_to_csv(subsamples_list)

check <- read.csv("subsample.AG.host_2.csv")
