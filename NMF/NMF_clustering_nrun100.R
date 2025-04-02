############# Functions to run and validate NMF clusters 

# Perform NMF and calculate reconstruction error
run_nmf <- function(data, rank, method = "brunet", nrun = 10) {
  nmf_result <- nmf(data, rank, method = method, nrun = nrun)
  basis_matrix <- basis(nmf_result)  # H matrix
  coefficient_matrix <- coef(nmf_result)  # W matrix
  
  return(list(nmf_result = nmf_result, basis_matrix = basis_matrix, coefficient_matrix = coefficient_matrix))
}

# Include cophenetic score calculation
run_consensus_clustering <- function(nmf_result, max_k = 6, reps = 100) {
  basis_matrix <- basis(nmf_result)  # H matrix from NMF
  
  # Consensus clustering with ConsensusClusterPlus
  consensus_result <- ConsensusClusterPlus(
    t(basis_matrix),  # Transpose because samples are in rows
    maxK = max_k,
    reps = reps,
    pItem = 0.8, 
    pFeature = 1,
    clusterAlg = "hc",  # Use hierarchical clustering on NMF basis matrix
    distance = "euclidean"
  )
  
  # Calculate cophenetic correlation score
  # Extract consensus matrix for the chosen K
  consensus_matrix <- consensus_result[[max_k]]$consensusMatrix
  
  # Perform hierarchical clustering on consensus matrix
  hc <- hclust(as.dist(1 - consensus_matrix), method = "average")  # 1 - consensus matrix to create a distance matrix
  
  # Compute cophenetic correlation score
  cophenetic_corr <- cor(cophenetic(hc), as.dist(1 - consensus_matrix))
  
  return(list(consensus_result = consensus_result, cophenetic_corr = cophenetic_corr))
}



library(tidyverse)
library(ConsensusClusterPlus)
library(cluster)
library(clusterProfiler)
library(NMF)

corrected_tpm <- read.delim("/path/to/corrected_log_denovobatch_tpm_withneg.tsv")

rownames(corrected_tpm) <- corrected_tpm$gene_id
numerical_tpm <- corrected_tpm %>%
  dplyr::select(-gene_id)
transcript_ids <- rownames(corrected_tpm)
unique_transcript_ids <- make.unique(transcript_ids)
rownames(corrected_tpm) <- unique_transcript_ids
allzero <- apply(numerical_tpm, 1, sum)
tpm <- numerical_tpm[allzero != 0, ]

## 2. Top 10% highest variant genes – 3075
variance <- data.frame(gene = rownames(tpm), mad=apply(tpm, 1, mad))
variance <- variance[order(-variance$mad), ]
# keep top 10%
keep_var <- variance[1:5500, ]
# keep top 10% high variable genes
tpm <- tpm[which(row.names(tpm) %in% keep_var$gene), ]

# Method 1: Set negative values to zero
tpm <- as.matrix(tpm)

# Method 2: Add a constant to make all values non-negative
# Choose a small constant based on your data
tpm <- tpm + abs(min(tpm))

data.frame(tpm)

# Initialize lists to store the times taken for each loop
nmf_time_taken <- list()
consensus_time_taken <- list()
frobenius_time_taken <- list()
nmf_results <- list()
consensus_results <- list()
frobenius_error <- list()

# NMF and Consensus Clustering Loop
for (k in 4:9) {
  print(k)
  # Record the start time
  nmf_start_time <- system.time({
    nmf_results[[k]] <- run_nmf(t(tpm), rank = k)
  })
  nmf_time_taken[[k]] <- nmf_start_time
  
  # Perform reconstruction and consensus clustering
  reconstruction <- nmf_results[[k]][[2]] %*% nmf_results[[k]][[3]]
  
  consensus_start_time <- system.time({
    consensus_results[[k]] <- run_consensus_clustering(nmf_results[[k]][[1]])
  })
  consensus_time_taken[[k]] <- consensus_start_time
  
  # Print the cophenetic correlation
  print(consensus_results[[k]]$cophenetic_corr)
}

# Frobenius Error Calculation Loop
for (k in 4:9) {
  print(k)
  frobenius_start_time <- system.time({
    reconstruction <- nmf_results[[k]][[2]] %*% nmf_results[[k]][[3]]
    # Calculate the Frobenius norm
    frobenius_error[[k]] <- norm(as.matrix(t(tpm)) - as.matrix(reconstruction), type = "F")
  })
  frobenius_time_taken[[k]] <- frobenius_start_time
}

# Print the recorded times
print(nmf_time_taken)
print(consensus_time_taken)
print(frobenius_time_taken)

save.image("NMF.RData")







