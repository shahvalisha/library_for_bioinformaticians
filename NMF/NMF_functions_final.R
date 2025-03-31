############# Functions to run and validate NMF clusters 

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






