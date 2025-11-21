# Step 1: Define function to compute Fisher's test for a single node
compute_fisher_for_node <- function(node, go_to_orthogroups, node_to_orthogroups, node_to_go, all_orthogroups) {
  node_go <- unique(node_to_go[[node]])
  node_orthogroups <- unique(node_to_orthogroups[[node]])
  results <- list()
  
  count <- 0  # Progress tracking for the node
  total_go_terms <- length(go_to_orthogroups)
  
  for (go_term in node_go) {
    go_orthogroups <- unique(go_to_orthogroups[[go_term]])
    # Create contingency table
    a <- length(intersect(node_orthogroups, go_orthogroups))  # At node & has GO term
    b <- length(setdiff(node_orthogroups, go_orthogroups))    # At node & no GO term
    c <- length(setdiff(go_orthogroups, node_orthogroups))    # Not at node & has GO term
    d <- length(setdiff(all_orthogroups, union(node_orthogroups, go_orthogroups)))  # Not at node & no GO term
    if (a >= 0) {
      # Fisher's test
      fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      
      # Store results
      results[[length(results) + 1]] <- list(
        Node = node,
        GO_Term = go_term,
        At_Node_Has_GO = a,
        At_Node_No_GO = b,
        Not_At_Node_Has_GO = c,
        Not_At_Node_No_GO = d,
        Odds_Ratio = fisher_result$estimate,
        P_Value = fisher_result$p.value
      )
      
      # Update progress
      count <- count + 1
      if (count %% 100 == 0 || count == total_go_terms) {
        cat(sprintf("Node %s: Processed %d/%d GO terms\n", node, count, total_go_terms))
      }
    }
  }
  
  return(rbindlist(results))
}

run_fast_fishers <- function(infile = NULL, outfile = NULL) {
  try(if(is.null(infile)) stop("infile needs to be a file name including a path"))
  require(tidyverse)
  require(data.table)
  require(parallel)
  
  data <- fread(infile)
  cat("Input file loaded...\n")
  
  # Step 2: Create mappings of nodes to orthogroups and GO terms to orthogroups
  node_to_orthogroups <- data %>% distinct(HOG, node) %>% with(., split(HOG, node))
  go_to_orthogroups <- data %>% distinct(HOG, value) %>% with(., split(HOG, value))
  node_to_go <- data %>% distinct(value, node) %>% with(., split(value, node))
  cat("Mappings created...\n")
  
  # Get all unique orthogroups
  all_orthogroups <- unique(data$HOG)
  
  # Step 3: Parallel computation for all nodes
  all_nodes <- names(node_to_orthogroups)
  num_cores <- 8  # Use one less core than available
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("node_to_orthogroups", "node_to_go", "go_to_orthogroups", "all_orthogroups", "compute_fisher_for_node"))
  clusterEvalQ(cl, library(data.table))
  
  # Apply the function to each node in parallel
  all_results <- parLapply(cl, all_nodes, function(node) {
    compute_fisher_for_node(node = node, node_to_orthogroups = node_to_orthogroups, go_to_orthogroups = go_to_orthogroups, node_to_go = node_to_go, all_orthogroups = all_orthogroups)
  })
  stopCluster(cl)
  
  # Combine all results
  final_results <- rbindlist(all_results)
  
  # Step 4: Save results to file
  if(is.null(outfile)) {
    frwite(final_results, "fishers_test_out.csv.gz", compress = "auto")
  } else {
    fwrite(final_results, outfile, compress = "auto")
  }
}


