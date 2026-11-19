
library(igraph)
library(furrr)
library(dplyr)
library(purrr) 
library(stats)

cohens_d <- function(Scores_i, Scores_j) {
  mu_i <- mean(Scores_i)
  mu_j <- mean(Scores_j)
  
  var_i <- var(Scores_i)
  var_j <- var(Scores_j)
  n_i <- length(Scores_i)
  n_j <- length(Scores_j)
  
  numerator_pooled_sd <- (n_i - 1) * var_i + (n_j - 1) * var_j
  denominator_pooled_sd <- n_i + n_j - 2
  
  s_pooled <- sqrt(numerator_pooled_sd / denominator_pooled_sd)
  
  cohens_d <- abs(mu_i - mu_j) / s_pooled
  
  return(cohens_d)
} 

compute_pairwise_d <- function(log_exp1_result) {
  K <- log_exp1_result$K
  pred <- log_exp1_result$pred
  GT <- log_exp1_result$ground_truth
  Cd_list <- list()
  
  pair_count <- 0

  M <- matrix(nrow = K * (K-1) / 2, ncol = 3)
  row <- 1
  
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      pair_count <- pair_count + 1
      pair_indices <- which(GT %in% c(i, j))
      scores_all <- pred[i, j, pair_indices]
      true_class_labels <- GT[pair_indices]
      Scores_i <- scores_all[true_class_labels == i]
      Scores_j <- scores_all[true_class_labels == j]
   
      stopifnot (length(Scores_i) > 1 && length(Scores_j) > 1)
      Cd_value <- cohens_d(Scores_i, Scores_j)
      M[row,] = c(i,j, Cd_value)
      row <- row + 1
    }
  }
  stopifnot(row - 1 == nrow(M))
  
  # Cd_vector <- unlist(Cd_list)
  # return(Cd_vector)
}



Cd_for_all<- function(dataset_name = "dna", K, max_tree = 15000, max_dia = 50) {
  plan("multisession", workers = 4)
  all_results <- future_map(0:19, function(run) {
    set.seed(123)
    
    e1 <- log_exp1(800, dataset = dataset_name, run)
    n <- e1$K
    sequences <- expand.grid(rep(list(1:n), n - 2))
    dp_vector <- compute_pairwise_d(e1) 
    Cd_df <- as.data.frame(as.list(dp_vector))
    
    run_results <- map(1:nrow(sequences), function(t) {
      ptree <- as.integer(sequences[t, ])
      g <- make_from_prufer(ptree) |> igraph::graph_from_edgelist(directed = FALSE)
      g_rad <- radius(g)
      g_dia <- diameter(g)
      g_apl <- mean_distance(g)  # Use igraph >= 2.0.0
      
      stree <- paste(sprintf("%d", ptree), collapse = ":")
      probs3 <- log_coupling(ptree, e1$pred)
      hard_decision3 <- apply(probs3, 1, which.max)
      acc <- mean(hard_decision3 == e1$ground_truth)
      
      result_row <- data.frame(
        tree = stree,
        rad = g_rad,
        dia = g_dia,
        apl = g_apl,
        acc = acc,
        run = run,
        dataset = dataset_name
      )
      
      cbind(result_row, Cd_df)
    }) |> list_rbind()
    run_results
  }, .options = furrr_options(seed = TRUE)) |> list_rbind()
  
  best_rows <- all_results %>%
    group_by(run) %>%
    slice_max(acc, n = 1, with_ties = FALSE) %>%
    ungroup()
  write.csv(best_rows, file = sprintf("data/best_dp_%s.csv", dataset_name), row.names = FALSE)
}