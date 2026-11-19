# Example 1 of the paper (Figure 3, Figure 5): three bivariate Gaussian classes
# with means mu1 = (0,0), mu2 = (0,3), mu3 = (4,0) and the common covariance
# Sigma_r = [[1, r], [r, 1]]. Copied from main.R of the paper repository.

mu1 <- c(0, 0)
mu2 <- c(0, 3)
mu3 <- c(4, 0)
mus <- rbind(mu1, mu2, mu3)
Sigma <- function(r) matrix(c(1, r, r, 1), ncol = 2)

# A dataset of N points per class, split into training and test parts
# (caret::createDataPartition, stratified by class).
dfs_for_r <- function(r, N, p_train = 0.8) {
  dataset <- purrr::map(1:3, function(i) {
    samples <- MASS::mvrnorm(n = N, mu = mus[i, ], Sigma = Sigma(r))
    df <- data.frame(x = samples[, 1], y = samples[, 2])
    df$class_id <- i
    df
  }) |> purrr::list_rbind()
  dataset$class_id <- as.factor(dataset$class_id)
  train_index <- caret::createDataPartition(dataset$class_id, p = p_train, list = FALSE)
  list(train = dataset[train_index, ], test = dataset[-train_index, ],
       r = r, mus = mus, sigma = Sigma(r))
}

# Population Fisher ratio (squared Mahalanobis distance of the class means) of
# each pair of classes.
compute_population_fisher <- function(r) {
  Sigma_inv <- solve(Sigma(r))
  pairs <- list(`12` = list(mu1, mu2), `13` = list(mu1, mu3), `23` = list(mu2, mu3))
  fisher_vals <- sapply(pairs, function(m) {
    d <- m[[1]] - m[[2]]
    as.numeric(t(d) %*% Sigma_inv %*% d)
  })
  names(fisher_vals) <- c("C1 and C2", "C1 and C3", "C2 and C3")
  fisher_vals
}
