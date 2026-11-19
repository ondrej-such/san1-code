# Minimal one-vs-one models on three bivariate Gaussian classes with LDA
# pairwise classifiers (Examples 1 and 2 of paper.qmd).
#
# Classes: means mu1 = (0,0), mu2 = (0,3), mu3 = (4,0), common covariance
# Sigma_r = [[1, r], [r, 1]].  For K = 3 every spanning tree is a star; the
# minimal model Mc is the star with centre c (Prufer code c, as in log_coupling):
#   M1 = edges 1-2, 1-3;  M2 = 1-2, 2-3;  M3 = 1-3, 2-3.

suppressPackageStartupMessages({
  library(MASS); library(dplyr); library(future); library(furrr)
})

lda_mus   <- rbind(c(0, 0), c(0, 3), c(4, 0))
lda_Sigma <- function(r) matrix(c(1, r, r, 1), 2)
lda_pairs <- list(c(1, 2), c(1, 3), c(2, 3))

# Pairwise LDA log-odds L[i, j, ] = log P(i | x) / P(j | x) on the test points,
# using the pooled covariance of the pair.  Computed in closed form so that
# predictions far from the boundary do not saturate to 0 or 1.
lda_pairwise_logodds <- function(train, Xte) {
  L <- array(0, c(3, 3, nrow(Xte)))
  for (p in lda_pairs) {
    i <- p[1]; j <- p[2]
    Xi <- train[[i]]; Xj <- train[[j]]
    mi <- colMeans(Xi); mj <- colMeans(Xj)
    Sp <- ((nrow(Xi) - 1) * cov(Xi) + (nrow(Xj) - 1) * cov(Xj)) /
          (nrow(Xi) + nrow(Xj) - 2)
    w  <- solve(Sp, mi - mj)
    b  <- -sum(w * (mi + mj) / 2) + log(nrow(Xi) / nrow(Xj))
    s  <- as.vector(Xte %*% w + b)
    L[i, j, ] <- s; L[j, i, ] <- -s
  }
  L
}

# Accuracy of M1, M2, M3 on one simulated dataset with n_train training and
# n_test test samples per class.
lda_one_rep <- function(r, n_train, n_test) {
  S <- lda_Sigma(r)
  train <- lapply(1:3, function(k) MASS::mvrnorm(n_train, lda_mus[k, ], S))
  test  <- lapply(1:3, function(k) MASS::mvrnorm(n_test,  lda_mus[k, ], S))
  Xte <- do.call(rbind, test)
  y   <- rep(1:3, each = n_test)
  L   <- lda_pairwise_logodds(train, Xte)
  # star with centre c: log p_k - log p_c = L[k, c, ]; predict the argmax
  sapply(1:3, function(c) {
    scores <- sapply(1:3, function(k) if (k == c) rep(0, nrow(Xte)) else L[k, c, ])
    mean(max.col(scores, ties.method = "first") == y)
  })
}

# Accuracies of M1-M3 over `reps` datasets at a single r (long format).
lda_example1 <- function(r = 0, reps = 1000, n_train = 50, n_test = 150, seed = 123) {
  set.seed(seed)
  a <- replicate(reps, lda_one_rep(r, n_train, n_test))      # 3 x reps
  data.frame(rep = rep(1:reps, each = 3), model = paste0("M", 1:3),
             accuracy = as.vector(a))
}

# Mean accuracy (and its standard error) of M1-M3 for each r, in parallel over r.
lda_example2 <- function(rs = seq(-0.9, 0.9, length.out = 50), reps = 1000,
                         n_train = 40, n_test = 10, seed = 123,
                         workers = max(1, parallel::detectCores() - 1)) {
  plan(multisession, workers = workers)
  on.exit(plan(sequential))
  future_map_dfr(rs, function(r) {
    a <- replicate(reps, lda_one_rep(r, n_train, n_test))
    data.frame(r = r, model = paste0("M", 1:3),
               accuracy = rowMeans(a), se = apply(a, 1, sd) / sqrt(reps))
  }, .options = furrr_options(seed = seed))
}

# Population Cohen's d (Mahalanobis distance between class means) of each pair.
lda_population_d <- function(rs) {
  bind_rows(lapply(rs, function(r) {
    Si <- solve(lda_Sigma(r))
    bind_rows(lapply(lda_pairs, function(p) {
      v <- lda_mus[p[1], ] - lda_mus[p[2], ]
      data.frame(r = r, pair = paste0(p[1], "-", p[2]),
                 d = sqrt(drop(t(v) %*% Si %*% v)))
    }))
  }))
}
