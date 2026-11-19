# Trees of diameter 3 for the third heuristic (Section 3.3).
#
# Every spanning tree of K_N of diameter 3 is a double star S_{p,q}: two
# adjacent centres c1, c2 with p and q leaves, p, q >= 1, p + q = N - 2.
# Its paths of length 3 join a leaf of c1 to a leaf of c2, so Omega_3 = p * q,
# and the average path length depends only on the split {p, q}:
#   apl = (2 * C - (N - 1) + p * q) / C,  where C = choose(N, 2).
# Uses the simulation code of the Section 3.2 app (vMF classes, classifiers,
# tree_predict) from ../sec32-app/R.

# All unordered splits {p, q} with p <= q, with their path statistics and the
# number of spanning trees of K_N that have this split.
diam3_splits <- function(N) {
  stopifnot(N >= 4)
  C <- choose(N, 2)
  p <- seq_len(floor((N - 2) / 2)); q <- N - 2 - p
  data.frame(p = p, q = q, omega3 = p * q,
             apl = (2 * C - (N - 1) + p * q) / C,
             trees = C * ifelse(p == q, choose(N - 2, p), 2 * choose(N - 2, p)))
}

# Edges of a random double star with split {p, q}: uniform over all such trees
# (random ordered pair of centres, random p of the others attached to c1).
random_double_star <- function(N, p) {
  centres <- sample.int(N, 2)
  others  <- setdiff(seq_len(N), centres)
  to_c1   <- if (length(others) == p) others else sample(others, p)
  to_c2   <- setdiff(others, to_c1)
  rbind(c(centres[1], centres[2]),
        cbind(centres[1], to_c1),
        cbind(centres[2], to_c2))
}

# One experiment: `reps` simulated datasets; in each, `per_split` random
# double stars for every split {p, q}.
run_diam3_experiment <- function(kappa, means, classifier = "lda", n_train = 100,
                                 n_test = 200, reps = 20, per_split = 50,
                                 seed = 123, calibrate = FALSE, progress = NULL) {
  set.seed(seed)
  N <- nrow(means)
  splits <- diam3_splits(N)
  far <- most_distant_pairs(means)
  out <- vector("list", reps)
  for (r in seq_len(reps)) {
    train <- simulate_classes(n_train, kappa, means)
    test  <- simulate_classes(n_test, kappa, means)
    L     <- pairwise_logodds(classifier, train, test$X, kappa, means)
    bayes <- max.col(test$X %*% t(means), ties.method = "first")
    # out-of-fold log-odds for temperature scaling (see sim.R); acc_cal is NA otherwise
    Loof  <- if (calibrate) pairwise_logodds_oof(classifier, train, kappa, means)
    out[[r]] <- do.call(rbind, lapply(seq_len(nrow(splits)), function(s) {
      p <- splits$p[s]
      do.call(rbind, lapply(seq_len(per_split), function(t) {
        e <- random_double_star(N, p)
        data.frame(p = p, q = splits$q[s], apl = splits$apl[s],
                   far_edges = sum(far[e]),
                   acc = mean(tree_predict(e, L) == test$y),
                   acc_cal = if (calibrate)
                     mean(tree_predict_scaled(e, L, fit_tree_scales_seq(e, Loof, train$y)) == test$y)
                   else NA_real_)
      }))
    }))
    out[[r]]$run <- r
    out[[r]]$bayes_acc <- mean(bayes == test$y)
    if (!is.null(progress)) progress(r / reps)
  }
  do.call(rbind, out)
}
