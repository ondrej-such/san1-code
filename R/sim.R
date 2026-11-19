# Simulation for the Section 3.2 replacement of satimage: N classes on the
# unit sphere S^2, von Mises-Fisher distributed around an optimal Tammes
# configuration of N points (R/tammes_configs.R), with a common concentration
# kappa. Minimal one-vs-one models are evaluated for spanning trees of K_N:
# all of them when there are few, otherwise a random sample, and then at most
# `max_per_dia` trees of each diameter per run (as in gam.R).

# Tree helpers shared with the benchmark code live in temperature_scaling.R,
# next to this file (R/ when the notebook runs from the repository root).
source(c("R/temperature_scaling.R", "temperature_scaling.R")[
  file.exists(c("R/temperature_scaling.R", "temperature_scaling.R"))][1])

# ---------------------------------------------------------------------------
# Configurations
# ---------------------------------------------------------------------------
tammes_means <- function(N) tammes_configs[[as.character(N)]]

# Angles (degrees) between all pairs of class means.
pair_angles <- function(means) {
  acos(pmax(pmin(means %*% t(means), 1), -1)) * 180 / pi
}

# Logical matrix: TRUE for the pairs of classes that are farthest apart.
most_distant_pairs <- function(means, tol = 0.01) {
  A <- pair_angles(means)
  far <- abs(A - max(A[upper.tri(A)])) < tol
  diag(far) <- FALSE
  far
}

# ---------------------------------------------------------------------------
# von Mises-Fisher sampling on S^2 (exact inversion for p = 3)
# ---------------------------------------------------------------------------
rvmf3 <- function(n, mu, kappa) {
  mu <- mu / sqrt(sum(mu^2))
  u  <- runif(n)
  # w = cosine of the angle to mu; numerically stable form of the inverse CDF
  w  <- 1 + log(u + (1 - u) * exp(-2 * kappa)) / kappa
  w  <- pmin(pmax(w, -1), 1)
  th <- runif(n, 0, 2 * pi)
  s  <- sqrt(1 - w^2)
  # orthonormal basis (a, b) of the plane orthogonal to mu
  helper <- if (abs(mu[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
  a <- helper - sum(helper * mu) * mu; a <- a / sqrt(sum(a^2))
  b <- c(mu[2] * a[3] - mu[3] * a[2], mu[3] * a[1] - mu[1] * a[3], mu[1] * a[2] - mu[2] * a[1])
  outer(s * cos(th), a) + outer(s * sin(th), b) + outer(w, mu)
}

simulate_classes <- function(n_per_class, kappa, means) {
  X <- do.call(rbind, lapply(seq_len(nrow(means)), function(k) rvmf3(n_per_class, means[k, ], kappa)))
  list(X = X, y = rep(seq_len(nrow(means)), each = n_per_class))
}

# ---------------------------------------------------------------------------
# Pairwise log-odds: L[n, i, j] = log P(i | x_n) / P(j | x_n), equal priors
# ---------------------------------------------------------------------------
logodds_bayes <- function(train, Xte, kappa, means) {
  K <- nrow(means); L <- array(0, c(nrow(Xte), K, K))
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    s <- kappa * as.vector(Xte %*% (means[i, ] - means[j, ]))
    L[, i, j] <- s; L[, j, i] <- -s
  }
  L
}

logodds_lda <- function(train, Xte, ...) {
  K <- max(train$y); L <- array(0, c(nrow(Xte), K, K))
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    Xi <- train$X[train$y == i, , drop = FALSE]; Xj <- train$X[train$y == j, , drop = FALSE]
    mi <- colMeans(Xi); mj <- colMeans(Xj)
    Sp <- ((nrow(Xi) - 1) * cov(Xi) + (nrow(Xj) - 1) * cov(Xj)) / (nrow(Xi) + nrow(Xj) - 2)
    w  <- solve(Sp, mi - mj)
    b  <- -sum(w * (mi + mj) / 2) + log(nrow(Xi) / nrow(Xj))
    s  <- as.vector(Xte %*% w + b)
    L[, i, j] <- s; L[, j, i] <- -s
  }
  L
}

# Firth (bias-reduced) logistic regression: maximises the log-likelihood
# penalised by the Jeffreys prior, so the estimates stay finite under linear
# separation (Firth 1993; Heinze and Schemper 2002). Newton iterations on the
# modified score X'(y - p + h (1/2 - p)), with step halving.
firth_fit <- function(x, y, maxit = 200, tol = 1e-9) {
  X <- cbind(1, x); b <- rep(0, ncol(X))
  for (it in seq_len(maxit)) {
    p <- plogis(drop(X %*% b)); W <- p * (1 - p)
    XtWXi <- solve(crossprod(X, X * W))
    h <- W * rowSums((X %*% XtWXi) * X)
    step <- drop(XtWXi %*% crossprod(X, y - p + h * (0.5 - p)))
    m <- max(abs(step)); if (m > 5) step <- step * 5 / m          # cap the step length
    b <- b + step
    if (max(abs(step)) < tol) break
  }
  b
}

logodds_firth <- function(train, Xte, ...) {
  K <- max(train$y); L <- array(0, c(nrow(Xte), K, K))
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    ix <- train$y %in% c(i, j)
    b  <- firth_fit(train$X[ix, , drop = FALSE], as.integer(train$y[ix] == i))
    s  <- drop(cbind(1, Xte) %*% b)
    L[, i, j] <- s; L[, j, i] <- -s
  }
  L
}

pairwise_logodds <- function(classifier, train, Xte, kappa, means) {
  switch(classifier,
         bayes    = logodds_bayes(train, Xte, kappa, means),
         lda      = logodds_lda(train, Xte),
         firth    = logodds_firth(train, Xte))
}

# ---------------------------------------------------------------------------
# Spanning trees via Pruefer codes
# ---------------------------------------------------------------------------

# Candidate spanning trees of K_K: all K^(K-2) of them when there are at most
# `max_enumerate`. Otherwise `n_sample` random Pruefer codes, plus all K stars
# (constant codes) and up to `n_dia3` random trees of diameter 3 (codes with
# exactly two distinct labels), which uniform sampling would almost never hit.
# `far` marks the most distant class pairs (counted per tree).
candidate_trees <- function(K, far, max_enumerate = 20000, n_sample = 20000,
                            n_dia3 = 2000, seed = 1) {
  if (K == 2) {
    codes <- matrix(integer(0), 1, 0)
  } else if (K^(K - 2) <= max_enumerate) {
    codes <- as.matrix(expand.grid(rep(list(1:K), K - 2)))
  } else {
    set.seed(seed)
    random <- matrix(sample.int(K, n_sample * (K - 2), replace = TRUE), ncol = K - 2)
    stars  <- matrix(rep(1:K, each = K - 2), ncol = K - 2, byrow = TRUE)
    dia3   <- t(replicate(n_dia3, {
      centres <- sample.int(K, 2)
      repeat { pick <- sample(centres, K - 2, replace = TRUE)
               if (length(unique(pick)) == 2) break }
      pick
    }))
    codes <- unique(rbind(stars, dia3, random))
  }
  trees <- lapply(seq_len(nrow(codes)), function(r) prufer_edges(as.integer(codes[r, ]), K))
  structure(
    data.frame(tree = if (K == 2) "1-2" else apply(codes, 1, paste, collapse = ":"),
               dia  = vapply(trees, tree_diameter, numeric(1), K = K),
               far_edges = vapply(trees, function(e) sum(far[e]), numeric(1)),
               id = seq_along(trees)),
    edges = trees, enumerated = K == 2 || K^(K - 2) <= max_enumerate)
}

# Predicted class of the minimal model of one tree: log p_k relative to class 1
# accumulated along tree paths, then argmax (equivalent to log_coupling()).
tree_predict <- function(edges, L) {
  K <- dim(L)[2]; n <- dim(L)[1]
  S <- matrix(NA_real_, n, K); S[, 1] <- 0
  adj <- lapply(1:K, function(v) c(edges[edges[, 1] == v, 2], edges[edges[, 2] == v, 1]))
  queue <- 1; seen <- 1
  while (length(queue)) {
    p <- queue[1]; queue <- queue[-1]
    for (c in setdiff(adj[[p]], seen)) {
      S[, c] <- S[, p] + L[, c, p]
      seen <- c(seen, c); queue <- c(queue, c)
    }
  }
  max.col(S, ties.method = "first")
}

# ---------------------------------------------------------------------------
# One experiment: `reps` simulated datasets; in each, up to `max_per_dia`
# random candidate trees of each diameter (as in gam.R), or all if NULL. With
# calibrate = TRUE the accuracy of the sequentially temperature-scaled minimal
# model is returned in `acc_cal` (NA otherwise).
# ---------------------------------------------------------------------------
run_experiment <- function(kappa, means, trees, classifier = "lda", n_train = 50,
                           n_test = 200, reps = 20, max_per_dia = 50, seed = 123,
                           calibrate = FALSE, progress = NULL) {
  set.seed(seed)
  edges <- attr(trees, "edges")
  out <- vector("list", reps)
  for (r in seq_len(reps)) {
    train <- simulate_classes(n_train, kappa, means)
    test  <- simulate_classes(n_test, kappa, means)
    L     <- pairwise_logodds(classifier, train, test$X, kappa, means)
    bayes <- max.col(test$X %*% t(means), ties.method = "first")   # Bayes rule
    Loof  <- if (calibrate) pairwise_logodds_oof(classifier, train, kappa, means)
    ids <- if (is.null(max_per_dia)) trees$id else
      unlist(lapply(split(trees$id, trees$dia), function(v)
        if (length(v) > max_per_dia) sample(v, max_per_dia) else v))
    acc <- vapply(ids, function(t) mean(tree_predict(edges[[t]], L) == test$y), numeric(1))
    acc_cal <- if (calibrate) vapply(ids, function(t) {
      a <- fit_tree_scales_seq(edges[[t]], Loof, train$y)
      mean(tree_predict_scaled(edges[[t]], L, a) == test$y)
    }, numeric(1)) else rep(NA_real_, length(ids))
    out[[r]] <- data.frame(trees[ids, c("tree", "dia", "far_edges")], acc = acc, acc_cal = acc_cal,
                           run = r, bayes_acc = mean(bayes == test$y))
    if (!is.null(progress)) progress(r / reps)
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# Per-edge temperature scaling (as in the preprint on multiplicative
# miscalibration). The log-odds of the classifier on edge e are multiplied by
# a_e = 1 / T_e. The scales are chosen, separately for each tree, to maximise
# the likelihood of the minimal model on OUT-OF-FOLD log-odds of the training
# data. Calibrating each binary classifier on its own cannot see
# multiplicative miscalibration, so the multiclass likelihood is used.
# ---------------------------------------------------------------------------

# Out-of-fold pairwise log-odds of the training points: Loof[n, i, j].
pairwise_logodds_oof <- function(classifier, train, kappa, means, folds = 5) {
  n <- length(train$y); K <- nrow(means)
  fold <- integer(n)
  for (k in seq_len(K)) { ix <- which(train$y == k); fold[ix] <- sample(rep_len(seq_len(folds), length(ix))) }
  Loof <- array(0, c(n, K, K))
  for (f in seq_len(folds)) {
    tr <- list(X = train$X[fold != f, , drop = FALSE], y = train$y[fold != f])
    Loof[fold == f, , ] <- pairwise_logodds(classifier, tr, train$X[fold == f, , drop = FALSE], kappa, means)
  }
  Loof
}

