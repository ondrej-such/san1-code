# Sequential temperature scaling of the pairwise classifiers of a minimal
# one-vs-one model (preprint on multiplicative miscalibration, Sections 5.5-5.6).
# Shared by the benchmark experiments (exp4.R) and the Shiny apps (sec32-app,
# sec33-app). All trees are given as edge matrices (K - 1 rows, 2 columns) and
# log-odds as arrays L[n, i, j] = log P(i | x_n) / P(j | x_n).
#
# The log-odds of the classifier on edge e are multiplied by a_e = 1 / T_e. The
# scales are chosen for each tree on OUT-OF-FOLD log-odds of the training data,
# from the multiclass likelihood, because calibrating each binary classifier on
# its own cannot see multiplicative miscalibration.

# Pruefer code -> edge matrix (K - 1 rows).
prufer_edges <- function(code, K = length(code) + 2) {
  if (K == 2) return(matrix(c(1L, 2L), 1, 2))
  degree <- rep(1, K); for (v in code) degree[v] <- degree[v] + 1
  edges <- matrix(0L, K - 1, 2)
  for (t in seq_along(code)) {
    leaf <- min(which(degree == 1))
    edges[t, ] <- c(leaf, code[t])
    degree[leaf] <- 0; degree[code[t]] <- degree[code[t]] - 1
  }
  edges[K - 1, ] <- which(degree == 1)
  edges
}

# Diameter of a tree: farthest vertex from any vertex, then farthest from it.
tree_diameter <- function(edges, K) {
  adj <- lapply(1:K, function(v) c(edges[edges[, 1] == v, 2], edges[edges[, 2] == v, 1]))
  bfs <- function(s) {
    d <- rep(-1, K); d[s] <- 0; q <- s
    while (length(q)) { v <- q[1]; q <- q[-1]
      for (u in adj[[v]]) if (d[u] < 0) { d[u] <- d[v] + 1; q <- c(q, u) } }
    d
  }
  d1 <- bfs(1)
  max(bfs(which.max(d1)))
}

# BFS structure of a tree rooted at vertex 1: children, parents, edge indices.
tree_structure <- function(edges) {
  K <- nrow(edges) + 1
  adj <- lapply(1:K, function(v) {
    r1 <- which(edges[, 1] == v); r2 <- which(edges[, 2] == v)
    list(nb = c(edges[r1, 2], edges[r2, 1]), eid = c(r1, r2)) })
  child <- integer(0); parent <- integer(0); eid <- integer(0)
  seen <- 1; queue <- 1
  while (length(queue)) {
    p <- queue[1]; queue <- queue[-1]
    for (m in seq_along(adj[[p]]$nb)) { c <- adj[[p]]$nb[m]
      if (!(c %in% seen)) { seen <- c(seen, c); queue <- c(queue, c)
        child <- c(child, c); parent <- c(parent, p); eid <- c(eid, adj[[p]]$eid[m]) } }
  }
  list(K = K, child = child, parent = parent, eid = eid)
}

tree_scores_scaled <- function(st, L, a) {
  S <- matrix(0, dim(L)[1], st$K)
  for (m in seq_along(st$child)) S[, st$child[m]] <- S[, st$parent[m]] + L[, st$child[m], st$parent[m]] * a[st$eid[m]]
  S
}

# Sequential temperature scaling (preprint, Sections 5.5-5.6).
#
# Two classifiers (u, x) and (x, w) that share the class x give the partial
# model on the classes {u, x, w}: with class x as the reference,
#   log p_u / p_x = a_ux * L[u, x],   log p_w / p_x = a_wx * L[w, x].
# Keeping the scale of the first classifier fixed, the scale a_wx of the second
# is chosen to maximise the likelihood of this three-class model on the
# out-of-fold log-odds of the training points of the classes u, x, w.
fit_pair_scale <- function(Loof, y, u, x, w, a_ux) {
  ix <- which(y %in% c(u, x, w)); n <- length(ix); lab <- match(y[ix], c(u, x, w))
  su <- Loof[ix, u, x] * a_ux; sw <- Loof[ix, w, x]
  nll <- function(la) {
    S <- cbind(su, 0, sw * exp(la))
    mx <- S[cbind(seq_len(n), max.col(S, ties.method = "first"))]
    -mean(S[cbind(seq_len(n), lab)] - mx - log(rowSums(exp(S - mx))))
  }
  exp(optimize(nll, c(-8, 3))$minimum)
}

# The edge that stays unscaled: for a tree of diameter 3 the central edge (both
# end points have degree > 1); otherwise the first edge that contains a leaf.
reference_edge <- function(edges) {
  deg <- tabulate(c(edges[, 1], edges[, 2]), nbins = max(edges))
  central <- which(deg[edges[, 1]] > 1 & deg[edges[, 2]] > 1)
  if (length(central) == 1 && tree_diameter(edges, length(deg)) == 3) central
  else {
    leaf_edges <- which(deg[edges[, 1]] == 1 | deg[edges[, 2]] == 1)
    lo <- pmin(edges[leaf_edges, 1], edges[leaf_edges, 2])
    hi <- pmax(edges[leaf_edges, 1], edges[leaf_edges, 2])
    leaf_edges[order(lo, hi)[1]]
  }
}

# Scales a_e for all edges of a tree. The reference edge has a_e = 1; every
# other edge is calibrated against the edge from which it is reached, starting
# at the reference edge and moving outwards (breadth first).
fit_tree_scales_seq <- function(edges, Loof, y) {
  m <- nrow(edges); a <- rep(NA_real_, m)
  ref <- reference_edge(edges); a[ref] <- 1
  queue <- ref
  while (length(queue)) {
    e <- queue[1]; queue <- queue[-1]
    for (x in edges[e, ]) {                       # shared class x
      u <- setdiff(edges[e, ], x)
      for (f in which(is.na(a) & (edges[, 1] == x | edges[, 2] == x))) {
        w <- setdiff(edges[f, ], x)
        a[f] <- fit_pair_scale(Loof, y, u, x, w, a[e])
        queue <- c(queue, f)
      }
    }
  }
  a
}

# Predicted classes of a tree with scaled log-odds.
tree_predict_scaled <- function(edges, L, a) {
  max.col(tree_scores_scaled(tree_structure(edges), L, a), ties.method = "first")
}

# Class probabilities of the temperature-scaled minimal model (softmax of the
# scaled tree scores), for the Brier score and the cross-entropy.
tree_probs_scaled <- function(edges, L, a) {
  S <- tree_scores_scaled(tree_structure(edges), L, a)
  S <- S - apply(S, 1, max)
  P <- exp(S)
  P / rowSums(P)
}

# ---- calibrated WLW coupling (full one-vs-one model) --------------------------
# Every pair (i, j) has its own scale a_ij that multiplies both L[, i, j] and
# L[, j, i]; the scaled log-odds are coupled by the method of Wu, Lin and Weng
# (wu2_ld from decode.cpp, which must be loaded). L is an array [n, K, K];
# `pairs` is a matrix of the K(K-1)/2 pairs (i < j) in the order of `a`.
wlw_probs_scaled <- function(L, a, pairs) {
  K <- dim(L)[2]
  A <- matrix(1, K, K); A[pairs] <- a; A[pairs[, 2:1, drop = FALSE]] <- a
  t(vapply(seq_len(dim(L)[1]), function(n) wu2_ld(L[n, , ] * A), numeric(K)))
}

# Scales of all K(K-1)/2 pairs jointly minimising the cross-entropy of the WLW
# coupling on out-of-fold log-odds `Loof` [n_train, K, K] of the training
# points with labels `y`. Optimised over log a_ij (L-BFGS-B, a_ij in
# [e^-5, e^3]), starting from the uncalibrated model (all a_ij = 1).
fit_wlw_scales <- function(Loof, y) {
  K <- dim(Loof)[2]
  pairs <- which(upper.tri(diag(K)), arr.ind = TRUE)
  idx <- cbind(seq_along(y), y)
  nll <- function(la) {
    P <- wlw_probs_scaled(Loof, exp(la), pairs)
    -mean(log(pmax(P[idx], 1e-15)))
  }
  o <- optim(rep(0, nrow(pairs)), nll, method = "L-BFGS-B", lower = -5, upper = 3)
  list(a = exp(o$par), pairs = pairs, convergence = o$convergence)
}

# --- Joint calibration of a tree on its multiclass likelihood ------------------
# With the tree rooted at vertex 1 and edge e joining parent(e) and child(e),
# the class scores are S[n, k] = sum over the edges e on the path root -> k of
#   a_e * L[n, child(e), parent(e)] + b_e,
# and p = softmax(S). The scores are linear in (a, b), so the cross-entropy on
# the out-of-fold log-odds is convex. Modes: "common" (one scale for all edges,
# b = 0), "scale" (K - 1 scales, b = 0; unlike fit_tree_scales_seq no edge is
# held fixed), "affine" (K - 1 scales and K - 1 intercepts, i.e. a Platt-type
# recalibration of every edge fitted jointly).
tree_path_matrix <- function(st) {
  P <- matrix(0, st$K, st$K - 1)
  for (m in seq_along(st$child)) { P[st$child[m], ] <- P[st$parent[m], ]; P[st$child[m], st$eid[m]] <- 1 }
  P
}

tree_edge_logodds <- function(st, L) {       # G[n, e] = L[n, child(e), parent(e)]
  G <- matrix(0, dim(L)[1], st$K - 1)
  for (m in seq_along(st$child)) G[, st$eid[m]] <- L[, st$child[m], st$parent[m]]
  G
}

softmax_rows <- function(S) { S <- S - apply(S, 1, max); E <- exp(S); E / rowSums(E) }

fit_tree_joint <- function(edges, Loof, y, mode = c("affine", "scale", "common")) {
  mode <- match.arg(mode)
  st <- tree_structure(edges); P <- tree_path_matrix(st); G <- tree_edge_logodds(st, Loof)
  E <- st$K - 1; n <- length(y); Y <- matrix(0, n, st$K); Y[cbind(seq_len(n), y)] <- 1
  unpack <- function(th) switch(mode,
    common = list(a = rep(th[1], E), b = rep(0, E)),
    scale  = list(a = th, b = rep(0, E)),
    affine = list(a = th[1:E], b = th[E + 1:E]))
  f <- function(th) {
    ab <- unpack(th); S <- (G * rep(ab$a, each = n) + rep(ab$b, each = n)) %*% t(P)
    mx <- apply(S, 1, max)
    -mean(rowSums(S * Y) - mx - log(rowSums(exp(S - mx))))
  }
  g <- function(th) {
    ab <- unpack(th); S <- (G * rep(ab$a, each = n) + rep(ab$b, each = n)) %*% t(P)
    RP <- (softmax_rows(S) - Y) %*% P / n             # n x E
    ga <- colSums(RP * G); gb <- colSums(RP)
    switch(mode, common = sum(ga), scale = ga, affine = c(ga, gb))
  }
  th0 <- switch(mode, common = 1, scale = rep(1, E), affine = c(rep(1, E), rep(0, E)))
  o <- optim(th0, f, g, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  c(unpack(o$par), list(edges = edges, mode = mode, nll = o$value, convergence = o$convergence))
}

tree_probs_joint <- function(fit, L) {
  st <- tree_structure(fit$edges); P <- tree_path_matrix(st); G <- tree_edge_logodds(st, L)
  n <- dim(L)[1]
  softmax_rows((G * rep(fit$a, each = n) + rep(fit$b, each = n)) %*% t(P))
}
