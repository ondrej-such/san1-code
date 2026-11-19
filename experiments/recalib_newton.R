# Fast fitting of the joint per-edge temperature scaling ("scale" mode of
# fit_tree_joint in temperature_scaling.R) by Newton's method.
#
# For a tree with path matrix P [K, K-1] and edge log-odds G [n, K-1] the
# calibrated scores are S = (G diag(a)) P^T and p = softmax(S). The
# cross-entropy is convex in the K-1 scales a (a conditional logit model with
# features X[n, k, e] = G[n, e] P[k, e]), with
#   gradient  g = (1/n) sum_n G_n o (E_p[P] - P_{y_n})
#   Hessian   H = (1/n) sum_n Cov_p(G_n o P_k),
# so Newton's method with a backtracking line search converges in a few steps.
# Same optimum as fit_tree_joint(..., mode = "scale"), about 10x faster.

softmax_fast <- function(S) {
  m <- S[, 1]
  for (k in seq_len(ncol(S))[-1]) m <- pmax(m, S[, k])
  E <- exp(S - m)
  E / rowSums(E)
}

# mean cross-entropy of the scores S [n, K] for labels y
xent_scores <- function(S, y) {
  m <- S[, 1]
  for (k in seq_len(ncol(S))[-1]) m <- pmax(m, S[, k])
  mean(m + log(rowSums(exp(S - m))) - S[cbind(seq_along(y), y)])
}

fit_scale_newton <- function(G, P, y, a0 = rep(1, ncol(G)), maxit = 50, tol = 1e-9) {
  n <- nrow(G); K <- nrow(P); E <- ncol(G)
  Py <- P[y, , drop = FALSE]                              # P_{y_n}, n x E
  tP <- t(P)
  # Hessian, first term: sum_n G_ne G_nf sum_k p_nk P_ke P_kf, all (e, f) at once
  ie <- rep(seq_len(E), E); jf <- rep(seq_len(E), each = E)
  PP <- P[, ie, drop = FALSE] * P[, jf, drop = FALSE]     # K x E^2
  GG <- G[, ie, drop = FALSE] * G[, jf, drop = FALSE]     # n x E^2
  obj <- function(a) xent_scores((G * rep(a, each = n)) %*% tP, y)
  a <- a0; f0 <- obj(a)
  for (it in seq_len(maxit)) {
    Pr <- softmax_fast((G * rep(a, each = n)) %*% tP)
    M <- Pr %*% P                                         # E_p[P_k] per point, n x E
    grad <- colSums(G * (M - Py)) / n
    if (max(abs(grad)) < tol) break
    H <- (matrix(colSums(GG * (Pr %*% PP)), E, E) - crossprod(G * M)) / n
    step <- tryCatch(solve(H + diag(1e-10, ncol(G)), grad), error = function(e) grad)
    t <- 1
    repeat {
      a1 <- a - t * step; f1 <- obj(a1)
      if (is.finite(f1) && f1 <= f0 - 1e-4 * t * sum(grad * step)) break
      t <- t / 2
      if (t < 1e-10) { a1 <- a; f1 <- f0; break }
    }
    done <- abs(f0 - f1) < 1e-14
    a <- a1; f0 <- f1
    if (done) break
  }
  a
}

# probabilities of a tree with scales a on edge log-odds G
probs_scaled <- function(G, P, a) softmax_fast((G * rep(a, each = nrow(G))) %*% t(P))
