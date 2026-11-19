# DAGSVM baseline (Platt, Cristianini & Shawe-Taylor 2000) for the setup of
# Section 6: the same pairwise RBF-SVMs, evaluated on the test set from the
# stored test log-odds (cache of eval_svm_alg3.R), so nothing is refitted.
#
# The DAG keeps a list of candidate classes; at each of its K - 1 steps it
# evaluates the classifier of the first and the last class of the list and
# removes the loser. Every test point thus uses K - 1 pairwise classifiers, as a
# minimal model does, but which ones depends on the point.
#
# Probabilities: the K - 1 comparisons of a test point form a spanning tree
# (every class except the winner is removed once, by a class still in the list),
# so the minimal-model coupling on this point-specific tree gives probabilities:
# log p_l = log p_w + r_lw for the loser l of a comparison with w, the winner of
# the DAG at the root. The class with the largest probability is the DAG winner,
# so the test error is unchanged.
#
#   DAG   list in the natural class order 1..K (as in Platt et al.), with the
#         probabilities above
#   DAGr  mean test error over N_ORD random list orders (sensitivity to order)
#   VOTE  max-wins voting over all K(K-1)/2 classifiers (reference)
# Output: data/dagsvm-svm.csv (dataset, run, method, metric, score per test point)
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
CACHE_DIR <- Sys.getenv("SVM_ALG3_CACHE", "cache/svm_alg3")
dss <- c("satimage", "segment", "usps", "mnist", "letter")
N_ORD <- 50

# P [K, K, n]: P[i, j, x] = log-odds of class i against class j (antisymmetric)
dag_predict <- function(P, ord, probs = FALSE) {
  K <- length(ord); n <- dim(P)[3]; lo <- rep(1L, n); hi <- rep(K, n)
  los <- win <- matrix(0L, n, K - 1)
  for (s in seq_len(K - 1)) {
    i <- ord[lo]; j <- ord[hi]
    iwins <- P[cbind(i, j, seq_len(n))] > 0
    los[, s] <- ifelse(iwins, j, i); win[, s] <- ifelse(iwins, i, j)
    hi <- hi - iwins; lo <- lo + !iwins
  }
  pred <- ord[lo]
  if (!probs) return(pred)
  L <- matrix(0, n, K); x <- seq_len(n)
  for (s in rev(seq_len(K - 1)))   # the winner of step s is the DAG winner or removed later
    L[cbind(x, los[, s])] <- L[cbind(x, win[, s])] + P[cbind(los[, s], win[, s], x)]
  E <- exp(L - apply(L, 1, max))
  list(pred = pred, probs = E / rowSums(E))
}
vote_predict <- function(P) {
  n <- dim(P)[3]
  W <- t(vapply(seq_len(n), function(x) rowSums(P[, , x] > 0), numeric(dim(P)[1])))
  max.col(W, ties.method = "first")
}
brier <- function(p, y) { T <- matrix(0, nrow(p), ncol(p)); T[cbind(seq_along(y), y)] <- 1; mean(rowSums((p - T)^2)) }
xent <- function(p, y) mean(-log(pmax(p[cbind(seq_along(y), y)], 1e-15)))   # as cross_entropy in exp4.R

if (sys.nframe() == 0) {
  set.seed(20261118)
  rows <- list()
  for (ds in dss) for (run in 0:19) {
    cc <- readRDS(sprintf("%s/%s_run_%02d.rds", CACHE_DIR, ds, run))
    P <- cc$le$pred2; y <- cc$le$GT2; K <- cc$le$K
    stopifnot(max(abs(P + aperm(P, c(2, 1, 3)))) < 1e-12)
    d <- dag_predict(P, seq_len(K), probs = TRUE)
    stopifnot(all(max.col(d$probs, ties.method = "first") == d$pred))
    e_rnd <- mean(vapply(seq_len(N_ORD), function(r) mean(dag_predict(P, sample(K)) != y), numeric(1)))
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, run = run,
      method = c("DAG", "DAG", "DAG", "DAGr", "VOTE"),
      metric = c("error_rate", "brier_score", "cross_entropy", "error_rate", "error_rate"),
      score = c(mean(d$pred != y), brier(d$probs, y), xent(d$probs, y), e_rnd, mean(vote_predict(P) != y)))
  }
  res <- bind_rows(rows)
  write.csv(res, "data/dagsvm-svm.csv", row.names = FALSE)

  # comparison with WLW and the minimal models of Sections 6 (U) and 7 (R)
  s6 <- purrr::map_dfr(sprintf("data/alg3svm-%s-7.csv", dss), read.csv) |>
    group_by(dataset, metric, run, alg) |> summarise(score = mean(test_score), .groups = "drop")
  fr <- sprintf("data/recal-svm-%s-7.csv", dss); fr <- fr[file.exists(fr)]
  s7 <- purrr::map_dfr(fr, read.csv) |> filter(method == "R", alg != "WLW") |>
    group_by(dataset, metric, run, alg) |> summarise(score = mean(test_score), .groups = "drop") |>
    mutate(alg = paste0(alg, "_R"))
  all <- bind_rows(res |> rename(alg = method), s6, s7)
  options(width = 200)
  for (m in c("error_rate", "brier_score", "cross_entropy")) {
    tab <- all |> filter(metric == m, alg %in% c("WLW", "VOTE", "DAG", "DAGr", "C", "G", "C_R", "G_R")) |>
      group_by(dataset, alg) |> summarise(v = mean(score) * ifelse(m == "error_rate", 100, 1), .groups = "drop") |>
      pivot_wider(names_from = alg, values_from = v) |> arrange(match(dataset, dss))
    cat("\nmean", m, if (m == "error_rate") "(%)", "over 20 runs\n"); print(as.data.frame(tab), digits = 3)
  }
  cmp <- function(a, b, m) all |> filter(metric == m, alg %in% c(a, b)) |>
    pivot_wider(names_from = alg, values_from = score) |>
    filter(!is.na(.data[[a]]), !is.na(.data[[b]])) |> group_by(dataset) |>
    summarise(diff = mean(.data[[a]] - .data[[b]]) * ifelse(m == "error_rate", 100, 1),
              a_wins = sum(.data[[a]] < .data[[b]]), ties = sum(.data[[a]] == .data[[b]]),
              p = t.test(.data[[a]], .data[[b]], paired = TRUE)$p.value, .groups = "drop") |>
    mutate(cmp = paste(a, "vs", b), metric = m)
  out <- bind_rows(lapply(c("error_rate", "brier_score", "cross_entropy"), function(m)
    bind_rows(cmp("DAG", "WLW", m), cmp("DAG", "C_R", m), cmp("DAG", "G_R", m)))) |>
    arrange(metric, cmp, match(dataset, dss))
  cat("\npaired comparisons over the 20 runs (diff: test error in percentage points; a_wins = runs in which DAG is better)\n")
  print(as.data.frame(out), digits = 3)
  cat("\nDAG vs DAGr (order sensitivity, test error)\n")
  print(as.data.frame(cmp("DAG", "DAGr", "error_rate")), digits = 3)
}
