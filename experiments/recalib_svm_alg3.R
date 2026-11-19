# Recalibration of minimal models (RBF-SVM, Section 6 setup). Uses the model
# cache of eval_svm_alg3.R: no SVM is fitted and no hyperparameter is selected.
#
# For every candidate tree of Algorithms A-G and every CV repetition r:
#   * the tree is recalibrated by joint per-edge temperature scaling (K - 1
#     scales a_e multiplying the log-odds of the edges, fitted by maximum
#     likelihood = minimal cross-entropy; recalib_newton.R);
#   * its CV score after recalibration is cross-fitted: in each of the 10 folds
#     of repetition r the scales are fitted on the out-of-fold log-odds of the
#     other 9 folds and the held-out fold is scored with them;
#   * every algorithm selects, per metric, the candidate with the smallest
#     recalibrated CV score; its scales are then fitted on all out-of-fold
#     log-odds of repetition r and applied to the test log-odds.
# The uncalibrated selection of eval_svm_alg3.R is replayed with the same random
# numbers (same candidate trees of the random Algorithms D and E, same ties), so
# method "U" reproduces data/alg3svm-<dataset>-7.csv exactly (a check).
#
# Output: cache/svm_recal/<dataset>_run_XX.rds and data/recal-svm-<dataset>-7.csv
#   columns rep, metric, alg, method (U = uncalibrated, R = recalibrated), tree,
#   cv_score, test_score, dataset, run; WLW (uncalibrated) as in eval_svm_alg3.R.
suppressPackageStartupMessages({ library(furrr); library(dplyr); library(purrr) })

R_REP <- 7; FOLDS <- 10
CACHE_DIR <- Sys.getenv("SVM_ALG3_CACHE", "cache/svm_alg3")
OUT_DIR <- "cache/svm_recal"
dss <- c("mnist", "usps", "letter", "satimage", "segment")   # as in eval_svm_alg3.R (seeds)
run_order <- c("letter", "mnist", "usps", "segment", "satimage")   # heaviest first

job <- function(ds, run) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  out <- sprintf("%s/%s_run_%02d.rds", OUT_DIR, ds, run)
  if (file.exists(out)) return(readRDS(out))
  for (k in 1:5) {
    if (exists("algs3") && exists("wu2_ld") && exists("fit_scale_newton")) break
    ok <- tryCatch({
      suppressPackageStartupMessages(source("exp4.R")); source("recalib_newton.R"); TRUE
    }, error = function(e) FALSE)
    if (!ok) Sys.sleep(runif(1, 2, 8))
  }
  t0 <- proc.time()[["elapsed"]]
  cc <- readRDS(sprintf("%s/%s_run_%02d.rds", CACHE_DIR, ds, run))
  le <- cc$le; K <- le$K; y <- as.integer(le$GT1); yt <- le$GT2; n <- length(y)
  L2 <- aperm(le$pred2, c(3, 1, 2))                        # test log-odds [n2, K, K]
  mets <- c("error_rate", "brier_score", "cross_entropy")
  seed <- 1000 * match(ds, dss) + run
  set.seed(seed + 500000)                                  # exactly as in eval_svm_alg3.R
  rows <- list()
  for (alg in c("A", "B", "C", "D", "E", "F", "G")) {
    trees <- algs3[[alg]](le)
    st <- lapply(trees, function(t) tree_structure(prufer_edges(as.integer(t), K)))
    Pm <- lapply(st, tree_path_matrix)
    for (r in 1:R_REP) {
      Lr <- aperm(cc$Loof[r, , , ], c(2, 3, 1))           # [K, K, n] for log_coupling
      Lr3 <- cc$Loof[r, , , ]                              # [n, K, K]
      fold <- cc$fold[r, ]
      # uncalibrated CV scores (as in eval_svm_alg3.R)
      cvU <- vapply(trees, function(t) {
        P <- log_coupling(t, Lr); vapply(mets, function(m) metrics[[m]](P, y), numeric(1))
      }, numeric(length(mets)))
      cvU <- matrix(cvU, nrow = length(mets), dimnames = list(mets, NULL))
      # recalibrated, cross-fitted CV scores, and the scales fitted on all points
      a_full <- vector("list", length(trees))
      cvR <- vapply(seq_along(trees), function(i) {
        G <- tree_edge_logodds(st[[i]], Lr3)
        a <- fit_scale_newton(G, Pm[[i]], y); a_full[[i]] <<- a
        s <- numeric(length(mets)); names(s) <- mets
        for (f in seq_len(FOLDS)) {
          tr <- fold != f
          af <- fit_scale_newton(G[tr, , drop = FALSE], Pm[[i]], y[tr], a0 = a)
          Pf <- probs_scaled(G[!tr, , drop = FALSE], Pm[[i]], af)
          for (m in mets) s[m] <- s[m] + metrics[[m]](Pf, y[!tr])
        }
        s
      }, numeric(length(mets)))
      cvR <- matrix(cvR, nrow = length(mets), dimnames = list(mets, NULL))
      for (m in mets) {
        imin <- sample(as.list(which(cvU[m, ] == min(cvU[m, ]))), 1)[[1]]   # same RNG use as eval_svm_alg3.R
        rows[[length(rows) + 1]] <- data.frame(
          rep = r, metric = m, alg = alg, method = "U", tree = str_prufer(trees[[imin]]),
          cv_score = cvU[m, imin] / n,
          test_score = metrics[[m]](log_coupling(trees[[imin]], le$pred2), yt) / length(yt))
        j <- which.min(cvR[m, ])
        G2 <- tree_edge_logodds(st[[j]], L2)
        rows[[length(rows) + 1]] <- data.frame(
          rep = r, metric = m, alg = alg, method = "R", tree = str_prufer(trees[[j]]),
          cv_score = cvR[m, j] / n,
          test_score = metrics[[m]](probs_scaled(G2, Pm[[j]], a_full[[j]]), yt) / length(yt))
      }
    }
  }
  W <- t(vapply(seq_along(yt), function(i) wu2_ld(le$pred2[, , i]), numeric(K)))
  for (m in mets) rows[[length(rows) + 1]] <- data.frame(
    rep = NA_integer_, metric = m, alg = "WLW", method = "U", tree = NA_character_,
    cv_score = NA_real_, test_score = metrics[[m]](W, yt) / length(yt))
  res <- bind_rows(rows) |> mutate(dataset = ds, run = run, minutes = (proc.time()[["elapsed"]] - t0) / 60)
  saveRDS(res, out)
  res
}

if (sys.nframe() == 0) {
  args <- commandArgs(TRUE)
  workers <- if (length(args) >= 1) as.integer(args[1]) else 10
  dss_run <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else run_order
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  jobs <- expand.grid(run = 0:19, ds = dss_run, stringsAsFactors = FALSE)
  plan(multisession, workers = workers)
  t0 <- Sys.time()
  res <- future_map2(jobs$ds, jobs$run, job,
                     .options = furrr_options(seed = TRUE, scheduling = Inf)) |> bind_rows()
  for (ds in unique(res$dataset))
    write.csv(res |> filter(dataset == ds), sprintf("data/recal-svm-%s-7.csv", ds), row.names = FALSE)
  message("ALL DONE in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
