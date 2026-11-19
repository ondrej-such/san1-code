# Section 6: Algorithms A-G and the full WLW model with RBF-SVM pairwise
# classifiers on satimage, segment, usps, mnist and letter (800/1000, runs 0-19).
#
# Per dataset and run:
#   1. full pairwise fit (log_exp1_dfs_svm: per-pair (C, gamma) and Platt
#      sigmoid by 5-fold CV log-loss; Cohen's d on the CV log-odds);
#   2. out-of-fold log-odds of the training points for 7 repeats of 10-fold CV
#      (fold refits reuse (C, gamma) and the sigmoid, as in eval3); the same
#      folds are shared by all algorithms, so their comparison is paired;
#   3. every algorithm proposes its candidate trees; per repeat and metric the
#      candidate with the smallest CV score (summed over the folds, i.e. scored
#      on the pooled out-of-fold log-odds) is selected and scored on the test set.
# Step 3 is equivalent to cvf() in exp4.R, but the CV refits are done once per
# run instead of once per algorithm.
#
# Cache (for recalibration without new hyperparameter selection):
#   cache/svm_alg3/<dataset>_run_XX.rds: le (K, dM, pred1, pred2, GT1, GT2,
#   L = C, gamma, Platt a, b; cvdv = the 5-fold CV decision values the sigmoids
#   were fitted on), Loof [repeat, n_train, K, K], fold [repeat, n_train].
#   The fitted SVMs are not stored (mnist: ~30 MB per run); they are refitted
#   from L without any selection.
# Results: data/alg3svm-<dataset>-7.csv (columns as written by eval3_all).
suppressPackageStartupMessages({ library(furrr); library(dplyr); library(purrr) })

R_REP <- 7; FOLDS <- 10
CACHE_DIR <- Sys.getenv("SVM_ALG3_CACHE", "cache/svm_alg3")
dss <- c("mnist", "usps", "letter", "satimage", "segment")   # heaviest first

job <- function(ds, run) {
  dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  res_file <- sprintf("%s/%s_run_%02d_res.rds", CACHE_DIR, ds, run)
  if (file.exists(res_file)) return(readRDS(res_file))
  # exp4.R -> msvm.R compiles decode.cpp into the shared cacheDir; retry on collisions
  for (k in 1:5) {
    if (exists("algs3") && exists("wu2_ld")) break
    ok <- tryCatch({ suppressPackageStartupMessages(source("exp4.R")); TRUE }, error = function(e) FALSE)
    if (!ok) Sys.sleep(runif(1, 2, 8))
  }
  seed <- 1000 * match(ds, dss) + run
  fit_file <- sprintf("%s/%s_run_%02d.rds", CACHE_DIR, ds, run)
  if (file.exists(fit_file)) {
    cc <- readRDS(fit_file)
  } else {
    set.seed(seed)
    t0 <- proc.time()[["elapsed"]]
    invisible(capture.output(dfs <- read_wlws(800, ds, run)))
    le <- log_exp1_dfs(dfs, method = "svm")
    GT <- le$GT1; n <- length(GT); K <- le$K
    Loof <- array(0, c(R_REP, n, K, K)); fold <- matrix(0L, R_REP, n)
    for (r in 1:R_REP) {
      fold[r, ] <- make5folds(GT, FOLDS)
      for (f in 1:FOLDS) {
        d2 <- dfs; d2$train <- dfs$train[fold[r, ] != f, ]; d2$test <- dfs$train[fold[r, ] == f, ]
        Loof[r, fold[r, ] == f, , ] <- aperm(log_exp1_dfs(d2, L = le$L, method = "svm")$pred2, c(3, 1, 2))
      }
    }
    cc <- list(dataset = ds, run = run, le = le[c("K", "dM", "pred1", "pred2", "GT1", "GT2", "L", "cvdv")],
               Loof = Loof, fold = fold, minutes = (proc.time()[["elapsed"]] - t0) / 60)
    saveRDS(cc, fit_file)
  }
  le <- cc$le; y <- as.integer(le$GT1); yt <- le$GT2; n <- length(y)
  mets <- c("error_rate", "brier_score", "cross_entropy")
  set.seed(seed + 500000)                      # random algorithms D, E and ties
  rows <- list()
  for (alg in c("A", "B", "C", "D", "E", "F", "G")) {
    trees <- algs3[[alg]](le)
    for (r in 1:R_REP) {
      Lr <- aperm(cc$Loof[r, , , ], c(2, 3, 1))           # [K, K, n]
      cv <- vapply(trees, function(t) {
        P <- log_coupling(t, Lr); vapply(mets, function(m) metrics[[m]](P, y), numeric(1))
      }, numeric(length(mets)))
      cv <- matrix(cv, nrow = length(mets), dimnames = list(mets, NULL))
      for (m in mets) {
        imin <- sample(as.list(which(cv[m, ] == min(cv[m, ]))), 1)[[1]]
        tree <- trees[[imin]]
        rows[[length(rows) + 1]] <- data.frame(
          rep = r, tree = str_prufer(tree), metric = m, cv = FOLDS, repeats = R_REP,
          train_score = cv[m, imin] / n,
          test_score = metrics[[m]](log_coupling(tree, le$pred2), yt) / length(yt),
          alg = alg, ncand = length(trees))
      }
    }
  }
  W <- t(vapply(seq_along(yt), function(i) wu2_ld(le$pred2[, , i]), numeric(le$K)))
  for (m in mets) rows[[length(rows) + 1]] <- data.frame(
    rep = NA_integer_, tree = NA_character_, metric = m, cv = 0, repeats = NA_integer_,
    train_score = NA_real_, test_score = metrics[[m]](W, yt) / length(yt), alg = "WLW", ncand = NA_integer_)
  res <- bind_rows(rows) |> mutate(dataset = ds, run = run, fit_minutes = cc$minutes)
  saveRDS(res, res_file)
  res
}

if (sys.nframe() == 0) {
  args <- commandArgs(TRUE)
  workers <- if (length(args) >= 1) as.integer(args[1]) else 10
  if (length(args) >= 2) dss_run <- strsplit(args[2], ",")[[1]] else dss_run <- dss
  dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create("data", showWarnings = FALSE)
  jobs <- expand.grid(run = 0:19, ds = dss_run, stringsAsFactors = FALSE)
  plan(multisession, workers = workers)
  t0 <- Sys.time()
  res <- future_map2(jobs$ds, jobs$run, job,
                     .options = furrr_options(seed = TRUE, scheduling = Inf)) |> bind_rows()
  for (ds in unique(res$dataset))
    write.csv(res |> filter(dataset == ds), sprintf("data/alg3svm-%s-7.csv", ds), row.names = FALSE)
  message("ALL DONE in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
