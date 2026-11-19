
suppressPackageStartupMessages(library(igraph))
source("exp1_log.R")

prufer_code <- function(g) {
  if (!is_connected(g) || gsize(g) != vcount(g) - 1) {
    stop("Graph must be a tree (connected and acyclic)")
  }
  
  # Ensure vertex names are numeric and sorted
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(V(g))
  }
  node_labels <- as.integer(V(g)$name)
  names(node_labels) <- V(g)
  
  prufer <- c()
  g_copy <- g
  
  while (vcount(g_copy) > 2) {
    degrees <- degree(g_copy)
    
    # Get the leaf with the smallest label
    leaves <- which(degrees == 1)
    leaf_names <- V(g_copy)[leaves]$name
    smallest_leaf <- as.character(min(as.integer(leaf_names)))
    
    #print(smallest_leaf)
    # Get the neighbor of the smallest leaf
    neighbor <- neighbors(g_copy, smallest_leaf)[[1]]$name
    
    # Append neighbor to Prufer code
    prufer <- c(prufer, as.integer(neighbor))
    
    # Remove the leaf
    g_copy <- delete_vertices(g_copy, smallest_leaf)
    #print(g_copy)
  }
  
  return(prufer)
}

# Get the minimum spanning tree based on Cohen's d
# 
get_A <- function(le) {
  dM <- le$dM
  g <- igraph::graph_from_edgelist(dM[,1:2], directed = FALSE)
  mst <- igraph::mst(g, weights = dM[,3])
  list(prufer_code(mst))
}

# Get the minimum spanning tree based on CV accuracy
# 
get_AC <- function(le) {
  dM <- le$dM
  g <- igraph::graph_from_edgelist(dM[,1:2], directed = FALSE)
  mst <- igraph::mst(g, weights = dM[,4])
  list(prufer_code(mst))
}

# Enumarate all stars
#
get_B <- function(le) {
 K <- le$K
 lapply (1:K, function(i) {
    rep(i, K -2)
  })
}

# Select K choose 2 trees of diameter 3
#
get_D <- function(le) {
  K <- le$K
  pairs <- expand.grid(i = 1:K, j = 1:K) |> filter(i < j)
  lapply(1:nrow(pairs), function(row) {
      i <- pairs$i[row]
      j <- pairs$j[row]
      rest <- (1:K)[-c(i,j) ]
      edges <- cbind(rest, if_else (rbinom (K-2, 1, 0.5) == 0, i, j))
      edges<- rbind(edges, c(i,j))
      g <- igraph::graph_from_edgelist(edges, directed = FALSE)
      prufer_code(g)
    }
  ) 
}

# Choose all trees of diameter 3 so that non-center vertices are connected to
# the closest center vertex
#
get_C <- function(le) {
  K <- le$K
  pairs <- expand.grid(i = 1:K, j = 1:K) |> filter(i < j)
  dS <- matrix(0, nrow = K, ncol = K)
  for ( i in 1:nrow(le$dM)) {
    row <- le$dM[i,]
    #print(row)
    i <- row["i"]
    j <- row["j"]
    d <- row["cohend"]
    dS[i,j] = d
    dS[j,i] = d
  }
  
  
  lapply(1:nrow(pairs), function(row) {
    i <- pairs$i[row]
    j <- pairs$j[row]
    rest <- (1:K)[-c(i,j) ]
    di <- dS[i,rest]
    dj <- dS[j,rest]
    edges <- cbind(rest, if_else (di < dj, i, j))
    edges<- rbind(edges, c(i,j))
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    prufer_code(g)
  }
  ) 
}

# Select (K choose 2) random trees
#
get_E <- function(le) {
  K <- le$K
  lapply(1:((K-1)* K/2), function (i) sample(1:K, K-2, replace = TRUE))
}


# Broom tree
get_F <- function(le) {
  K <- le$K
  dS <- matrix(0, nrow = K, ncol = K)

  for (i in 1:nrow(le$dM)) {
    row <- le$dM[i, ]
    ii <- row["i"]
    jj <- row["j"]
    d <- row["cohend"]
    dS[ii, jj] <- d
    dS[jj, ii] <- d
  }
  
  broom_trees <- lapply(1:K, function(center) {
    rest <- setdiff(1:K, center)
    max_idx <- which.max(dS[center, rest])
    remaining <- rest[max_idx]
    leaves <- rest[-max_idx]
    leaf_dists <- dS[remaining, leaves]
    attach_leaf <- leaves[which.min(leaf_dists)]
    edges <- rbind(
      cbind(leaves, center),
      c(remaining, attach_leaf)
    )
    
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    prufer_code(g)
  })
}

# Brooms with targeted handle
#
get_G <- function(le) {
  K <- le$K
  dS <- matrix(0, nrow = K, ncol = K)
  
  for (i in 1:nrow(le$dM)) {
    row <- le$dM[i, ]
    ii <- row["i"]
    jj <- row["j"]
    d <- row["cohend"]
    dS[ii, jj] <- d
    dS[jj, ii] <- d
  }
  pairs <- expand.grid(i = 1:K, j = 1:K) |> dplyr::filter(i != j)
  lapply(1:nrow(pairs), function(r) {
    i <- pairs$i[r]
    j <- pairs$j[r]
    rest = (1:K)[-c(i,j)]
    ds <- dS[j, rest]
    n <- which.min(ds)
    v = rep(i, K-2)
    v[n] = j
    v
    })
}

# Enumerate all brooms
#
get_BR <- function(le) {
  K <- le$K
  triples <- expand.grid(i = 1:K, j = 1:K,  k = 1:(K-2)) |>
            filter(i != j) 
  
  sapply (1:nrow(triples), function(r) {
    i <- triples$i[r]
    j <- triples$j[r]
    k <- triples$k[r]
    v <- rep(i, K - 2)
    v[k] = j
    #print(v)
    #print(class(v))
    v
  }) |> t() |> apply(1, identity, simplify = F)
}

# Enumarate all trees of diameter 3
#
get_D3 <- function(le) {
  K <- le$K
  num_columns <- K-2
  pairs <- expand.grid(i = 1:K, j = 1:K) |> filter(i < j)
  df <- map(1:nrow(pairs), function (k) {
    i <- pairs$i[k]
    j <- pairs$j[k]
    base_vector <- list(c(i,j ))
    input_list <- rep(base_vector, num_columns)
    
    # 4. Use expand.grid() to generate all combinations
    result_df <- expand.grid(input_list)
    stopifnot(sum(result_df[1,]) == i *num_columns)
    stopifnot(sum(result_df[nrow(result_df),]) == j *num_columns)
    result_df
  }) |> list_rbind()
  df[-c(1, nrow(df)), ] |> apply(1, identity, simplify = F)
}

# Enumarate all trees
#
get_TR <- function(le) {
  K <- le$K
  gr <- expand.grid(rep(list(1:K), K-2)) |>
    apply(1, identity, simplify = F)
}

# Creates stratified f=5 folds by a given group. 
#
make5folds <- function(grp, f = 5) {
  gs <- unique(grp)
  K <- length(gs)
  res <- vector("integer", length(grp))
  for (g in gs) {
    N <- sum(grp == g) 
    N5 <- N %/% f
    ubds <- N5 * seq(1, f)
    R5 <- N %% f
    v <- as.integer(rep(0, f))
    if (R5 > 0) {
      pos <- sample(1:f, size = R5, replace = F)
      v[pos] <- 1
    }
    cumv <- cumsum(v)
    cubds <- ubds + cumv
    idx <- grp == g
    perm <- sample(1:N, replace = F, size = N)  - 1
    folds <- sapply(perm, function(x) sum(x >= cubds)) + 1
    res[idx] <- folds
  }
  res
}

eval_tree <- function(tree, pred, GT) {
  soft <- log_coupling(tree, pred)
  hard <- apply(soft, 1, which.max)
  sum(hard == GT)
}

brier_score <- function(probs, GT) {
  K <- ncol(probs)
  n <- nrow(probs)
  targets <- matrix(0, n, K)
  targets[cbind(1:n, GT)] <- 1
  sum(rowSums((probs - targets)^2))
}

cross_entropy <- function(probs, GT) {
  n <- nrow(probs)
  p_true <- probs[cbind(1:n, GT)]
  -sum(log(pmax(p_true, 1e-15)))
}

error_rate <- function(probs, GT) {
  # max.col(ties.method = "first") is the vectorised which.max of every row
  sum(max.col(probs, ties.method = "first") != GT)
}

metrics <- new.env(parent = emptyenv())
metrics[["brier_score"]]  = brier_score
metrics[["cross_entropy"]] = cross_entropy
metrics[["error_rate"]] = error_rate


source("temperature_scaling.R")

# ---- timing of all actions ---------------------------------------------------
# Every worker process keeps its own log in .timing_env. timed(stage, expr)
# evaluates expr, records its wall-clock time with the current algorithm
# (.timing_env$alg) and returns the value. eval3_all() writes the log of each
# run to data/timing-<file>-<dataset>-<repeats>.csv.
.timing_env <- new.env()
.timing_env$log <- list()
.timing_env$alg <- NA_character_

timed <- function(stage, expr, detail = NA_character_) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  .timing_env$log[[length(.timing_env$log) + 1]] <- data.frame(
    stage = stage, alg = .timing_env$alg, detail = detail,
    seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
  value
}

timing_table <- function() {
  if (length(.timing_env$log) == 0) return(NULL)
  do.call(rbind, .timing_env$log)
}

# Test score of the temperature-scaled minimal model of `tree` (a Pruefer code).
# The scales are fitted sequentially on the out-of-fold log-odds `Loof`
# (array [n_train, K, K]) of the training points `y_train`; see
# temperature_scaling.R. Returns the metric per test point, like test_score.
calibrated_test_score <- function(tree, met, pred_test, GT_test, Loof, y_train) {
  e <- prufer_edges(as.integer(tree), dim(pred_test)[1])
  a <- fit_tree_scales_seq(e, Loof, as.integer(y_train))
  P <- tree_probs_scaled(e, aperm(pred_test, c(3, 1, 2)), a)
  metrics[[met]](P, GT_test) / length(GT_test)
}

# With calibrate = TRUE every selected tree is also scored after temperature
# scaling (column test_score_cal, NA otherwise). The out-of-fold log-odds of the
# first repeat are returned as attribute "Loof" for the non-CV baseline.
cvf <- function(trees, dfs, le, folds, repeats, method, calibrate = FALSE) {
  Loof_first <- NULL
  fold_first <- NULL   # fold assignment of the first repeat (attribute "fold")
  res <- map(1:repeats, function(j) {
    fold <- make5folds(dfs$train$class_id, f = folds)
    if (j == 1) fold_first <<- fold

    scores <- matrix(0, nrow = length(ls(metrics)), ncol = length(trees))
    rownames(scores) <- ls(metrics)
    Loof <- if (calibrate) array(0, c(nrow(dfs$train), le$K, le$K))

    for (i in 1:folds) {
      dfs2 <- dfs
      dfs2$train = dfs$train[fold != i,]
      dfs2$test = dfs$train[fold == i, ]

      le2 <- timed("cv_fold_fit", log_exp1_dfs(dfs2, L = le$L, method = method),
                   detail = sprintf("rep%d fold%d", j, i))
      if (calibrate) Loof[fold == i, , ] <- aperm(le2$pred2, c(3, 1, 2))

      # score every candidate tree in this fold, per metric
      timed("cv_tree_scoring", {
        for (t in seq_along(trees)) {
          probs <- log_coupling(trees[[t]], le2$pred2)
          for (met in ls(metrics)) {
            scores[met, t] <- scores[met, t] + metrics[[met]](probs, le2$GT2)
          }
        }
      }, detail = sprintf("rep%d fold%d (%d trees)", j, i, length(trees)))
    }

    if (calibrate && j == 1) Loof_first <<- Loof

    map(ls(metrics), function(met) {
      # brier_score / cross_entropy / error_rate are all losses -> pick the minimum
      imin <- sample(as.list(which(scores[met, ] == min(scores[met, ]))), 1)[[1]]
      tree <- trees[[imin]]
      test_probs <- log_coupling(tree, le$pred2)

      data.frame(rep = j, tree = str_prufer(tree), metric = met, cv = folds,
                 repeats = repeats,
                 train_score = scores[met, imin] / nrow(dfs$train),
                 test_score = metrics[[met]](test_probs, le$GT2) / length(le$GT2),
                 test_score_cal = if (calibrate)
                   timed("calibration",
                         calibrated_test_score(tree, met, le$pred2, le$GT2, Loof, dfs$train$class_id),
                         detail = sprintf("cv rep%d %s", j, met))
                 else NA_real_)
    }) |> list_rbind()
  }) |> list_rbind()

  attr(res, "Loof") <- Loof_first
  attr(res, "fold") <- fold_first
  return(res)
  

  # if(length(trees) == 1) {
  #   res_a <- matrix(res$correct, nrow = 1)
  #   res_b <- matrix(res$brier, nrow = 1)
  #   res_c <- matrix(res$cent, nrow = 1)
  # } else {
  #   res_a <- sapply(res, function(r) r$acc, simplify = T) 
  #   res_b <- sapply(res, function(r) r$brier, simplify = T)
  #   res_c <- sapply(res, function(r) r$cent, simplify = T)
  # }
  
  #print(nrow(res))
  #print(sprintf("Length of trees = %d", length(trees)))
  stopifnot(nrow(res) == length(trees))
  correct <- apply(res, 1, median)
  
  imax <- sample(as.list(which(correct == max(correct))), 1)[[1]]
  tree <- trees[[imax]]
  
  cv_acc = median(res[imax,]) / nrow(dfs$train)
  
  data.frame(tree = str_prufer(tree), agree = 1,  
             cv_acc = cv_acc, 
             test_acc = eval_tree(tree, le$pred2, le$GT2) / nrow(dfs$test),
             cv = folds, repeats = repeats)
}

algs <- new.env(parent = emptyenv())
algs[["A"]] = get_A
algs[["B"]] = get_B
algs[["C"]] = get_C
algs[["D"]] = get_D
algs[["E"]] = get_E
algs[["F"]] = get_F

algs1 = new.env(parent = emptyenv())
algs1[["A"]] = get_A
algs1[["B"]] = get_B

algs4 = new.env(parent = emptyenv())
algs4[["A"]] = get_A
algs4[["B"]] = get_B
algs4[["AC"]] = get_AC

algs2 <- new.env(parent = emptyenv())
algs2[["A"]] = get_A
algs2[["B"]] = get_B
algs2[["BR"]] = get_BR
algs2[["D3"]] = get_D3
algs2[["TR"]] = get_TR

algs3 <-  new.env(parent = emptyenv())
algs3[["A"]] = get_A
algs3[["B"]] = get_B
algs3[["C"]] = get_C
algs3[["D"]] = get_D
algs3[["E"]] = get_E
algs3[["F"]] = get_F
algs3[["G"]] = get_G

eval3 <- function(dataset = "satimage",
                  method = "glmnet",
                  algs = algs,
                  run = 0,
                  folds = 10,
                  repeats = 1,
                  calibrate = TRUE,
                  cache_file = NULL) {

  dfs <- timed("read_data", read_wlws(800, dataset, run))
  eval3_dfs (dfs, method = method, algs = algs, run = run, folds = folds, repeats = repeats,
             calibrate = calibrate, cache_file = cache_file)
}

eval3_dfs <- function(dfs,
                      method,
                      algs,
                      run = 0,
                      folds = 10,
                      repeats = 1, calibrate = TRUE, cache_file = NULL, ...) {
	.timing_env$alg <- NA_character_
	le <- timed("pairwise_fit_full", log_exp1_dfs(dfs, method = method, ...))
  eval3_dfs_le(dfs, le, method = method, algs = algs, run = run, folds = folds, repeats = repeats,
               calibrate = calibrate, cache_file = cache_file, ...)
}

eval3_dfs_le <- function(dfs, le, 
                      method,
                      algs, 
                      run = 0,
                      folds = 10,
                      repeats = 1, calibrate = TRUE,
                      # the WLW scales are fitted jointly (K(K-1)/2 parameters),
                      # so only for small K
                      calibrate_wlw = calibrate && le$K <= 7,
                      # if set, the pairwise models and the out-of-fold log-odds
                      # of every algorithm are saved there (saveRDS)
                      cache_file = NULL, ...) {
  print(method)
	Loof_wlw <- NULL   # out-of-fold log-odds of the first algorithm, reused for WLW
	oof <- list()      # per algorithm: out-of-fold log-odds and folds (for cache_file)
	full_df <- map(ls(algs), function(alg) {
	  print(sprintf("Doing %s", alg))
	  .timing_env$alg <- alg
	  trees <- timed("tree_generation", algs[[alg]](le), detail = NA_character_)
    print("tree")
    print(trees[[1]])
	  timed("alg_total", {

	  # cross-validated selection first: it also yields the out-of-fold log-odds
	  # needed to temperature-scale the non-CV baseline (cvf only consumes random
	  # numbers, the baseline below none, so raw results do not depend on the order)
	  dcv <- if (folds > 0) timed("cv_total", cvf(trees, dfs, le, folds, repeats, method, calibrate = calibrate)) else NULL
	  Loof <- if (!is.null(dcv)) attr(dcv, "Loof")
	  if (is.null(Loof_wlw)) Loof_wlw <<- Loof
	  if (!is.null(dcv)) oof[[alg]] <<- list(Loof = Loof, fold = attr(dcv, "fold"))

	  # non-CV baseline for the 3 loss metrics: pick the tree that minimizes
	  # each metric on the full training set (no CV)
	  # The coupling probabilities of a tree on the training set do not depend on
	  # the metric: compute them once per tree and evaluate all metrics on them
	  # (previously log_coupling ran once per metric, 3 times per tree).
	  mets <- ls(metrics)
	  dd <- timed("noncv_baseline", {
	  train_all <- vapply(trees, function(t) {
	    probs1 <- log_coupling(t, le$pred1)
	    vapply(mets, function(met) metrics[[met]](probs1, le$GT1), numeric(1))
	  }, numeric(length(mets)))
	  map(mets, function(met) {
	    train_scores <- train_all[match(met, mets), ]
	    imin <- which.min(train_scores)
	    tree_m <- trees[[imin]]
	    test_probs <- log_coupling(tree_m, le$pred2)

	    data.frame(rep = NA_integer_, tree = str_prufer(tree_m), metric = met,
	               cv = 0, repeats = NA_integer_,
	               train_score = train_scores[imin] / nrow(dfs$train),
	               test_score = metrics[[met]](test_probs, le$GT2) / nrow(dfs$test),
	               test_score_cal = if (calibrate && !is.null(Loof))
	                 timed("calibration", calibrated_test_score(tree_m, met, le$pred2, le$GT2, Loof, dfs$train$class_id),
	                       detail = paste("noncv", met))
	               else NA_real_)
	  }) |> list_rbind()
	  })

	  df <- if (!is.null(dcv)) rbind(dcv, dd) else dd
	  df$alg <- alg
	  df
	  }, detail = "all stages for this algorithm")
	}) |> list_rbind()

	# WLW coupling (full maximal one-vs-one model) — no tree structure, no CV,
	# and (like compute_soft_metrics.R) no train-side evaluation since there's
	# no selection-from-training step to justify one; only test_score is filled.
	if (!exists("wu2_ld")) {
	  library(Rcpp)
	  sourceCpp("decode.cpp")
	}
	.timing_env$alg <- "WLW"
	n2 <- dim(le$pred2)[3]
	wlw_probs2 <- timed("wlw_coupling", sapply(1:n2, function(i) wu2_ld(le$pred2[,,i])) |> t())

	# calibrated WLW: one scale per pair, fitted jointly by minimising the WLW
	# cross-entropy on the out-of-fold log-odds of the training data (the same
	# ones used to calibrate the trees of the first algorithm)
	wlw_cal2 <- if (calibrate_wlw && !is.null(Loof_wlw)) timed("wlw_calibration", {
	  s <- fit_wlw_scales(Loof_wlw, as.integer(dfs$train$class_id))
	  wlw_probs_scaled(aperm(le$pred2, c(3, 1, 2)), s$a, s$pairs)
	})

	dd_wlw <- map(ls(metrics), function(met) {
	  data.frame(rep = NA_integer_, tree = NA_character_, metric = met,
	             cv = 0, repeats = NA_integer_,
	             train_score = NA_real_,
	             test_score = metrics[[met]](wlw_probs2, le$GT2) / nrow(dfs$test),
	             test_score_cal = if (!is.null(wlw_cal2))
	               metrics[[met]](wlw_cal2, le$GT2) / nrow(dfs$test) else NA_real_)
	}) |> list_rbind()
	dd_wlw$alg <- "WLW"

	full_df <- rbind(full_df, dd_wlw)

	full_df$dataset = dfs$dataset
	full_df$run = run

	# le: K, dM (Cohen's d), pred1 / pred2 (training / test log-odds [K, K, n]),
	# GT1 / GT2, L (hyperparameters, e.g. C, gamma and Platt a, b) and, for the
	# SVM, the fitted models; oof[[alg]]: Loof [n_train, K, K] and fold of the
	# first CV repeat of that algorithm
	if (!is.null(cache_file))
	  saveRDS(list(dataset = dfs$dataset, run = run, method = method,
	               le = le[intersect(c("K", "dM", "pred1", "pred2", "GT1", "GT2", "L", "models", "cvdv"), names(le))],
	               oof = oof, results = full_df), cache_file)
	full_df
}

eval3_all <- function(dataset = "satimage",
                     method = "glmnet",
                      file = "alg1",
                      folds = 10,
                      repeats = 1,
                      seed = 123,
                      calibrate = TRUE,
                      suffix = "",
                      workers = getOption("eval3.workers", 10),
                      # if set, every run saves its pairwise models and
                      # out-of-fold log-odds to <cache_dir>/<dataset>_run_XX.rds
                      cache_dir = NULL) {
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  # Self-seeding: previously this function relied on the caller having
  # called set.seed() immediately beforehand, which is fragile (easy to
  # forget, and breaks if the caller does it once before a loop of
  # multiple eval3_all() calls). furrr_options(seed = TRUE) below derives
  # independent per-run L'Ecuyer-CMRG streams from whatever the RNG state
  # is at the time future_map() is called, so fixing that state here makes
  # every eval3_all() call reproducible on its own.
  set.seed(seed)
  algs <- switch(file,
                     "alg1" = algs1,
                     "alg2" = algs2,
                     
                 "alg3" = algs3)
  
  my_globals <- c(
    # Functions defined directly in exp1_log.R
    "log_exp1", 
    "dataset",
    "repeats",
    "folds",
    "log_exp1_dfs", 
    "log_exp1_dfs_glmnet", 
    "log_exp1_dfs_nnet", "log_exp1_dfs_firth", "firth_fit",
    "log_exp1_dfs_svm", "svm_grid_C", "svm_grid_g", "platt_fit", "platt_logloss",
    "svm_dv", "svm_fit01", "inner_folds", "svm_cv_dv",
    "eval_dataset",
    "eval3",
    "eval3_dfs",
    "eval3_dfs_le",
    "read_wlw",
    "prufer_code",
    "eval_tree",
    "log_coupling_internal",
    "cvf",
    "metrics",
    "make5folds",
    # temperature scaling (temperature_scaling.R) and timing helpers
    "calibrate", "calibrated_test_score", "prufer_edges", "tree_diameter",
    "tree_structure", "tree_scores_scaled", "fit_pair_scale", "reference_edge",
    "fit_tree_scales_seq", "tree_predict_scaled", "tree_probs_scaled",
    "wlw_probs_scaled", "fit_wlw_scales",
    ".timing_env", "timed", "timing_table",
    # Functions sourced from external files (CRITICAL)
    "read_wlws",      # From msvm.R? Or somewhere else? Must be accessible.
    "cohens_d",       # From cohen.R
    "do_coupling",    # From log_coupling.R
    "log_coupling",   # From log_coupling.R
    "str_prufer"     # From Prufer.R
    
    # ... Plus all your functions from the main script (prufer_code, eval3, etc.) ...
  )
  
  my_packages = c(
    "igraph",
    "dplyr",
    "purrr",
    "glmnet", # New addition from exp1_log.R
    "e1071",  # RBF-SVM pairwise classifiers (method = "svm")
    "Rcpp"    # needed for wu2_ld (WLW coupling), sourced fresh per worker
    # If "if_else" is used without dplyr::, you may need "dplyr"
  )
  plan("multisession", workers = workers)
  #plan("multisession", workers = 7,
  #                 future.packages = my_packages,
  #               future.globals = "eval3_all")
  #plan(sequential)
  t_start <- Sys.time()
  runs <- future_map(0:19, function(i) {
    .timing_env$log <- list()                      # fresh log for this run
    t0 <- proc.time()[["elapsed"]]
    res <- eval3(dataset = dataset, method = method, algs = algs, run = i, folds = folds,
                 repeats = repeats, calibrate = calibrate,
                 cache_file = if (!is.null(cache_dir)) sprintf("%s/%s_run_%02d.rds", cache_dir, dataset, i))
    tm <- rbind(timing_table(), data.frame(stage = "run_total", alg = NA_character_,
                detail = NA_character_, seconds = proc.time()[["elapsed"]] - t0))
    tm$dataset <- dataset; tm$run <- i
    attr(res, "timing") <- tm                      # collected below
    res
  },
  .options = furrr::furrr_options(seed = TRUE,
                                  packages = my_packages,
                                  globals = my_globals), 
  .progress = T)
  
  df <- list_rbind(runs)
  timing <- list_rbind(lapply(runs, function(r) attr(r, "timing")))
  timing$wall_clock_all_runs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  timing$workers <- workers
  write.csv(timing, file = sprintf("data/timing-%s%s-%s-%d.csv", file, suffix, dataset, repeats),
            row.names = F)

  write.csv(df, 
            file = sprintf("data/%s%s-%s-%d.csv", file, suffix, dataset, repeats),
            quote = F, 
            row.names = F)
}
