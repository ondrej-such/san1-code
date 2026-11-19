
source("msvm.R")
library(glmnet)
source("cohen.R")
library(caret)

log_exp1 <- function(N = 300, dataset = "mnist", i = 0, method, normalize = NA, ...) {
  dfs <- read_wlws(N, dataset, i)
  log_exp1_dfs(dfs, method = method, normalize = normalize, ...)
}


log_exp1_dfs_nnet <- function(dfs, 
                            L = NULL, 
                            normalize = NA, 
                            tune_grid = expand.grid(size = c(3, 15), 
                                    decay = c(0.01, 0.1)), ...) {
# Arguments:
#   L 
#       environment specifying decay and size for each pair of classes, maybe NULL
#   normalize 
#       NULL if no normalization is provided, but we want one
#       NA if we don't want to normalize data (i.e. already normalized, e.g. WLW data)
#       list with center_scaled and scale_scaled attributes specifying normalization
  if (is.null(L)) {
    L <- new.env(parent = emptyenv())
  }

  dataset <- dfs$dataset
  N <- dfs$N
  i <- dfs$run
  
  # K <- max(dfs$train$class_id)
  K <- dfs$K
  
  pred1 <- array(0, dim = c(K, K, nrow(dfs$train)))
  pred2 <- array(0, dim = c(K, K, nrow(dfs$test)))
  dM <- matrix(nrow = K * (K-1) / 2, ncol = 4)
  row <- 0
  
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      
      # first create the model
      ix <- dfs$train$class_id %in% c(i,j)
      x <- as.matrix(dfs$train[ix, setdiff(names(dfs$train), "class_id")])
      y <- matrix(dfs$train$class_id[ix] == i, ncol = 1 ) 
      # print(dim(x))
      # print(dim(y))
      # print(length(y))
      idx <- sprintf("%d:%d", i,j)

      if (is.null(normalize)) {
         x_norm <- scale(x)
         ij_norm <- list(scaled_center = attr(x_norm, "scaled:center"), 
                           scaled_scale = attr(x_norm, "scaled:scale"))
      } else if (is.na(normalize)) {
        # do nothing
         x_norm <- x
      } else {
         x_norm <- scale(x, normalize$scaled_center, normalize$scaled_scale)
         ij_norm <- normalize
      }
      df <- as.data.frame(x_norm)
      df$class_id = as.factor(dfs$train$class_id[ix]) |> droplevels()
      if (is.null(L[[idx]])) {
         # Set up cross-validation (5-fold)
        ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
        # tune_grid <- expand.grid(size = c(3, 15), decay = c(0.01, 0.1))
        seed <- sample(1:100000, 1)
        set.seed(seed)
        cv_model <- train(
              class_id ~ ., 
              data = df,
              method = "nnet",
              tuneGrid = tune_grid,
              trControl = ctrl,
              trace = FALSE,
              maxit = 200
             # Standardize features (crucial for MLPs)
              #)  #BUGBUG
        )
        # cv.result <- cv.glmnet(x, y, family = "binomial", alpha = 0)
        L[[idx]] <- list(decay = cv_model$bestTune$decay, size = cv_model$bestTune$size, seed = seed)
        model <- cv_model$finalModel
      }  else {
        hparam <- L[[idx]] 
        set.seed(hparam$seed)
        model <- nnet(class_id  ~. , data = df, size = hparam$size, decay = hparam$decay, trace = FALSE)
        # print(sprintf("using supplied lambda %f", L[i,j]))
      }
      
      # let's compute cohen's d's on TRAIN data
      # x1 <- as.matrix(dfs$train[, setdiff(names(dfs$train), "class_id")])
      p1 <- predict(model, dfs$train, type = "raw")
      scores <- qlogis(pmin(pmax(p1, 1e-15), 1 - 1e-15))
      pred1[i,j,] <- -scores
      pred1[j,i,] <- scores
      GT <- dfs$train$class_id
      Scores_i <- scores[GT[ix] == i]
      Scores_j <- scores[GT[ix] == j]
      
      stopifnot (length(Scores_i) > 1 && length(Scores_j) > 1)
      Cd_value <- cohens_d(Scores_i, Scores_j)
      row <- row + 1
      dM[row,] = c(i,j, Cd_value, if(exists("cv_model")) mean(cv_model$resample$Accuracy) else NA)
       
      x2 <- as.matrix(dfs$test[ , setdiff(names(dfs$test), "class_id")])
      if (exists("ij_norm")) {
         if (!is.null(ij_norm) && (length(ij_norm) == 2)) {
            x2 <- scale(x2, center = ij_norm$scaled_center, scale = ij_norm$scaled_scale)
            ij_norm <- NULL
         }
      } 
      df <- as.data.frame(x2)
      df$class_id = as.factor(dfs$test$class_id)
      p2 <- predict(model, x2, type = "raw") 
      p <- qlogis(pmin(pmax(p2, 1e-15), 1 - 1e-15))
      # stopifnot(all(p >= 0))
      pred2[i,j,] <- -p # in i-th row is the probability of class i
      pred2[j,i,] <- p # in j-th row is the probability of class j

      ix2 <- dfs$test$class_id %in% c(i,j)
      pc2 <- predict(model, df[ix2,], type = "class")
      # print(head(pc2))
      print(sprintf("%d:%d acc = %f size = %d decay = %f", i,j, mean(dfs$test$class_id[ix2] == pc2), L[[idx]]$size,
      L[[idx]]$decay))
    }
  }
  stopifnot(row == nrow(dM) )
  colnames(dM) <- c("i", "j", "cohend", "Accuracy")
  return(list(K = K, 
              dM = dM, 
              pred1 = pred1, 
              pred2 = pred2,
              GT1 = dfs$train$class_id,
              GT2 = dfs$test$class_id, 
              L = L))
}


log_exp1_dfs <- function(dfs, L = NULL, method, normalize = NA, ...) {
  stopifnot(method %in% c("glmnet", "nnet", "firth", "svm"))
  if (method == "svm") {
     res <- log_exp1_dfs_svm(dfs, L = L, normalize = normalize, ...)
  } else if (method == "firth") {
     res <- log_exp1_dfs_firth(dfs, L = L, normalize = normalize, ...)
  } else if (method == "nnet") {
     res <- log_exp1_dfs_nnet(dfs, L = L , normalize = normalize , ...)
  } else if (method == "glmnet")  {
     res <-log_exp1_dfs_glmnet(dfs, L = L , normalize = normalize , ...)
  } else {
    stopifnot(FALSE)
  }
  res$method = method
  res
}

log_exp1_dfs_glmnet <- function(dfs, L = NULL , normalize = NA, ...) {
  if (is.null(L)) {
    L <- matrix(NA, nrow = dfs$K, ncol = dfs$K)
  }
  dataset <- dfs$dataset
  N <- dfs$N
  i <- dfs$run
  
  K <- max(dfs$train$class_id)
  
  #tree1 <- sample(K, K - 2, replace = TRUE)
  
  pred1 <- array(0, dim = c(K, K, nrow(dfs$train)))
  pred2 <- array(0, dim = c(K, K, nrow(dfs$test)))
  dM <- matrix(nrow = K * (K-1) / 2, ncol = 3)
  row <- 0
  
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      
      # first create the model
      ix <- dfs$train$class_id %in% c(i,j)
      x <- as.matrix(dfs$train[ix, setdiff(names(dfs$train), "class_id")])
      y <- matrix(dfs$train$class_id[ix] == i, ncol = 1 ) 
      # print(dim(x))
      # print(dim(y))
      # print(length(y))
      if (is.na(L[i,j])) {
        cv.result <- cv.glmnet(x, y, family = "binomial", alpha = 0)
        L[i,j] <- cv.result$lambda.min
      }  else {
        # print(sprintf("using supplied lambda %f", L[i,j]))
      }
      lambda <- L[i,j]
      # print(lambda)
      stopifnot(!is.na(lambda) & !is.nan(lambda))
      model <- glmnet(x, y, lambda = lambda, family = "binomial", alpha = 0)
      coeffs <- coef(model)
      if (sum(coeffs != 0) <= 1) {
        print(model)
        print(dim(x))
        print(dim(y))
        print(table(y))
        print(head(x))
        print(svd(x)$d)
        print(svd(x[y, ])$d)
        print(svd(x[!y, ])$d)
        stop("cannot continue with empty model")
      }
      stopifnot(!any(is.na(coeffs)))
      stopifnot(!any(is.nan(coeffs)))
      
      # let's compute cohen's d's on TRAIN data
      x1 <- as.matrix(dfs$train[, setdiff(names(dfs$train), "class_id")])
      scores <- predict(model, x1, type = "link")
      pred1[i,j,] <- scores
      pred1[j,i,] <- -scores
      GT <- dfs$train$class_id
      Scores_i <- scores[GT == i]
      Scores_j <- scores[GT == j]
      
      stopifnot (length(Scores_i) > 1 && length(Scores_j) > 1)
      Cd_value <- cohens_d(Scores_i, Scores_j)
      row <- row + 1
      dM[row,] = c(i,j, Cd_value)
       
      x2 <- as.matrix(dfs$test[ , setdiff(names(dfs$test), "class_id")])
      p <- predict(model, x2, type = "link") 
      # stopifnot(all(p >= 0))
      pred2[i,j,] <- p # in i-th row is the probability of class i
      pred2[j,i,] <- -p # in j-th row is the probability of class j
    }
  }
  stopifnot(row == nrow(dM) )
  colnames(dM) <- c("i", "j", "cohend")
  return(list(K = K, 
              dM = dM, 
              pred1 = pred1, 
              pred2 = pred2,
              GT1 = dfs$train$class_id,
              GT2 = dfs$test$class_id, 
              L = L))
}



# we are ready for coupling

source("Prufer.R")

source("log_coupling.R")

eval_dataset <- function(N = 800, dataset = "dna", i = 0) {
  pairwise <- exp1(N, dataset, i)
  K <- pairwise$K
  
  for (i in 1:min(10, K^(K-2))) {
    tree1 <- sample(K, K-2, replace = TRUE)
    p <- do_coupling(tree1, pairwise$pred)
    pv <- apply(p, 1, which.max)
    accuracy = mean(pv == pairwise$ground_truth)
    print(sprintf("Accuracy is %f", accuracy))
  }
}



# ---- Firth (Jeffreys-penalised) logistic regression as pairwise classifier ----
# Has no tuning parameter, so L (the lambda matrix) is returned unchanged.
firth_fit <- function(x, y, maxit = 200, tol = 1e-9) {
  keep <- apply(x, 2, function(v) length(unique(v)) > 1)     # drop columns constant within the pair
  if (!all(keep)) {
    b <- rep(0, ncol(x) + 1); b[c(TRUE, keep)] <- firth_fit(x[, keep, drop = FALSE], y, maxit, tol); return(b)
  }
  X <- cbind(1, x); b <- rep(0, ncol(X))
  for (it in seq_len(maxit)) {
    p <- plogis(drop(X %*% b)); W <- p * (1 - p)
    XtWX <- crossprod(X, X * W)
    XtWXi <- tryCatch(solve(XtWX), error = function(e) solve(XtWX + 1e-8 * mean(diag(XtWX)) * diag(ncol(X))))
    h <- W * rowSums((X %*% XtWXi) * X)
    step <- drop(XtWXi %*% crossprod(X, y - p + h * (0.5 - p)))
    m <- max(abs(step)); if (m > 5) step <- step * 5 / m
    b <- b + step
    if (max(abs(step)) < tol) break
  }
  b
}

# --- RBF-SVM pairwise classifiers with Platt log-odds -------------------------
# As in Wu, Lin & Weng (2004, Sec. 7.1): decision values are turned into
# log-odds by a Platt sigmoid fitted on cross-validated decision values. For
# each pair, (C, gamma) is chosen on a grid by 5-fold CV within the pair: the
# pooled CV decision values of every grid point get a Platt sigmoid and the
# grid point with the smallest CV log-loss wins (its sigmoid is kept). When L
# supplies (C, gamma) and the sigmoid (the CV folds of eval3, like glmnet's
# lambda), they are reused and only the SVM is refitted; pred1 and dM, which
# the CV folds do not use, are then left as NA. Cohen's d (dM) is computed on
# the pooled CV log-odds of the selected grid point, not on the training fit.
# 8 x 10 grid. The first version used C = 2^(-1, 1, ..., 7) and gamma =
# 2^(-2:2) / m; its edges were selected for up to a third of the pairs
# (results kept as *-grid25.csv). Wu, Lin & Weng (2004) search C, gamma in
# 2^(-5, -3, ..., 15).
svm_grid_C <- 2^seq(-5, 9, by = 2)
svm_grid_g <- 2^(-5:4)            # multiples of 1 / median squared distance

# Platt (2000) with smoothed targets; returns (a, b) with log-odds = a * f + b
platt_fit <- function(f, y) {
  np <- sum(y); nn <- length(y) - np
  t <- ifelse(y == 1, (np + 1) / (np + 2), 1 / (nn + 2))
  nll <- function(ab) { z <- ab[1] * f + ab[2]; sum(pmax(z, 0) + log1p(exp(-abs(z))) - t * z) }
  gr <- function(ab) { r <- plogis(ab[1] * f + ab[2]) - t; c(sum(r * f), sum(r)) }
  optim(c(1, 0), nll, gr, method = "BFGS", control = list(reltol = 1e-12, maxit = 500))$par
}

platt_logloss <- function(ab, f, y) {
  z <- ab[1] * f + ab[2]
  mean(pmax(z, 0) + log1p(exp(-abs(z))) - y * z)
}

# decision values oriented so that positive means y = 1
svm_dv <- function(m, x) {
  d <- attr(predict(m, x, decision.values = TRUE), "decision.values")
  if (colnames(d)[1] == "1/0") drop(d) else -drop(d)
}

svm_fit01 <- function(x, y, C, g)
  e1071::svm(x, factor(y, levels = c(0, 1)), type = "C-classification", kernel = "radial",
             cost = C, gamma = g, scale = FALSE)

# stratified inner folds
inner_folds <- function(y, k = 5) {
  f <- integer(length(y))
  for (v in unique(y)) { ix <- which(y == v); f[ix] <- sample(rep_len(1:k, length(ix))) }
  f
}

svm_cv_dv <- function(x, y, fold, C, g) {
  f <- numeric(length(y))
  for (k in unique(fold)) f[fold == k] <- svm_dv(svm_fit01(x[fold != k, , drop = FALSE], y[fold != k], C, g),
                                                 x[fold == k, , drop = FALSE])
  f
}

log_exp1_dfs_svm <- function(dfs, L = NULL, normalize = NA, ...) {
  K <- max(dfs$train$class_id)
  full <- is.null(L)
  if (full) L <- list(C = matrix(NA, K, K), gamma = matrix(NA, K, K), a = matrix(NA, K, K), b = matrix(NA, K, K))
  fn <- setdiff(names(dfs$train), "class_id")
  X1 <- as.matrix(dfs$train[, fn]); X2 <- as.matrix(dfs$test[, fn])
  pred1 <- array(0, dim = c(K, K, nrow(X1))); pred2 <- array(0, dim = c(K, K, nrow(X2)))
  dM <- matrix(nrow = K * (K - 1) / 2, ncol = 3); row <- 0
  GT <- dfs$train$class_id
  # full fit only: the fitted SVMs and the pooled 5-fold CV decision values of
  # the selected (C, gamma), on which the Platt sigmoid was fitted (rows of the
  # training set in `rows`); named "i-j"
  models <- list(); cvdv <- list()
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    ix <- GT %in% c(i, j)
    x <- X1[ix, , drop = FALSE]; y <- as.integer(GT[ix] == i)
    if (full) {
      fold <- inner_folds(y)
      g0 <- 1 / median(as.vector(dist(x))^2)
      best <- list(loss = Inf)
      for (C in svm_grid_C) for (g in g0 * svm_grid_g) {
        f <- svm_cv_dv(x, y, fold, C, g)
        ab <- platt_fit(f, y); loss <- platt_logloss(ab, f, y)
        if (loss < best$loss) best <- list(loss = loss, C = C, g = g, ab = ab, f = f)
      }
      L$C[i, j] <- best$C; L$gamma[i, j] <- best$g; L$a[i, j] <- best$ab[1]; L$b[i, j] <- best$ab[2]
    }
    m <- svm_fit01(x, y, L$C[i, j], L$gamma[i, j])
    s2 <- L$a[i, j] * svm_dv(m, X2) + L$b[i, j]
    pred2[i, j, ] <- s2; pred2[j, i, ] <- -s2
    row <- row + 1
    if (full) {
      s1 <- L$a[i, j] * svm_dv(m, X1) + L$b[i, j]
      pred1[i, j, ] <- s1; pred1[j, i, ] <- -s1
      # Cohen's d from the cross-validated log-odds of the selected (C, gamma):
      # the in-sample decision values of an RBF-SVM are overfitted and their d
      # carries no information about how hard the pair is (check_A_svm.R)
      s_cv <- best$ab[1] * best$f + best$ab[2]
      dM[row, ] <- c(i, j, cohens_d(s_cv[y == 1], s_cv[y == 0]))
      models[[paste(i, j, sep = "-")]] <- m
      cvdv[[paste(i, j, sep = "-")]] <- list(rows = which(ix), f = best$f, y = y)
    } else {
      pred1[i, j, ] <- pred1[j, i, ] <- NA; dM[row, ] <- c(i, j, NA)
    }
  }
  colnames(dM) <- c("i", "j", "cohend")
  out <- list(K = K, dM = dM, pred1 = pred1, pred2 = pred2, GT1 = GT, GT2 = dfs$test$class_id, L = L)
  if (full) { out$models <- models; out$cvdv <- cvdv }
  out
}

log_exp1_dfs_firth <- function(dfs, L = NULL, normalize = NA, ...) {
  K <- max(dfs$train$class_id)
  if (is.null(L)) L <- matrix(NA, K, K)
  fn <- setdiff(names(dfs$train), "class_id")
  X1 <- as.matrix(dfs$train[, fn]); X2 <- as.matrix(dfs$test[, fn])
  pred1 <- array(0, dim = c(K, K, nrow(X1))); pred2 <- array(0, dim = c(K, K, nrow(X2)))
  dM <- matrix(nrow = K * (K - 1) / 2, ncol = 3); row <- 0
  GT <- dfs$train$class_id
  for (i in 1:(K - 1)) for (j in (i + 1):K) {
    ix <- GT %in% c(i, j)
    b <- firth_fit(X1[ix, , drop = FALSE], as.integer(GT[ix] == i))
    s1 <- drop(cbind(1, X1) %*% b); s2 <- drop(cbind(1, X2) %*% b)
    pred1[i, j, ] <- s1; pred1[j, i, ] <- -s1
    pred2[i, j, ] <- s2; pred2[j, i, ] <- -s2
    row <- row + 1
    dM[row, ] <- c(i, j, cohens_d(s1[GT == i], s1[GT == j]))
  }
  colnames(dM) <- c("i", "j", "cohend")
  list(K = K, dM = dM, pred1 = pred1, pred2 = pred2, GT1 = GT, GT2 = dfs$test$class_id, L = L)
}
