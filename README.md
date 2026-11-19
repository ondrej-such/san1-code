# san1-code

Code accompanying the paper *Speeding inference using minimal one-vs-one classification models* (Ondrej Such, Muhammad Azeem and Martin Bachratý).

## Contents

| Path | Description |
|---|---|
| `reproduce.ipynb` | Recreates the figures and tables of the paper that are computed in R: Figures 3–9 and Tables 1, 2, 3, 5 and 6 (R kernel). |
| `R/` | Functions used by the notebook: the simulated examples of Section 3 (Gaussian classes with LDA, von Mises–Fisher classes on the sphere, Tammes configurations, diameter-3 trees) and the plotting helpers of Figures 7–9. |
| `data/` | Results of the RBF-SVM experiments of Sections 5–7 (per run, dataset, algorithm and metric), and of the simulation of Figure 6. |
| `experiments/` | The code that produced the result files in `data/`. |

Figures 1 and 2 (drawings of trees) and Table 4 (overview of the algorithms) are not computed.

## Running the notebook

The notebook needs R with the [IRkernel](https://irkernel.github.io/) Jupyter kernel and the packages `dplyr`, `tidyr`, `purrr`, `ggplot2`, `patchwork`, `MASS`, `caret`, `future` and `furrr` (its first cell installs the missing ones). Start Jupyter in this folder, so that `R/` and `data/` are found. It runs in under a minute.

* The simulations of Section 3 (Figures 3–5, Tables 1 and 2) are rerun with the seeds of the paper.
* Figure 6 averages 1000 repetitions for each of 50 values of the correlation parameter. By default it reads `data/example2_lda.csv`; set `RECOMPUTE_EXAMPLE2 <- TRUE` in its cell to rerun the simulation (a few minutes, in parallel).
* Figures 7–9 and Tables 3, 5 and 6 are computed from the result files in `data/`.

## Result files

| File | Paper | Produced by |
|---|---|---|
| `alg1svm-{dna,waveform}-1.csv` | Figure 7 | `experiments/run_svm_k3.R` |
| `alg2svm-{satimage,segment}-1.csv` | Figure 8, Table 3 | `eval3_all(file = "alg2")` in `experiments/exp4.R` |
| `alg3svm-<dataset>-7.csv` | Table 5, Figure 9 | `experiments/eval_svm_alg3.R` |
| `dagsvm-svm.csv` | Table 5 (DAGSVM) | `experiments/dagsvm_svm.R` |
| `recal-svm-<dataset>-7.csv` | Table 6 | `experiments/recalib_svm_alg3.R` |
| `example2_lda.csv` | Figure 6 | `lda_example2()` in `R/lda_minimal.R` |

The `alg*` and `recal*` files have one row per run (0–19), algorithm, metric (`error_rate`, `brier_score`, `cross_entropy`) and, where trees are selected by cross-validation, repetition of the cross-validation (`rep`). `test_score` is the score of the selected tree on the test set; `tree` is its Prüfer code. In the `alg1`/`alg2` files the rows with `cv = 10` are selected by 10-fold cross-validation and those with `cv = 0` on the whole training set; WLW (the full one-vs-one model with the coupling of Wu, Lin and Weng) has no selection. In the `recal` files, `method` is `U` for the selection of Section 6 and `R` for the selection among recalibrated trees.

## Rerunning the experiments

The experiments use the benchmark datasets of Wu, Lin and Weng (2004), *Probability estimates for multi-class classification by pairwise coupling*, JMLR 5, in their 20 random splits of 800 training and 1000 test points (LIBSVM format). They are not included here (about 220 MB). The scripts read them from `experiments/unzips/n800/`: the training set of run `i` from `<dataset>.scale-<i>` and the test set from `<dataset>.scale.t-<i>`, for `i = 0, ..., 19` and the datasets `dna`, `waveform`, `satimage`, `segment`, `usps`, `mnist` and `letter`.

They need the R packages `e1071`, `Rcpp` (and a C++ compiler), `igraph`, `glmnet`, `caret`, `dplyr`, `purrr`, `tidyr`, `readr`, `magrittr`, `future`, `furrr` and `ParBayesianOptimization`. Run them in `experiments/`; they write their results to `experiments/data/` (copy them to `data/` for the notebook) and their model caches to `experiments/cache/`.

```bash
cd experiments
mkdir -p data
# Figure 7 (a few minutes)
Rscript run_svm_k3.R
# Figure 8 and Table 3 (slow: all spanning trees are scored in every fold)
Rscript -e "source('exp4.R'); for (ds in c('satimage', 'segment')) eval3_all(dataset = ds, method = 'svm', file = 'alg2', repeats = 1, suffix = 'svm')"
# Table 5 and Figure 9; the argument is the number of parallel workers
# (about 15 minutes with 30 workers); also stores the pairwise SVM log-odds in cache/svm_alg3
Rscript eval_svm_alg3.R 10
# Table 5, DAGSVM (from cache/svm_alg3)
Rscript dagsvm_svm.R
# Table 6 (from cache/svm_alg3; about 1.5 hours with 30 workers)
Rscript recalib_svm_alg3.R 10
```

The pairwise RBF-SVMs are fitted with e1071 (LIBSVM). Their parameters $(C, \gamma)$ and the Platt sigmoid are selected separately for every pair of classes by 5-fold cross-validation; see Section 8 of the paper.
