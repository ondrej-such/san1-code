# RBF-SVM pairwise classifiers on the three-class datasets (Section 5.1):
# Algorithms A, B and the full WLW model, 20 runs, 10-fold CV, 1 repeat.
# Output: data/alg1svm-{dna,waveform}-1.csv
suppressPackageStartupMessages(source("exp4.R"))
for (ds in c("dna", "waveform")) {
  t0 <- Sys.time()
  eval3_all(dataset = ds, method = "svm", file = "alg1", repeats = 1, suffix = "svm")
  message(ds, " done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
message("ALL DONE")
