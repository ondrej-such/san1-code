# Shared look of the boxplot figures of Sections 5 and 6 (Figures 7-9): the full
# one-vs-one model (WLW) is the baseline, drawn blue with a blue axis label;
# light green marks the best algorithm by mean and every algorithm that is not
# significantly different from it; the others are grey.

baseline_blue <- "#2166ac"
baseline_fill <- ggplot2::scale_fill_manual(
  values = c(other = "grey85", best = "#a1d99b", baseline = "#6baed6"), guide = "none")

# d: the observations of one panel (columns run, alg, score). The best of the
# algorithms `among` by mean, and every algorithm of `among` whose paired t-test
# against it (over the runs) has p >= alpha, get grp = "best".
best_group <- function(d, among, alpha = 0.05) {
  m <- d |> dplyr::filter(alg %in% among) |> dplyr::group_by(alg) |>
    dplyr::summarise(m = mean(score), .groups = "drop")
  b <- as.character(m$alg[which.min(m$m)])
  sb <- d |> dplyr::filter(alg == b) |> dplyr::arrange(run) |> dplyr::pull(score)
  same <- vapply(as.character(m$alg), function(a) {
    sa <- d |> dplyr::filter(alg == a) |> dplyr::arrange(run) |> dplyr::pull(score)
    dd <- sa - sb
    if (sd(dd) == 0) return(all(dd == 0))          # identical (or constant shift)
    t.test(sa, sb, paired = TRUE)$p.value >= alpha
  }, logical(1))
  green <- as.character(m$alg)[same]
  d |> dplyr::mutate(grp = ifelse(alg == "WLW", "baseline",
                                  ifelse(as.character(alg) %in% green, "best", "other")))
}

# x-axis labels: WLW bold and blue
baseline_axis <- function(lv)
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    colour = ifelse(lv == "WLW", baseline_blue, "grey30"),
    face = ifelse(lv == "WLW", "bold", "plain")))

# the three metrics in the order and with the names of the paper
metric_factor <- function(m)
  factor(m, levels = c("error_rate", "brier_score", "cross_entropy"),
         labels = c("test error", "Brier score", "cross-entropy"))
