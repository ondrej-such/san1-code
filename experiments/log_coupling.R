# function to compute posterior probabilies for a given spanning tree
# where we are given pairwise log-odds
#
log_coupling <- function(tree, pred) {
  probs <- log_coupling_internal(tree, pred, 1:(length(tree) + 2))
  stopifnot(sum(is.nan(probs)) == 0)
  # softmax of each row (vectorised; same operations as row by row)
  m <- probs[, 1]
  for (k in seq_len(ncol(probs))[-1]) m <- pmax(m, probs[, k])
  e <- exp(probs - m)
  e / rowSums(e)
}

log_coupling_internal <- function(tree, pred, vertices = 1:(length(tree) + 2)) {
  # tree is spanning tree represented as Prufer code
  # pred is 3-dimensional array of predictions
  # vertices indicate, which numbers occur as labels of the tree
  stopifnot(length(tree) + 2 == length(vertices))
  if (length(tree) >= 1){
    stopifnot(all(tree %in% vertices))
  }
  else {
    stopifnot(length(vertices) == 2)
    # need to be careful, whether rows or columns
    i1 <- vertices[1]
    i2 <- vertices[2]
    M <- cbind(pred[i1,i2,], 0)
    return(M)
  }

  # let's remove vertex from the first edge
  leaf_vertex <- vertices[min(which(!(vertices %in% tree)))]
  new_vertices <- vertices[vertices != leaf_vertex]
  new_tree <- tree[-1]

  prob <- log_coupling_internal(new_tree, pred, new_vertices)
  # prob is a matrix with number of rows equal to number of samples

  K <- length(vertices)
  normalized <- matrix(0, nrow = dim(pred)[3], ncol = K)
  index_leaf_vertex <- which(leaf_vertex == vertices)
  index_joined <- which(tree[1] == vertices)

  # all samples at once (previously a loop over samples): the unnormalised
  # log-score of the leaf is that of the vertex it is joined to plus the
  # log-odds of the leaf against it
  normalized[, -index_leaf_vertex] <- prob
  normalized[, index_leaf_vertex] <- normalized[, index_joined] + pred[leaf_vertex, tree[1], ]
  stopifnot(!any(is.nan(normalized)))
  # we return the matrix of predictions
  return(normalized)
}
