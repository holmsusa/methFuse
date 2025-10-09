#' @useDynLib fuseR, .registration = TRUE
NULL


#' Find Optimal Number of Clusters
#'
#' @description
#' Determines the optimal number of clusters to cut a hierarchical clustering tree,
#' based on the selected information criterion (e.g., BIC or AIC).
#'
#' @param tree A sorted clustering tree (matrix or data.frame).
#' @param n Number of samples in the original data.
#' @param method Information criterion method. One of `"BIC"` or `"AIC"`.
#'
#' @return
#' An integer representing the optimal number of clusters.
#'
#' @examples
#' # Example: Determine number of clusters in dummy data set
#' K0 <- matrix(sample(1:200, 80, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(1:200, 80, replace = TRUE), ncol = 5)
#' tree <- fuse.cluster(K0, K1, sort = TRUE)
#' k <- number.of.clusters(tree, ncol(K0), 'BIC')
#' k
#'
#' @export
number.of.clusters <- function(tree, n, method = c("BIC", "AIC")) {

  # --- Input checks ---
  if (!is.matrix(tree) && !is.data.frame(tree)) {
    stop("`tree` must be a matrix or data.frame.")
  }

  if (ncol(tree) < 3) {
    stop("`tree` must have at least 3 columns.")
  }

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a single positive numeric value.")
  }

  # Match and validate method
  method <- match.arg(method)

  # Define constants
  m <- nrow(tree) + 1
  k <- seq(nrow(tree), 1)

  # Set penalty based on method
  penalty <- switch(method,
                    "AIC" = 2 * n,
                    "BIC" = n * log(m * n)
  )

  # Compute information criterion scores
  ics <- penalty * k + 2 * tree[, 3]

  # Determine optimal number of clusters
  optimal_clusters <- as.integer(length(ics) - (which.min(ics) - 1))

  return(optimal_clusters)
}




#' Cut Hierarchical Clustering Tree into Clusters
#'
#' @description Divides the clustering tree into a specified number of clusters.
#' @param tree Sorted clustering tree
#' @param k Number of clusters
#' @return A vector indicating which cluster each element in the original data frame belonged to
#' @examples
#' # Example: Cutting small tree in 2 segments
#' tree <- matrix(c(
#' -1, -2,  49.53106,  49.53106,  1.14473,
#' -3, -4,  78.49604,  78.49604,  1.14473,
#' -5, -6, 147.07154, 147.07154,  1.14473,
#' 1,  2,  72.98287, 201.00997,  1.14473,
#' 4,  3, 106.38879, 454.47029,  1.14473
#' ), ncol = 5, byrow = TRUE)
#'
#' segments <- fuse.cutree(tree, 2)
#' segments
#'
#' @export
fuse.cut.tree <- function (tree, k) {
  # Cuts the given tree into k clusters.
  # Input: matrix tree, numeric integer k
  # Output: integer vector
  return(.Call('cuttree_R', tree, as.integer(k)))
}
