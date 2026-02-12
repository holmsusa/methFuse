#' Find Optimal Number of Clusters
#'
#' @description
#' Determines the optimal number of clusters to cut a hierarchical clustering tree,
#' based on the selected information criterion (e.g., BIC or AIC).
#'
#' @param tree Clustering tree of class \code{hclust}.
#' @param n Number of samples in the original data.
#' @param method Information criterion method. One of `"BIC"` or `"AIC"`.
#'
#' @return
#' An integer representing the optimal number of clusters.
#'
#' @examples
#' # Example: Determine number of clusters in dummy data set
#' set.seed(1234)
#' K0 <- matrix(
#'   rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
#'   nrow = 100, byrow = TRUE
#' )
#' K1 <- matrix(
#'   rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
#'   nrow = 100, byrow = TRUE
#' )
#' tree <- fuse.cluster(K0, K1)
#' k <- number.of.clusters(tree, ncol(K0), 'BIC')
#' k
#'
#' @export
number.of.clusters <- function(tree, n, method = c("BIC", "AIC")) {

  # --- Input checks ---
  if(class(tree)[1] != "hclust") {
    stop("Input must be hclust.")
  }

  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a single positive numeric value.")
  }

  # Match and validate method
  method <- match.arg(method)

  # Define constants
  m <- nrow(tree$merge) + 1
  k <- seq(nrow(tree$merge), 1)

  # Set penalty based on method
  penalty <- switch(method,
                    "AIC" = 2 * n,
                    "BIC" = n * log(m * n)
  )

  # Compute information criterion scores
  ics <- penalty * k + 2 * tree$height

  # Determine optimal number of clusters
  optimal_clusters <- as.integer(length(ics) - (which.min(ics) - 1))

  return(optimal_clusters)
}




#' Cut Hierarchical Clustering Tree into Clusters
#'
#' @description Divides the clustering tree into a specified number of clusters.
#' @param tree Clustering tree of class \code{hclust}.
#' @param k Number of clusters
#' @return A vector indicating which cluster each element in the original data frame belonged to
#' @examples
#' # Example: Cutting small tree in 2 segments
#' set.seed(1234)
#' K0 <- matrix(
#'  rep(c(sample(0:20, 200, replace = TRUE), sample(20:40, 200, replace = TRUE)), 2),
#'  nrow = 100, byrow = TRUE
#' )
#' K1 <- matrix(
#'  rep(c(sample(20:40, 200, replace = TRUE), sample(0:20, 200, replace = TRUE)), 2),
#'  nrow = 100, byrow = TRUE
#' )
#' tree <- fuse.cluster(K0, K1)
#'
#' segments <- fuse.cut.tree(tree, 2)
#' segments
#'
#' @export
fuse.cut.tree <- function (tree, k) {
  # Cuts the given tree into k clusters.
  # Input: matrix tree, numeric integer k
  # Output: integer vector
  parsed_tree <- as.matrix(cbind(tree$merge, tree$height))

  return(.Call('cuttree_R', parsed_tree, as.integer(k)))
}
