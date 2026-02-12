# number.of.clusters
test_that("number.of.clusters works correctly on valid input", {
  # Create a dummy sorted clustering tree
  # Columns: assume [merge1, merge2, criterion_value, ...]
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)

  chr <- rep("chr3", nrow(K0))
  pos <- seq_len(nrow(K0))

  tree <- fuse.cluster(K0, K1, chr, pos)

  n <- 10

  # Should return a single integer
  result_bic <- number.of.clusters(tree, n, method = "BIC")
  result_aic <- number.of.clusters(tree, n, method = "AIC")

  expect_type(result_bic, "integer")
  expect_type(result_aic, "integer")
  expect_true(result_bic > 0 && result_bic <= nrow(tree$merge) + 1)
  expect_true(result_aic > 0 && result_aic <= nrow(tree$merge) + 1)

  # The two methods can give different results, but both must be valid integers
  expect_true(is.finite(result_bic))
  expect_true(is.finite(result_aic))
})

test_that("number.of.clusters validates inputs", {
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)

  chr <- rep("chr3", nrow(K0))
  pos <- seq_len(nrow(K0))

  tree <- fuse.cluster(K0, K1, chr, pos)


  # Invalid tree type
  expect_error(number.of.clusters(list(1, 2, 3), 5, "BIC"),
               "Input must be hclust.")

  # Not enough columns
  expect_error(number.of.clusters(matrix(1:6, ncol = 2), 5, "BIC"),
               "Input must be hclust.")

  # Invalid n
  expect_error(number.of.clusters(tree, -1, "BIC"),
               "`n` must be a single positive numeric value")
  expect_error(number.of.clusters(tree, c(1, 2), "BIC"),
               "`n` must be a single positive numeric value")
  expect_error(number.of.clusters(tree, "a", "BIC"),
               "`n` must be a single positive numeric value")

  # Invalid method
  expect_error(number.of.clusters(tree, 5, "XYZ"),
               "should be one of")
})

test_that("number.of.clusters edge behavior is consistent", {
  # tree with one row
  tree <- list(
    merge = matrix(c(-1, -1), nrow = 1),
    height = 5,
    order = 1,
    labels = "chr3.43",
    call = "call",
    method = "method",
    dist.method = "dist.method"
  )
  attr(tree, "class") <- "hclust"

  result <- number.of.clusters(tree, 5, "BIC")
  expect_true(result %in% 1:2)
})


# fuse.cut.tree
test_that("fuse.cut.tree returns integer vector of correct length", {
  tree <- list(
    merge = matrix(c(-1, -2,
                     -3, -4,
                     -5, -6,
                     1, 2,
                     3, 4), ncol = 2, byrow = T),
    height = c(2, 4, 5, 10, 13),
    order = seq_len(6),
    labels = paste("chr3", seq_len(6), sep = "."),
    call = "call",
    method = "method",
    dist.method = "dist.method"
  )
  attr(tree, "class") <- "hclust"

  segments <- fuse.cut.tree(tree, 3)

  expect_true(is.integer(segments))
  expect_equal(length(segments), 6)
})


