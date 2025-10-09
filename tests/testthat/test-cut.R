# number.of.clusters
test_that("number.of.clusters works correctly on valid input", {
  # Create a dummy sorted clustering tree
  # Columns: assume [merge1, merge2, criterion_value, ...]
  tree <- matrix(c(
    1, 2, 5,
    3, 4, 3,
    5, 6, 10,
    7, 8, 7
  ), ncol = 3, byrow = TRUE)

  n <- 10

  # Should return a single integer
  result_bic <- number.of.clusters(tree, n, method = "BIC")
  result_aic <- number.of.clusters(tree, n, method = "AIC")

  expect_type(result_bic, "integer")
  expect_type(result_aic, "integer")
  expect_true(result_bic > 0 && result_bic <= nrow(tree) + 1)
  expect_true(result_aic > 0 && result_aic <= nrow(tree) + 1)

  # The two methods can give different results, but both must be valid integers
  expect_true(is.finite(result_bic))
  expect_true(is.finite(result_aic))
})

test_that("number.of.clusters validates inputs", {
  tree <- matrix(1:9, ncol = 3)

  # Invalid tree type
  expect_error(number.of.clusters(list(1, 2, 3), 5, "BIC"),
               "`tree` must be a matrix or data.frame")

  # Not enough columns
  expect_error(number.of.clusters(matrix(1:6, ncol = 2), 5, "BIC"),
               "`tree` must have at least 3 columns")

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
  tree <- matrix(c(1, 2, 0), ncol = 3)
  result <- number.of.clusters(tree, 5, "BIC")
  expect_true(result %in% 1:2)

  # Data frame input works the same
  df_tree <- as.data.frame(tree)
  expect_equal(number.of.clusters(tree, 5, "AIC"),
               number.of.clusters(df_tree, 5, "AIC"))
})


# fuse.cut.tree
test_that("fuse.cut.tree returns integer vector of correct length", {
  tree <- matrix(c(
    -1, -2,  49.5,  49.5,  1.1,
    -3, -4,  78.5,  78.5,  1.1,
    -5, -6, 147.0, 147.0,  1.1,
    1,  2,  72.9, 201.0,  1.1,
    4,  3, 106.3, 454.4,  1.1
  ), ncol = 5, byrow = TRUE)

  segments <- fuse.cut.tree(tree, 3)

  expect_true(is.integer(segments))
  expect_equal(length(segments), 6)
})

#
