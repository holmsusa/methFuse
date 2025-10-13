# fuse.cluster
test_that("fuse.cluster returns matrix with correct column names", {
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)

  result <- fuse.cluster(K0, K1)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), 5)
  expect_equal(
    colnames(result),
    c("m1", "m2", "logl_tot", "logl_merge", "genomic_dist")
  )
})

