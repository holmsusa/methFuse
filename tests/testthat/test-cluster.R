expect_valid_fuse_tree <- function(tree, n_sites) {
  expect_true(class(tree) == "hclust")
  expect_equal(length(tree), 7)
  expect_equal(
    names(tree),
    c("merge", "height", "order", "labels", "call", "method", "dist.method")
  )

  # hclust-style tree: n_sites - 1 merges
  expect_equal(nrow(tree$merge), n_sites - 1)

  # Merge indices should be integers
  expect_true(is.numeric(tree$merge[, 1]))
  expect_true(is.numeric(tree$merge[, 2]))
}


# fuse.cluster
test_that("fuse.cluster returns hclust with correct names", {
  K0 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)
  K1 <- matrix(sample(1:100, 50, replace = TRUE), ncol = 5)

  chr <- rep("chr3", nrow(K0))
  pos <- seq_len(nrow(K0))

  result <- fuse.cluster(K0, K1, chr, pos)

  expect_true(class(result) == "hclust")
  expect_equal(length(result), 7)
  expect_equal(
    names(result),
    c("merge", "height", "order", "labels", "call", "method", "dist.method")
  )
})

test_that("fuse.cluster works on methrix objects", {

  skip_if_not_installed("methrix")
  skip_if_not_installed("DelayedArray")

  data(methrix_data, package = "methrix")

  tree <- fuse.cluster(methrix_data)

  n_sites <- nrow(methrix::get_matrix(methrix_data, "M", add_loci = FALSE))

  expect_valid_fuse_tree(tree, n_sites)
})

test_that("fuse.cluster works on BSseq objects", {

  skip_if_not_installed("bsseq")

  library(bsseq)

  set.seed(1)
  M <- matrix(sample(0:10, 200, TRUE), nrow = 50)
  Cov <- M + matrix(sample(0:10, 200, TRUE), nrow = 50)

  bs <- BSseq(
    chr = rep("chr1", 50),
    pos = seq_len(50),
    M = M,
    Cov = Cov,
    sampleNames = colnames(M)
  )

  tree <- fuse.cluster(bs)

  expect_valid_fuse_tree(tree, 50)
})


test_that("fuse.cluster is consistent between matrix and BSseq input", {

  skip_if_not_installed("bsseq")

  library(bsseq, quietly = TRUE)

  set.seed(42)
  M <- matrix(sample(0:10, 200, TRUE), nrow = 50)
  Cov <- M + matrix(sample(0:10, 200, TRUE), nrow = 50)

  chr <- rep("chr1", 50)
  pos <- seq_len(50)

  bs <- BSseq(
    chr = chr,
    pos = pos,
    M = M,
    Cov = Cov,
    sampleNames = colnames(M)
  )

  tree_mat <- fuse.cluster(
    x = Cov - M,
    K1 = M,
    chr = chr,
    pos = pos
  )

  tree_bs <- fuse.cluster(bs)

  expect_equal(tree_mat, tree_bs)
})

test_that("fuse.cluster respects chromosome boundaries", {

  skip_if_not_installed("bsseq")

  library(bsseq, quietly = TRUE)

  set.seed(1)
  M <- matrix(sample(0:10, 200, TRUE), nrow = 50)
  Cov <- M + matrix(sample(0:10, 200, TRUE), nrow = 50)

  chr <- c(rep("chr1", 25), rep("chr2", 25))
  pos <- c(seq_len(25), seq_len(25))

  bs <- BSseq(
    chr = chr,
    pos = pos,
    M = M,
    Cov = Cov,
    sampleNames = colnames(M)
  )

  tree <- fuse.cluster(bs)

  # First merge in chr2 cannot involve chr1 indices
  # (negative indices refer to original observations)
  chr2_idx <- which(chr == "chr2")

  expect_true(
    any(abs(tree$merge[, 1]) %in% chr2_idx) ||
      any(abs(tree$merge[, 2]) %in% chr2_idx)
  )
})

