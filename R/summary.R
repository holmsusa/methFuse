#' @useDynLib fuseR, .registration = TRUE
NULL

#' Summarize FUSE Segmentation Results
#'
#' @description
#' Summarizes FUSE segmentation results into one row per segment,
#' including genomic coordinates, CpG count, segment length, average
#' methylation (beta), and stability flag based on likelihood testing.
#' Also returns per-sample methylation estimates for each segment.
#'
#' @param K0 Integer matrix of unmethylated counts.
#' @param K1 Integer matrix of methylated counts.
#' @param chr Character vector giving the chromosome for each site.
#' @param pos Numeric vector giving genomic coordinates for each site.
#' @param segments Integer vector (length = nrow(K0)) giving segment IDs for each site.
#'
#' @return
#' A list with two elements:
#' \describe{
#'   \item{summary}{A data frame with one row per segment and columns:}
#'     \describe{
#'       \item{Segment}{Segment ID}
#'       \item{Chr}{Chromosome}
#'       \item{Start}{Start genomic coordinate}
#'       \item{End}{End genomic coordinate}
#'       \item{CpGs}{Number of CpGs in the segment}
#'       \item{Length}{Genomic length (End - Start + 1)}
#'       \item{Beta}{Average methylation across samples and CpGs}
#'       \item{Stable}{Logical indicator (TRUE if segment is stable, else FALSE)}
#'     }
#'   \item{betas_per_segment}{Matrix of per-sample methylation estimates for each segment
#'   (rows = segments, columns = samples).}
#' }
#'
#' @examples
#' K0 <- matrix(sample(1:200, 125, replace = TRUE), ncol = 5)
#' K1 <- matrix(sample(1:200, 125, replace = TRUE), ncol = 5)
#' tree <- fuse.cluster(K0, K1)
#' segments <- fuse.cut.tree(tree, 40)
#' res <- fuse.summary(K0, K1, rep("chr1", nrow(K0)), 1:nrow(K0), segments)
#' head(res$summary)
#' head(res$betas_per_segment)
#'
#' @export
fuse.summary <- function(K0, K1, chr, pos, segments) {
  # --- Input validation ---
  stopifnot(
    is.matrix(K0),
    is.matrix(K1),
    all(dim(K0) == dim(K1)),
    is.character(chr),
    length(chr) == nrow(K0),
    is.numeric(pos),
    length(pos) == nrow(K0),
    (is.numeric(segments) || is.integer(segments)),
    length(segments) == nrow(K0)
  )


  # --- Initialization ---
  N <- ncol(K0)
  n_segments <- length(unique(segments))

  # --- Helper function ---
  safe_sum <- function(x) sum(x, na.rm = TRUE)

  # --- Model 1: Independent CpG model ---
  phat_ind <- .Call("bino_div_R", K1, K1 + K0)
  logL_ind <- -(.Call("bino_xlogy_R", K1, phat_ind) + .Call("bino_xlogy_R", K0, 1 - phat_ind))
  LL1 <- tapply(rowSums(logL_ind), segments, safe_sum)
  df1 <- tapply(rep(N, nrow(K0)), segments, safe_sum)

  # --- Model 2: Common segment model ---
  K0_seg <- rowsum(K0, group = segments, na.rm = TRUE)
  K1_seg <- rowsum(K1, group = segments, na.rm = TRUE)
  group_sizes <- as.numeric(table(segments))

  phat_seg <- .Call("bino_div_R", K1_seg, K1_seg + K0_seg)
  logL_seg <- -(.Call("bino_xlogy_R", K1_seg, phat_seg) + .Call("bino_xlogy_R", K0_seg, 1 - phat_seg))
  LL2 <- rowSums(logL_seg)
  df2 <- rep(N, n_segments)

  # --- Compute p-values and stability flag (using base R) ---
  pvals <- 1 - pchisq(-2 * (LL1 - LL2), df = df1 - df2)
  stable_flag <- pvals < 0.05

  # --- Compute average beta per segment ---
  total_K0 <- rowSums(K0_seg, na.rm = TRUE)
  total_K1 <- rowSums(K1_seg, na.rm = TRUE)
  beta_seg <- .Call("bino_div_R", total_K1, total_K0 + total_K1)

  # --- Compute genomic summaries ---
  segment_start <- tapply(pos, segments, min)
  segment_end   <- tapply(pos, segments, max)
  segment_chr   <- tapply(chr, segments, unique)
  segment_len   <- segment_end - segment_start + 1
  segment_cpgs  <- as.numeric(table(segments))

  # --- Combine into final summary data frame ---
  summary_df <- data.frame(
    Segment = as.integer(names(segment_start)),
    Chr = unlist(segment_chr),
    Start = unlist(segment_start),
    End = unlist(segment_end),
    CpGs = segment_cpgs,
    Length = segment_len,
    Beta = beta_seg,
    Stable = stable_flag
  )

  rownames(summary_df) <- NULL

  # --- Return both summary and per-segment betas ---
  return(list(
    summary = summary_df,
    betas_per_segment = phat_seg
  ))
}



