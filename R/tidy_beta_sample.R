#' tidy_beta_sample
#'
#' @param x matris de abundancia o presencia/ausencia.
#' @param index.family "jaccard" o "sorensen"
#' @param sites	 number of sites for which multiple-site dissimilarities will be computed. If not specified, default is all sites.
#' @param samples	 number of repetitions. If not specified, default is 1.

tidy.beta.sample <- function (x, index.family = "sorensen", sites = nrow(x$data),
                              samples = 1)
{ requireNamespace("tidyverse")
  requireNamespace("vegan")
  requireNamespace("betapart")
  requireNamespace("utils")
  x <- vegan::decostand(x, method = "pa")
  if (!inherits(x, "betapart")) {
    x <- betapart::betapart.core(x)
  }
  pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  if (sites > nrow(x$data))
    stop("More sites requested for sample than are in the dataset")
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  results.n <- as.data.frame(matrix(nrow = samples, ncol = 3))
  subset.betapart.core <- function(bpc, rows) {
    bpc$data <- bpc$data[rows, ]
    bpc$shared <- bpc$shared[rows, rows]
    bpc$not.shared <- bpc$not.shared[rows, rows]
    bpc$sumSi <- sum(diag(bpc$shared))
    bpc$St <- sum(colSums(bpc$data) > 0)
    bpc$a <- bpc$sumSi - bpc$St
    bpc$sum.not.shared <- bpc$sum.not.shared[rows, rows]
    bpc$max.not.shared <- bpc$max.not.shared[rows, rows]
    bpc$min.not.shared <- bpc$min.not.shared[rows, rows]
    return(bpc)
  }
  for (i in 1:samples) {
    position <- as.vector(1:nrow(x$data))
    sample.position <- sample(position, sites)
    x.sample <- subset.betapart.core(x, sample.position)
    x.beta <- betapart::beta.multi(x.sample, index.family)
    results.n[i, ] <- unlist(x.beta)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  names(results.n) <- names(x.beta)
  result <- dplyr::as_data_frame(results.n)

  return(result)
}
