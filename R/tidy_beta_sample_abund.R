#' tidy_beta_sample_abund
#'
#' @param x matris de abundancia de especies.
#' @param index.family "bray" o "ruzicka".
#' @param sites  number of sites to be tested.
#' @param samples number of replicates.
#' @return
#' @export
#'
#' @examples
tidy_beta_sample_abund <- function (x, index.family = "bray", sites = nrow(x), samples = 1)
{
  requireNamespace("tidyverse")
  requireNamespace("vegan")
  requireNamespace("betapart")

  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  if (sites > nrow(x))
    stop("More sites requested for sample than are in the dataset")
  pb <- txtProgressBar(min = 0, max = samples, style = 3)
  results.n <- as.data.frame(matrix(nrow = samples, ncol = 3))
  for (i in 1:samples) {
    position <- as.vector(1:nrow(x))
    sample.position <- sample(position, sites)
    x.sample <- x[sample.position, ]
    x.beta <- betapart::beta.multi.abund(x.sample, index.family)
    results.n[i, ] <- unlist(x.beta)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(results.n) <- names(x.beta)
  result <- dplyr::as_data_frame(results.n)
  return(result)
}
