#' tidy_beta_temp
#' @param x matris de abundancia inicial.
#' @param y matris de abundancia a contrastar.
#' @param index.family "jaccard" o "sorensen".
#'
#' @return
#' @export
#'
#' @examples

tidy_beta_temp <- function (x, y, index.family = "sorensen")
{
  requireNamespace("tidyverse")
  requireNamespace("betapart")
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  if (!inherits(x, "betapart")) {
    x <- betapart::betapart.core(x)
  }
  if (!inherits(y, "betapart")) {
    y <- betapart::betapart.core(y)
  }
  if (!identical(dim(x$data), dim(y$data)))
    stop("The two data matrices do not have the same dimensions.")
  if (!identical(rownames(x$data), rownames(y$data)))
    stop("The two data matrices do not have the same site names.")
  if (!identical(colnames(x$data), colnames(y$data)))
    stop("The two data matrices do not have the same species names .")
  ai <- apply(x$data & y$data, 1, sum)
  bi <- apply(x$data & !y$data, 1, sum)
  ci <- apply(!x$data & y$data, 1, sum)
  switch(index.family, sorensen = {
    beta.sor <- (bi + ci)/(2 * ai + bi + ci)
    beta.sim <- pmin(bi, ci)/(ai + pmin(bi, ci))
    beta.sne <- beta.sor - beta.sim
    result <- data.frame(beta.sim, beta.sne, beta.sor) %>%
     dplyr::rownames_to_column("site") %>%
      dplyr::as_data_frame()
  }, jaccard = {
    beta.jac <- (bi + ci)/(ai + bi + ci)
    beta.jtu <- 2 * pmin(bi, ci)/(ai + (2 * pmin(bi, ci)))
    beta.jne <- beta.jac - beta.jtu
    result <- data.frame(beta.jtu, beta.jne, beta.jac) %>%
      dplyr::rownames_to_column("site") %>%
      dplyr::as_data_frame()
  })
  return(result)
}
