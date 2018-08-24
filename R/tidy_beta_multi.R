#' tidy_beta_multi
#'
#' @param x matriz de laabundancia de especies.
#' @param index.family pude trabajar con "sorensen" o "jaccard".
#' @return A a data_frame
#' @export


tidy_beta_multi <- function (x, index.family = "sorensen")
{
  requireNamespace("tidyverse")
  requireNamespace("vegan")
  requireNamespace("betapart")


  x <- vegan::decostand(x, method = "pa")
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  if (!inherits(x, "betapart")) {
    x <- betapart::betapart.core(x)
  }
  maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
  minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
  switch(index.family, sorensen = {
    beta.sim <- minbibj/(minbibj + x$a)
    beta.sne <- (x$a/(minbibj + x$a)) * ((maxbibj - minbibj)/((2 *
                                                                 x$a) + maxbibj + minbibj))
    beta.sor <- (minbibj + maxbibj)/(minbibj + maxbibj +
                                       (2 * x$a))
    multi <-dplyr::data_frame(beta.SIM = beta.sim, beta.SNE = beta.sne,
                        beta.SOR = beta.sor)
  }, jaccard = {
    beta.jtu <- (2 * minbibj)/((2 * minbibj) + x$a)
    beta.jne <- (x$a/((2 * minbibj) + x$a)) * ((maxbibj -
                                                  minbibj)/((x$a) + maxbibj + minbibj))
    beta.jac <- (minbibj + maxbibj)/(minbibj + maxbibj +
                                       x$a)
    multi <- dplyr::data_frame(beta.JTU = beta.jtu, beta.JNE = beta.jne,
                        beta.JAC = beta.jac)
  })
  return(multi)
}
