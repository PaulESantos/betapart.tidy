#' tidy_functional_beta_multi
#'
#' @param x matris de abundancia.
#' @param traits traits info data.frame
#' @param index.family "jaccard" o "sorensen".
#' @param warning.time progres bar.
#' @return data_frame
#' @export
#'
#' @examples
tidy_functional_beta_multi <- function (x, traits, index.family = "sorensen", warning.time = TRUE)
{
  requireNamespace("tidyverse")
  requireNamespace("betapart")

  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  fbc <- x
  if (!inherits(x, "functional.betapart")) {
    fbc <- betapart::functional.betapart.core(x, traits, multi = TRUE,
                                    warning.time = warning.time, return.details = FALSE)
  }
  maxbibj <- sum(fbc$max.not.shared[lower.tri(fbc$max.not.shared)])
  minbibj <- sum(fbc$min.not.shared[lower.tri(fbc$min.not.shared)])
  switch(index.family, sorensen = {
    funct.beta.sim <- minbibj/(minbibj + fbc$a)
    funct.beta.sne <- (fbc$a/(minbibj + fbc$a)) * ((maxbibj -
                                                      minbibj)/((2 * fbc$a) + maxbibj + minbibj))
    funct.beta.sor <- (minbibj + maxbibj)/(minbibj + maxbibj +
                                             (2 * fbc$a))
    functional.multi <- dplyr::data_frame(
      funct.beta.SIM = funct.beta.sim,
      funct.beta.SNE = funct.beta.sne,
      funct.beta.SOR = funct.beta.sor)
  }, jaccard = {
    funct.beta.jtu <- (2 * minbibj)/((2 * minbibj) + fbc$a)
    funct.beta.jne <- (fbc$a/((2 * minbibj) + fbc$a)) *
      ((maxbibj - minbibj)/((fbc$a) + maxbibj + minbibj))
    funct.beta.jac <- (minbibj + maxbibj)/(minbibj + maxbibj +
                                             fbc$a)
    functional.multi <- dplyr::data_frame(
      funct.beta.JTU = funct.beta.jtu,
      funct.beta.JNE = funct.beta.jne,
      funct.beta.JAC = funct.beta.jac)
  })
  return(functional.multi)
}
