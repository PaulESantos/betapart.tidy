#' tidy_functional_beta_pair
#'
#' @param x species matrix
#' @param traits traits dat.frame
#' @param index.family "jaccard" o "sorensen".
#'
#' @return
#' @export
#'
#' @examples
tidy_functional_beta_pair <- function (x, traits, index.family = "sorensen")
{
  requireNamespace("tidyverse")
  requireNamespace("betapart")
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  fbc <- x
  if (!inherits(x, "functional.betapart")) {
    fbc <- betapart::functional.betapart.core(x, traits, multi = FALSE,
                                    warning.time = FALSE, return.details = FALSE)
  }
  switch(index.family, sorensen = {
    funct.beta.sim <- fbc$min.not.shared/(fbc$min.not.shared +
                                            fbc$shared)
    funct.beta.sne <- ((fbc$max.not.shared - fbc$min.not.shared)/((2 *
                                                                     fbc$shared) + fbc$sum.not.shared)) * (fbc$shared/(fbc$min.not.shared +
                                                                                                                         fbc$shared))
    funct.beta.sor <- fbc$sum.not.shared/(2 * fbc$shared +
                                            fbc$sum.not.shared)
    functional.pairwise <- dplyr::bind_cols(
      funct.beta.sim = as.dist(funct.beta.sim) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      funct.beta.sne = as.dist(funct.beta.sne) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      funct.beta.sor = as.dist(funct.beta.sor) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
      dplyr::select(x, y,
                    funct.beta.sim = r,
                    funct.beta.sne = r1,
                    funct.beta.sor= r2)
  }, jaccard = {
    funct.beta.jtu <- (2 * fbc$min.not.shared)/((2 * fbc$min.not.shared) +
                                                  fbc$shared)
    funct.beta.jne <- ((fbc$max.not.shared - fbc$min.not.shared)/(fbc$shared +
                                                                    fbc$sum.not.shared)) * (fbc$shared/((2 * fbc$min.not.shared) +
                                                                                                          fbc$shared))
    funct.beta.jac <- fbc$sum.not.shared/(fbc$shared + fbc$sum.not.shared)
    functional.pairwise <- dplyr::bind_cols(
      funct.beta.jtu = as.dist(funct.beta.jtu) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      funct.beta.jne = as.dist(funct.beta.jne) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      funct.beta.jac = as.dist(funct.beta.jac) %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
      dplyr::select(x, y,
                    funct.beta.jtu = r,
                    funct.beta.jne = r1,
                    funct.beta.jac= r2)
  })
  return(functional.pairwise)
}
