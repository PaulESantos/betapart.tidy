#' tidy_beta_pair
#'
#'
#' @param x matriz de p/a o abundancia de especies.
#' @param index.family pude trabajar con "jaccard" o "sorensen".
#' @return A data_frame.
#' @export
#'
tidy_beta_pair <- function (x, index.family = "sorensen")
{ requireNamespace("tidyverse")
  requireNamespace("vegan")
  requireNamespace("corrr")
  requireNamespace("stats")

  x <- vegan::decostand(x, method = "pa")

  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  if (!inherits(x, "betapart")) {
    x <- betapart::betapart.core(x)
  }
  switch(index.family, sorensen = {
    beta.sim <- x$min.not.shared/(x$min.not.shared + x$shared)
    beta.sne <- ((x$max.not.shared - x$min.not.shared)/((2 *
                                                           x$shared) + x$sum.not.shared)) * (x$shared/(x$min.not.shared +
                                                                                                         x$shared))
    beta.sor <- x$sum.not.shared/(2 * x$shared + x$sum.not.shared)
    pairwise <- dplyr::bind_cols(beta.sim = stats::as.dist(beta.sim) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.sne = stats::as.dist(beta.sne) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.sor = stats::as.dist(beta.sor) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
      dplyr::select(x, y,beta.sim = r, beta.sne = r1, beta.sor= r2)
  }, jaccard = {
    beta.jtu <- (2 * x$min.not.shared)/((2 * x$min.not.shared) +
                                          x$shared)
    beta.jne <- ((x$max.not.shared - x$min.not.shared)/(x$shared +
                                                          x$sum.not.shared)) * (x$shared/((2 * x$min.not.shared) +
                                                                                            x$shared))
    beta.jac <- x$sum.not.shared/(x$shared + x$sum.not.shared)
    pairwise <- dplyr::bind_cols(beta.jtu = stats::as.dist(beta.jtu) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.jne = stats::as.dist(beta.jne) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.jac = stats::as.dist(beta.jac)%>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA")
    ) %>%
      dplyr::select(x, y,beta.jtu = r, beta.jne = r1, beta.jac= r2)
  })
  return(pairwise)
}
