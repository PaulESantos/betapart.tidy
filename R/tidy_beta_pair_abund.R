
#' tidy_beta_pair_abund
#'
#' @param x matris de p/a o abundancia de especies.
#' @param index.family indices "bray" o "ruzicka"
#'
#' @return data_frame
#' @export
#'
#' @examples

tidy_beta_pair_abund <- function (x, index.family = "bray")
{ requireNamespace("tidyverse")
  requireNamespace("corrr")
  requireNamespace("betapart")
  index.family <- match.arg(index.family, c("bray", "ruzicka"))
  if (!inherits(x, "betapart.abund")) {
    x <- betapart::betapart.core.abund(x)
  }
  switch(index.family, bray = {
    beta.bray.bal <- x$pair.min.not.shared.abund/(x$pair.min.not.shared.abund +
                                                    x$pair.shared.abund)
    beta.bray.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/((2 *
                                                                                      x$pair.shared.abund) + x$pair.sum.not.shared.abund)) *
      (x$pair.shared.abund/(x$pair.min.not.shared.abund +
                              x$pair.shared.abund))
    beta.bray <- x$pair.sum.not.shared.abund/(2 * x$pair.shared.abund +
                                                x$pair.sum.not.shared.abund)
    pairwise <- dplyr::bind_cols(beta.bray.bal = as.dist(beta.bray.bal) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.bray.gra = as.dist(beta.bray.gra) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                           corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.bray = as.dist(beta.bray) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
      dplyr::select(x, y, beta.bray.bal = r, beta.bray.gra = r1, beta.bray= r2)
  }, ruzicka = {
    beta.ruz.bal <- (2 * x$pair.min.not.shared.abund)/((2 *
                                                          x$pair.min.not.shared.abund) + x$pair.shared.abund)
    beta.ruz.gra <- ((x$pair.max.not.shared.abund - x$pair.min.not.shared.abund)/(x$pair.shared.abund +
                                                                                    x$pair.sum.not.shared.abund)) * (x$pair.shared.abund/((2 *
                                                                                                                                             x$pair.min.not.shared.abund) + x$pair.shared.abund))
    beta.ruz <- x$pair.sum.not.shared.abund/(x$pair.shared.abund +
                                               x$pair.sum.not.shared.abund)
    pairwise <- dplyr::bind_cols(beta.ruz.bal = as.dist(beta.ruz.bal) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.ruz.gra = as.dist(beta.ruz.gra) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA"),
                          beta.ruz = as.dist(beta.ruz) %>%
                            as.matrix() %>%
                            corrr::as_cordf() %>%
                            corrr::shave() %>%
                            corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
     dplyr::select(x, y, beta.ruz.bal = r, beta.ruz.gra = r1, beta.ruz= r2)
  })
  return(pairwise)
}
