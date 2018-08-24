#' tidy_phylo_beta_pair
#'
#' @param x matrix species.
#' @param tree dendrogram.
#' @param index.family "jaccard" or "sorensen".
#'
#' @return data_frame
#' @export
#'
#' @examples
tidy_phylo_beta_pair <- function (x, tree, index.family = "sorensen")
{
  requireNamespace("tidyverse")
  requireNamespace("corrr")
  requireNamespace("betapart")

  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  pbc <- x
  if (!inherits(x, "phylo.betapart")) {
    pbc <- betapart::phylo.betapart.core(x, tree)
  }
  switch(index.family, sorensen = {
    phylo.beta.sim <- pbc$min.not.shared/(pbc$min.not.shared +
                                            pbc$shared)
    phylo.beta.sne <- ((pbc$max.not.shared - pbc$min.not.shared)/((2 *
                                                                     pbc$shared) + pbc$sum.not.shared)) * (pbc$shared/(pbc$min.not.shared +
                                                                                                                         pbc$shared))
    phylo.beta.sor <- pbc$sum.not.shared/(2 * pbc$shared +
                                            pbc$sum.not.shared)
    phylo.pairwise <- dplyr::bind_cols(
      phylo.beta.sim = phylo.beta.sim %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      phylo.beta.sne = phylo.beta.sne %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      phylo.beta.sor = phylo.beta.sor %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"))%>%
      dplyr::select(x , y,
                    phylo.beta.sim = r,
                    phylo.beta.sne = r1,
                    phylo.beta.sor = r2)
  }, jaccard = {
    phylo.beta.jtu <- (2 * pbc$min.not.shared)/((2 * pbc$min.not.shared) +
                                                  pbc$shared)
    phylo.beta.jne <- ((pbc$max.not.shared - pbc$min.not.shared)/(pbc$shared +
                                                                    pbc$sum.not.shared)) * (pbc$shared/((2 * pbc$min.not.shared) +
                                                                                                          pbc$shared))
    phylo.beta.jac <- pbc$sum.not.shared/(pbc$shared + pbc$sum.not.shared)
    phylo.pairwise <- dplyr::bind_cols(
      phylo.beta.jtu = phylo.beta.jtu %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      phylo.beta.jne = phylo.beta.jne %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA"),
      phylo.beta.jac = phylo.beta.jac %>%
        as.matrix() %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
      dplyr::select(x, y,
                    phylo.beta.jtu = r,
                    phylo.beta.jne = r1,
                    phylo.beta.jac = r2 )
  })
  return(phylo.pairwise)
}
