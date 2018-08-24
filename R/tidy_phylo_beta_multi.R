#' tidy_phylo_beta_multi
#' @param x species matrix.
#' @param tree dendrogram
#' @param index.family "jaccard" or "sorensen".
#'
#' @return a data_frame.
#' @export
#'
#' @examples
tidy_phylo_beta_multi <- function (x, tree, index.family = "sorensen")
{ requireNamespace("betapart")
  requireNamespace("tidyverse")
  requireNamespace("corrr")

  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  pbc <- x
  if (!inherits(x, "phylo.betapart")) {
    pbc <- betapart::phylo.betapart.core(x, tree)
  }
  switch(index.family, sorensen = {
    phylo.beta.SIM <- sum(pbc$min.not.shared)/(pbc$sumSi -
                                                 pbc$St + sum(pbc$min.not.shared))
    phylo.beta.SNE <- ((sum(pbc$max.not.shared) - sum(pbc$min.not.shared))/(2 *
                                                                              (pbc$sumSi - pbc$St) + sum(pbc$min.not.shared) +
                                                                              sum(pbc$max.not.shared))) * ((pbc$sumSi - pbc$St)/(pbc$sumSi -
                                                                                                                                   pbc$St + sum(pbc$min.not.shared)))
    phylo.beta.SOR <- (sum(pbc$min.not.shared) + sum(pbc$max.not.shared))/(2 *
                                                                             (pbc$sumSi - pbc$St) + sum(pbc$min.not.shared) +
                                                                             sum(pbc$max.not.shared))
    phylo.multi <- dplyr::data_frame(phylo.beta.SIM = phylo.beta.SIM,
                                     phylo.beta.SNE = phylo.beta.SNE,
                                     phylo.beta.SOR = phylo.beta.SOR)
  }, jaccard = {
    phylo.beta.JTU <- (2 * sum(pbc$min.not.shared))/((2 *
                                                        sum(pbc$min.not.shared)) + pbc$sumSi - pbc$St)
    phylo.beta.JNE <- ((sum(pbc$max.not.shared) - sum(pbc$min.not.shared))/(pbc$sumSi -
                                                                              pbc$St + sum(pbc$max.not.shared) + sum(pbc$min.not.shared))) *
      ((pbc$sumSi - pbc$St)/(2 * sum(pbc$min.not.shared) +
                               pbc$sumSi - pbc$St))
    phylo.beta.JAC <- (sum(pbc$min.not.shared) + sum(pbc$max.not.shared))/(pbc$sumSi -
                                                                             pbc$St + sum(pbc$min.not.shared) + sum(pbc$max.not.shared))
    phylo.multi <- dplyr::data_frame(phylo.beta.JTU = phylo.beta.JTU,
                                     phylo.beta.JNE = phylo.beta.JNE,
                                     phylo.beta.JAC = phylo.beta.JAC)
  })
  return(phylo.multi)
}
