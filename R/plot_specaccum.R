#' plot_specaccum
#'
#' @param comm  Community data set.
#' @param method The function use the same methods that specaccum function from
#' vegan package "random", "exact", "coleman", "rarefaction".
#' Except for "collector".
#' @param permutations Number of permutations with method = "random".
#'
#' @return A graph with ggplot2 characters.
#' @details For more info about the calculations checks help(specaccum).
#' @export
#'
#' @examples
#' require(vegan)
#' data("BCI")
#'
#' plot_specaccum(BCI)
#' plot_specaccum(BCI, method = "random", permutations = 1000)
#'
#' Modify the layers:
#'
#' plot_specaccum(BCI)+
#' theme_classic()+
#' theme(text = element_text(face = "bold"))
#'
plot_specaccum <- function(comm, method = "exact", permutations = 100){
  df <-vegan::specaccum(comm, method = method, permutations = permutations)

  plot <- dplyr::data_frame(sites = df$sites,
                            richness = df$richness,
                            sd = df$sd) %>%
    ggplot2::ggplot(aes(sites, richness))+
    ggplot2::geom_line(linetype = 1, size = .5, color = "red") +
    ggplot2::geom_errorbar(aes(ymin = richness - 2 * sd,
                               ymax = richness + 2 * sd),
                           width = 0, size = 1)+
    ggplot2::geom_point(color = "red", size = 1.5)+
    ggplot2::theme_bw()+
    ggplot2::scale_y_continuous(name = "Number of species")+
    ggplot2::scale_x_continuous(name = "Sites")+
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 0,
                                                        vjust = 0.5, size = 14),
                   axis.text.y.left  = ggplot2::element_text(angle = 0,
                                                             vjust = 0.5, size = 14),
                   axis.title.y = ggplot2::element_text(size = 18),
                   axis.title.x = ggplot2::element_text(size = 18),
                   legend.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size=12, face="bold"),
                   plot.title =ggplot2::element_text(hjust = .5,
                                                     face = "bold", size = 21),
                   legend.justification=c(0,1),
                   legend.position=c(.8, 1))+
    ggplot2::labs( color = "", fill="")+
    ggplot2::ggtitle("Species Accumulation Curve")

  return(plot)
}
