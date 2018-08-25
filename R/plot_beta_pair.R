#' plot_bet_pair
#'
#' Representacion grafica de las matrices de disimilaridad resultado de
#' la funcion "beta.pair". Visualizacion grafica de Turnover y Nestedness.
#'
#'
#' @param x matriz de p/a o abundancia de especies.
#' @param index.family pude trabajar con "sorensen" o "jaccard".
#' @return A plot.
#' @export
#'
#'
plot_beta_pair <- function(x, index.family= "sorensen" ){
  requireNamespace("magrittr")
  requireNamespace("betapart")
  requireNamespace("tidyverse")
  requireNamespace("vegan")
  requireNamespace("lazyeval")
  #devtools::use_package("tidyverse")
  #devtools::use_package("forcats")



  dist_plot<- function(rdf,
                       legend = TRUE,
                       shape = 16,
                       colours = c("white", "black"),
                       print_cor = FALSE,
                       colors) {
    requireNamespace("magrittr")
    requireNamespace("tidyverse")
    requireNamespace("lazyeval")
    if (!missing(colors))
      colours <- colors

    # Store order for factoring the variables
    row_order <- rdf$rowname

    # Prep dots for mutate_
    dots <- stats::setNames(list(lazyeval::interp(~ factor(x, levels = row_order),
                                                  x = quote(x)),
                                 lazyeval::interp(~ factor(y, levels = rev(row_order)),
                                                  y = quote(y)),
                                 lazyeval::interp(~ abs(r),
                                                  r = quote(r)),
                                 lazyeval::interp(~ as.character(corrr::fashion(r)),
                                                  r = quote(r))
    ),
    list("x", "y", "size", "label"))

    # Convert data to relevant format for plotting
    pd <- rdf %>%
      # Convert to wide
      corrr::stretch(na.rm = TRUE) %>%
      # Factor x and y to correct order
      # and add text column to fill diagonal
      # See dots above
      dplyr::mutate_(.dots = dots)

    plot_ <- list(
      # Geoms
      ggplot2::geom_point(shape = shape),
      if (print_cor) ggplot2::geom_text(color = "black", size = 3, show.legend = FALSE),
      ggplot2::scale_colour_gradientn(limits = c(0, 1), colors = colours),
      # Theme, labels, and legends
      ggplot2::theme_classic(),
      ggplot2::labs(x = "", y =""),
      ggplot2::guides(size = "none", alpha = "none"),
      if (legend)  ggplot2::labs(colour = NULL),
      if (!legend) ggplot2::theme(legend.position = "none")
    )

    ggplot2::ggplot(pd, ggplot2::aes_string(x = "x", y = "y", color = "r",
                                            size = "size", alpha = "size",
                                            label = "label")) +
      plot_
  }

#data procesing
  x <- vegan::decostand(x, method = "pa")
  df <- betapart::beta.pair(x, index.family = index.family)
##################################################################################
  plot_a <- df[[1]] %>%
    as.matrix() %>%
    corrr::as_cordf() %>%
    corrr::shave() %>%
    dist_plot()+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 0,
                                               vjust = 0.5, size = 14),
                   axis.text.y.left  = ggplot2::element_text(angle = 0,
                                                    vjust = 0.5, size = 14),
                   axis.title.y = ggplot2::element_text(size = 18),
                   axis.title.x = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = .5,
                                             face = "bold", size = 15))
##################################################################################
  plot_b <- df[[2]] %>%
    as.matrix() %>%
    corrr::as_cordf() %>%
    corrr::shave() %>%
    dist_plot()+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 0,
                                                        vjust = 0.5, size = 14),
                   axis.text.y.left  = ggplot2::element_text(angle = 0,
                                                             vjust = 0.5, size = 14),
                   axis.title.y = ggplot2::element_text(size = 18),
                   axis.title.x = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(hjust = .5,
                                                      face = "bold", size = 15))
################################################################################
 cowplot::plot_grid(plot_a, plot_b, labels = c("Turnover", "Nestedness"))

  }



