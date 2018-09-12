#' plot_rarecurve
#'
#' @param comm Community data, a matrix.
#'
#' @return A plot class **ggplot**.
#' @export
#'
#' @examples
#' require(vegan)
#' data("dune")
#'
#' plot_rarecurve(dune)
#' Manipulate the output with tidyverse.
#' plot_rarecurve(dune) +
#'   facet_wrap(.~sites) +
#'   theme(legend.position = "none")

plot_rarecurve <- function(comm){

  rare <- function (comm, step = 1, sample) {
    x <- as.matrix(comm)
    if (!identical(all.equal(x, round(x)), TRUE))
      stop("function accepts only integers (counts)")
    tot <- rowSums(x)
    S <- vegan::specnumber(x)
    if (any(S <= 0)) {
      message("empty rows removed")
      x <- x[S > 0, , drop = FALSE]
      tot <- tot[S > 0]
      S <- S[S > 0]
    }
    nr <- nrow(x)
    out <- lapply(seq_len(nr), function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i])
        n <- c(n, tot[i])
      drop(vegan::rarefy(x[i, ], n))
    })
    return(out)
  }

  dat <- rare(comm)

  output <- dat %>%
    reshape2::melt() %>%
    dplyr::as_data_frame() %>%
    dplyr::rename( "species" = "value","sites" = "L1") %>%
    dplyr::mutate(sites = as.character(sites)) %>%
    dplyr::group_by(sites) %>%
    dplyr::mutate(size = seq(1, length(sites), 1)) %>%
    ggplot2::ggplot(aes(size, species, color = sites))+
    ggplot2::geom_line(aes(group = sites), size = 1)+
    ggplot2::theme_bw()+
    ggplot2::scale_y_continuous(breaks = seq(1, 100, 1))+
    ggplot2::theme(legend.position = "bottom",
                   strip.text.x = element_text(size = 18, face = "bold"),
                   strip.background = element_rect(fill = "white"),
                   text = element_text(face = "bold", size = 10),
                   legend.text = element_text(face = "bold", size = 12))+
    ggplot2::labs(x = "Sample Size",
                  y = "Species",
                  color = "Sites") +
    ggplot2::guides(color = guide_legend(ncol = 10))
  class(output) <- c("gg", "ggplot")
  return(output)
}
