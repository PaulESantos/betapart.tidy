#' plot_rarecurve
#' @description
#'
#'   Draws a rarefaction curve for each row of the input data.
#' @param comm Community data.
#' @param facet.var Index for changing the output plot:
#'  "none" the default output.
#'  "sites" the output is a multi-panel plot.
#'
#' @return
#' @export
#'
#' @examples
#' require(vegan)
#' data(dune)
#' plot_rarecurve(dune)
#' plot_rarecurve(dune, facet.var = "sites")
#'
#'
plot_rarecurve <-  function (comm, facet.var = "none")
  {
    SPLIT <- c("none","sites")
    if (is.na(pmatch(facet.var, SPLIT)) | pmatch(facet.var,
                                                 SPLIT) == -1)
      stop("invalid facet variable")
    facet.var <- match.arg(facet.var, SPLIT)

    rare <- function(comm, step = 1, sample) {
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
    meta <- dplyr::data_frame(name = rownames(comm),
                              sites = seq(1,length(name)))

    dat <- rare(comm)

    dat <- dat %>%
      reshape2::melt() %>%
      dplyr::as_data_frame() %>%
      dplyr::rename(species = "value", sites = "L1")

    dat1 <- meta %>%
      dplyr::full_join(dat, by = "sites")%>%
      dplyr::group_by(name) %>%
      dplyr::mutate(size = seq(1,length(sites), 1)) %>%
      dplyr::ungroup()

    output <- dat1  %>%
      ggplot2::ggplot(aes(size, species)) +
      ggplot2::geom_line(aes(group = name), size = 1, color = "grey") +
      ggplot2::theme_bw() +
      ggplot2::geom_text(data = dat1 %>% filter(size == last(size)),
                         aes(
                           label = name,
                           x = size + 0.5,
                           y = species,
                           color = sites
                         ), size = 5)+
      ggplot2::scale_y_continuous(breaks = seq(1, 1000, 10)) +
      ggplot2::theme(
        text = element_text(face = "bold", size = 10),
        axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 12)) +
      ggplot2::labs(x = "Sample Size", y = "Species") +
      ggplot2::guides(color = FALSE)


    if (facet.var == "none") {
      return(output)
    }
    else (facet.var == "sites")
    {
      output <- dat1 %>%
        ggplot2::ggplot(aes(size, species)) +
        ggplot2::geom_line(aes(group = name), size = .8, color = "green") +
        ggplot2::theme_bw() +
        ggplot2::scale_y_continuous(breaks = seq(1, 1000, 10)) +
        ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                       strip.text.x = element_text(size = 12, face = "bold"),
                       strip.background = element_rect(fill = "white"),
                       text = element_text(face = "bold", size = 10),
                       axis.text = element_text(face = "bold", size = 10),
                       panel.grid = element_blank()) +
        ggplot2::facet_wrap(.~name)+
        ggplot2::labs(x = "Sample Size", y = "Species", color = "Sites")+
        ggplot2::theme(legend.position = "none")
      return(output)
    }
    class(output) <- c("gg", "ggplot")
    return(output)
  }
