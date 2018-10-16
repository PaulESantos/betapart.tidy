#' Title
#'
#' @param comm
#' @param method
#'
#' @return
#' @export
#'
#' @examples
plot_rankabund <- function(comm, method = "abundance"){
  rang_abun <- function(x){

    SPLIT <- c("none", "sites")
    if (is.na(pmatch(group, SPLIT)) | pmatch(group,
                                             SPLIT) == -1)
      stop("invalid group variable")
    group <- match.arg(group, SPLIT)

    df <- x %>%
      dplyr::as_data_frame() %>%
      tibble::rownames_to_column("sites") %>%
      tidyr::gather(species, abundance, -sites)

    output <- df %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(abundance = sum(abundance)) %>%
      dplyr::arrange(desc(abundance)) %>%
      dplyr::mutate(
        rank = seq(1, length(species)),
        proportion = (abundance / sum(abundance)) * 100,
        acumfreq = cumsum(proportion),
        logabun = log10(abundance)) %>%
      dplyr::ungroup()

      return(output)

  }

  SPLIT <- c("abundance", "logabun")
  if (is.na(pmatch(method, SPLIT)) | pmatch(method,
                                            SPLIT) == -1)
    stop("invalid method")
  method <- match.arg(method, SPLIT)


  df <-  rankabund(comm)

  themes <- ggplot2::theme(
    text = element_text(face = "bold", size = 15),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(size = rel(0.5)),
    strip.background = element_rect(fill = "grey85", colour = "grey20"))


  output <- df %>%
    ggplot2::ggplot(aes(rank, abundance))+
    ggplot2::geom_line(color = "grey", size = 1.2)+
    ggplot2::geom_point( size = 2)+
    ggplot2::labs(x = "Species Rank",
                  y = "Abundance")+
    themes+
    ggplot2::annotate("text",
                      x = as.vector(df$rank[1:5]),
                      y = as.vector(df$abundance[1:5]),
                      label = as.vector(df$species[1:5]),
                      hjust = -.2)

  if(method == "abundance"){
    return(output)
  }
  else(method == "logabun")
  {
    return(df %>%
             ggplot2::ggplot(aes(rank, logabun))+
             ggplot2::geom_line(color = "grey", size = 1.2)+
             ggplot2::geom_point( size = 2)+
             ggplot2::labs(x = "Species Rank",
                           y = "Log10(Abundance)")+
             themes)
  }
}
