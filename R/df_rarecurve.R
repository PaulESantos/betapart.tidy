#' df_rarecurve
#'
#' @param comm Community data, a matrix.
#'
#' @return A data_frame with the output of function **rarefy**.
#' @export
#'
#' @examples
#' require(vegan)
#' require(tidyverse)
#'
#' data("dune")
#' data("BCI")
#'
#' df_rarecurve(dune)
#'
#' df_rarecurve(BCI) %>%
#' filter(sites == 10)
#'
df_rarecurve <- function(comm){

  rare <- function (comm, step = 1, sample)
  {
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
    dplyr::mutate(size = seq(1, length(sites), 1))
  return(output)
}
plot_rarecurve(dune)
