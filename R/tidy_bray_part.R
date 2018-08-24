#' tidy_bray_part
#'
#' @param x matris de abundancia de especies
#'
#' @return data_frame
#' @export
#'
#' @examples

tidy_bray_part <- function (x)
{requireNamespace("tidyverse")
  requireNamespace("corrr")
  x <- as.matrix(x)
  result <- matrix(nrow = nrow(x), ncol = nrow(x))
  rownames(result) <- rownames(x)
  colnames(result) <- rownames(x)
  for (i in 1:nrow(x)) {
    for (j in i:nrow(x)) {
      A <- sum(pmin(x[i, ], x[j, ]))
      B <- sum(x[i, ]) - sum(pmin(x[i, ], x[j, ]))
      C <- sum(x[j, ]) - sum(pmin(x[i, ], x[j, ]))
      result[i, j] <- min(B, C)/(A + min(B, C))
      result[j, i] <- (B + C)/(2 * A + B + C)
    }
  }
  bray <- as.dist(result)
  bray.bal <- as.dist(t(upper.tri(result) * result))
  bray.gra = bray - bray.bal
  results <- dplyr::bind_cols(bray.bal = bray.bal%>%
                         as.matrix() %>%
                         corrr::as_cordf() %>%
                        corrr::shave() %>%
                         corrr::stretch() %>% dplyr::filter(r != "NA"),
                       bray.gra = bray.gra%>%
                         as.matrix() %>%
                         corrr::as_cordf() %>%
                         corrr::shave() %>%
                        corrr::stretch() %>% dplyr::filter(r != "NA"),
                       bray = bray%>%
                         as.matrix() %>%
                         corrr::as_cordf() %>%
                         corrr::shave() %>%
                         corrr::stretch() %>% dplyr::filter(r != "NA")) %>%
    dplyr::select(x, y, bray.bal = r, bray.gra = r1, bray= r2)
  return(results)
}
