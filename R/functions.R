#' A function that claculates the mean and standard deviation for each metabolite and rounds them to 1 decimal.
#'
#' @param data
#'
#' @return A data.frame/tibble.
descriptive_stats <- function(data) {
    data |>
        dplyr::group_by(metabolite) |>
        dplyr::summarise(dplyr::across(
            value,
            list(
                mean = mean,
                sd = sd
            )
        )) |>
        dplyr::mutate(dplyr::across(
            tidyselect::where(is.numeric),
            ~ round(.x, digits = 1)
        ))
}
