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

#' A histogram of metabolite distributions
#'
#' @param data
#'
#' @return A plot object
plot_distributions <- function(data) {
    ggplot2::ggplot(data, aes(x = value)) +
        ggplot2::geom_histogram() +
        ggplot2::facet_wrap(vars(metabolite), scales = "free")
}
