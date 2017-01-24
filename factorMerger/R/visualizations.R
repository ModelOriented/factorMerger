#' Visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot
#'

plot.generatedSample <- function(data) {
    data %>% ggplot2::ggplot() + ggplot2::geom_boxplot(ggplot2::aes(y = numericVec,
                                               x = factorVec,
                                               group = factorVec)) +
        ggplot2::xlab("groups") + ggplot2::ylab("numeric")
}
