#' Visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot
#' @importFrom ggplot2 ggplot geom_boxplot aes ylab xlab stat_summary labs
#' @importFrom magrittr %>%
#'
#' @export

plot.generatedSample <- function(data) {
    data %>% ggplot(aes(y = numericVec,
                        x = factorVec,
                        group = factorVec)) +
        geom_boxplot() +
        xlab("groups") + ylab("numeric") +
        stat_summary(fun.y = mean, geom = "point", shape = 23, size = 3, fill = "violet") +
        labs(title = "Generated sample - boxplot",
             subtitle = "Groups sorted by their means (means - violet rectangles, medians - black lines)")
}
