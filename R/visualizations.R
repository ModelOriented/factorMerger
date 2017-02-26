#' Visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot
#' @importFrom ggplot2 ggplot geom_boxplot aes ylab xlab stat_summary labs
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize left_join
#'
#' @export

plot.generatedSample <- function(data) {
    colnames(data) <- c("y", "x")
    data %>% left_join(
        data %>% group_by(x) %>% summarize(y0 = min(y),
                                           y25 = quantile(y, 0.25),
                                           y50 = mean(y),
                                           y75 = quantile(y, 0.75),
                                           y100 = max(y)), by = "x") %>%
        ggplot(aes(y = y, x = x, group = x)) +
        geom_boxplot(aes(ymin = y0,
                         lower = y25,
                         middle = y50,
                         upper = y75,
                         ymax = y100), stat = "identity") +

        xlab("Groups") + ylab("Response") +
        labs(title = "Generated sample - boxplot",
             subtitle = "Groups are sorted by means.
             Displayed summary statistic - mean in a group.")
}

plot.subsequentFactorMerger <- function(fm) {
    fmMeans <- means(fm)
    fmList <- NULL
}
