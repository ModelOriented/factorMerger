#' Generated sample visualization
#'
#' @param x generatesSample object
#'
#' @rdname plot.generatedSample
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
             subtitle = "Groups are sorted by means. Displayed summary statistic - mean in a group.")
}

#' Plot heatmap
#'
#' @param x generatesSample object
#'
#' @importFrom ggplot2 ggplot geom_tile aes ylab xlab stat_summary labs theme_bw
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#'

plotHeatmap <- function(factorMerger) {
    factorMerger$factor <- factor(factorMerger$factor, levels = getFinalOrder(factorMerger))
    data.frame(cbind(response = factorMerger$response), factor = factorMerger$factor) %>%
        melt(id.vars = "factor") %>% ggplot() +
        geom_tile(aes(x = factor, y = variable, fill = value)) +
        theme_bw()
}

#' @importFrom ggtree read.tree ggtree gheatmap theme_tree2 geom_tiplab scale_x_ggtree
plotTree <- function(factorMerger) {
    tr <- getTree(factorMerger)
    df <- data.frame(factorMerger$factor, factorMerger$response)
    df <- aggregate(df[, -1], list(df[, 1]), mean)
    rownames(df) <- df[, 1]
    df <- df[, -1]
    (ggtree(read.tree(text = tr)) + geom_tiplab(align = TRUE) + theme_tree2()) %>%
        gheatmap(df, offset = 4, width = 0.5, colnames = FALSE) %>%
        scale_x_ggtree()
}

plot.subsequentFactorMerger <- function(fm) {
    fmMeans <- means(fm)
    fmList <- NULL
}
