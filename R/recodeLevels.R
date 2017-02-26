#' @title Recode levels
#'
#'
#' @param y numeric vector
#' @param c factor vector (same length as y).
#'
#' @rdname recodeLevels
#'
#' @export
#'

recodeLevels <- function(y, c) {
    data <- data.table::data.table(y = y, c = c)
    newOrder <- data[, mean(y), by = c] %>% splyr::arrange(V1)
    c <- factor(c, levels = as.character(newOrder[, 1]))
    cat("fdsf")
    c
}
