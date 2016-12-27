#' @title Recode levels
#'
#' todo.
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
    newOrder <- data[, mean(y), by = c] %>% arrange(V1)
    c <- factor(c, levels = as.character(newOrder[, 1]))
    c
}
