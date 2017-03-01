#' Set increasing order
#'
#' Given a numeric vector and a factor vector changes levels order
#' in the factor. Means of the numeric variable in groups
#' appointed by the new order of the factor are increasing.
#'
#' @param numericVec numeric vector
#' @param factorVec factor vector (same length as numbericVec).
#'
#' @rdname setIncreasingOrder
#'

setIncreasingOrder <- function(numericVec, factorVec) {
    if (length(numericVec) != length(factorVec)) {
        stop("Vectors' lengths do not match.")
    }
    stopifnot(is.vector(numericVec))
    stopifnot(is.vector(numericVec))
    factorVec <- as.factor(as.character(factorVec))
    data <- data.table::data.table(y = numericVec, c = factorVec)
    newOrder <- aggregate(numericVec ~ c, mean, data = data) %>%
        dplyr::arrange(numericVec)
    factor(factorVec, levels = as.character(newOrder[, 1]))
}

#' Merge factor
#'
mergeLevels <- function(factor, groupA, groupB, groupAB) {
    whichLevels <- which(levels(factor) %in% c(groupA, groupB))
    levels(factor)[whichLevels] <- groupAB
    factor
}

#' Calculate means by factor
#'
#' @rdname calculateMeans
#' @importFrom dplyr group_by summarize arrange
#'
calculateMeans <- function(numericVec, factorVec) {
    if (!is.null(dim(numericVec))) {
        return(NA)
    }
    data.frame(num = numericVec, level = factorVec) %>%
        group_by(level) %>% summarize(mean = mean(num)) %>%
        arrange(mean)
}

#' Filter groups - ...
#'
filterGroups <- function(response, factor, groupA, groupB) {
    response <- as.matrix(response)
    xA <- response[factor == groupA, ]
    xB <- response[factor == groupB, ]
    return(list(xA, xB))
}

#' Get final order - ...
#'
getFinalOrder <- function(factorMerger) {
    lastLevel <- factorMerger$mergingList[[length(factorMerger$mergingList)]]$groups
    splitted <- strsplit(lastLevel, split = "\\(|\\)|\\,")[[1]]
    splitted[nchar(splitted) > 0]
}

bindLevels <- function(groups, groupVec) {
    groupLabel <- paste(groupVec, sep = ":", collapse = ":")
    groups[groups %in% groupVec] <- groupLabel
    groups
}

getTree <- function(factorMerger) {
    steps <- length(factorMerger$mergingList)
    return(paste0(factorMerger$mergingList[[steps]]$groups, ";"))
}
