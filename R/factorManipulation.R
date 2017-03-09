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
mergeLevels <- function(factor, groupA, groupB, groupAB = NULL) {
    if (is.null(groupAB)) {
        groupAB <- paste0(groupA, groupB)
    }
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
    df <- data.frame(num = numericVec, level = factorVec)
    aggregate(num ~ level, mean, data = df) %>%
        rename(mean = num) %>%
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

bindLevels <- function(groups, groupVec) {
    groupLabel <- paste(groupVec, sep = ":", collapse = ":")
    groups[groups %in% groupVec] <- groupLabel
    groups
}

getTree <- function(factorMerger) {
    steps <- length(factorMerger$mergingList)
    return(paste0(factorMerger$mergingList[[steps]]$groups, ";"))
}

#' @importFrom reshape2 melt
#' @importFrom dplyr rename
calculateMeansByFactor <- function(response, factor) {
    means <- apply(as.data.frame(response), 2, function(x) {
        aggregate(x ~ level, mean, data = data.frame(x = x, level = factor))
    })

    means <- lapply(means, function(x) {
        df <- x %>% arrange(x)
        df$rank <- ave(df$x, FUN = rank)
        df
    })

    return(melt(means, id.vars = c("level", "rank")) %>%
               subset(select = -variable) %>%
               rename(mean = value,
                      variable = L1))
}
