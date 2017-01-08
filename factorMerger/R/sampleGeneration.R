#' Generate sample
#'
#' Function generateSample() produces a random sample consisting
#' of k groups drawn from the same distribution with different
#' parameters.
#'
#' @param N sample size
#' @param k number of groups
#' @param distr type of distribution from \code{c("norm", "exp", "beta")}
#'
#' @examples
#' generateSample(100, 2)
#' generateSample(100, 2, "exp")
#'
#' @return
#' \code{data.table} with two columns: numeric variable
#' and factor variable.
#'
#' @rdname generateSample
#'
#'
#'
#'
generateSample <- function(N, k, distr = "norm") {
    numericVec <- switch(distr,
        "norm" = rnorm(N),
        "exp" = rexp(N, 1),
        "beta" = rbeta(N, 1, 1),
        stop("Unknown distribution."))

    factorVec <- as.factor(sample(LETTERS[1:k],
                              size = N,
                              replace = TRUE))

    for (i in 1:k) {
        randomShift <- sample(seq(0, 1, 0.1), size = 1)
        let <- LETTERS[i]
        numericVec[factorVec == let] <- numericVec[factorVec == let] + randomShift
    }

    generatedSample <-
        data.table::data.table(numericVec = numericVec,
                               factorVec = setIncreasingOrder(numericVec, factorVec))
    class(generatedSample) <- append("generatedSample", class(generatedSample))
    generatedSample
}
