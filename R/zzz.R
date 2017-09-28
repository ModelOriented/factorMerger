.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste0("Welcome to factorMerger ",
               packageVersion("factorMerger"), "!\n"),
        "In this version arguments: 'method' and 'successive' of
        the mergeFactors() function are merged in one argument ('method') with
        its new values:
        'fast-fixed' (method = 'hclust', successive = TRUE),
        'fixed' (method = 'hclust', successive = FALSE),
        'fast-adapive' (method = 'LRT', successive = TRUE),
        'adaptive' (method = 'LRT', successive = FALSE).")
}

.onLoad <- function(libname, pkgname) {
    op <- options()
    op.devtools <- list(
        devtools.path = "~/R-dev",
        devtools.install.args = "",
        devtools.name = "Durszlaczek",
        devtools.desc.author =
            '"Agnieszka Sitko <ag.agnieszka.sitko@gmail.cpm> [aut, cre]"',
        devtools.desc.license = "GPL",
        devtools.desc.suggests = NULL,
        devtools.desc = list()
    )
    toset <- !(names(op.devtools) %in% names(op))
    if(any(toset)) options(op.devtools[toset])

    invisible()
}
