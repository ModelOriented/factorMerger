.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Welcome to FactorMerger!")
}

.onLoad <- function(libname, pkgname) {
    op <- options()
    op.devtools <- list(
        devtools.path = "~/R-dev",
        devtools.install.args = "",
        devtools.name = "Durszlaczek",
        devtools.desc.author = '"Agnieszka Sitko <ag.agnieszka.sitko@gmail.cpm> [aut, cre]"',
        devtools.desc.license = "What license is it under?",
        devtools.desc.suggests = NULL,
        devtools.desc = list()
    )
    toset <- !(names(op.devtools) %in% names(op))
    if(any(toset)) options(op.devtools[toset])

    invisible()
}
