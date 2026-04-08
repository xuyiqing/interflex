## print.interflex
##
## When an interflex object carries a $figure ggplot, auto-print it. This
## preserves the historical behavior of `interflex(estimator = "raw", ...)`,
## which used to return a bare ggplot directly and could be auto-printed at
## the REPL or in a knitr chunk. With raw now wrapped in an interflex object,
## the print method routes back to the figure so vignettes that rely on
## auto-print continue to work without modification.
##
## For interflex objects that do not carry a $figure (e.g. estimator-only
## fits where the user requested no figure construction), fall through to
## the default printer.

#' @export
print.interflex <- function(x, ...) {
    if (!is.null(x$figure)) {
        print(x$figure, ...)
        invisible(x)
    } else {
        cat("interflex object (", paste(class(x), collapse = "/"), ")\n", sep = "")
        if (!is.null(x$type)) cat("  type:", x$type, "\n")
        invisible(x)
    }
}
