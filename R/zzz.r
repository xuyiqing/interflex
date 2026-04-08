.onAttach <-
   function(libname, pkgname) {
      # echo output to screen
      packageStartupMessage("## Major upgrade since v.1.3.0\n")
      packageStartupMessage("## See http://bit.ly/interflex for more info\n## Comments and suggestions -> yiqingxu@stanford.edu\n")
   }

## Silence R CMD check NOTEs about non-standard-evaluation symbols used
## inside aes() / dplyr / subset() calls. These are column names resolved
## at evaluation time, not package-level objects.
utils::globalVariables(c(
    ".data",
    "CI_uniform_lower",
    "CI_uniform_upper",
    "excluded.iv"
))
