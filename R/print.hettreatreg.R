#' @method print hettreatreg
#' @export
# Define new print function
#-------------------------------------------------------------------------------
print.hettreatreg <- function(x, ...) {

  # ---- Display results
  cat('\n')
  cat('"OLS" is the estimated regression coefficient on', paste(x$treatment_name, '.', sep = ""), '\n', '\n')
  cat('   OLS =', signif(x$OLS1, 4), '\n', '\n')
  cat('P(d=1) =', round(x$'P(d=1)', 3), '\n')
  cat('P(d=0) =', round(x$'P(d=0)', 3), '\n', '\n')
  cat('    w1 =', round(x$w1, 3), '\n')
  cat('    w0 =', round(x$w0, 3), '\n')
  cat(' delta =', round(x$delta, 3), '\n', '\n')
  cat('   ATE =', signif(x$ATE, 4), '\n')
  cat('   ATT =', signif(x$ATT, 4), '\n')
  cat('   ATU =', signif(x$ATU, 4), '\n', '\n')
  cat('OLS = w1*ATT + w0*ATU =', signif(x$OLS2, 4), '\n', '\n')
}
