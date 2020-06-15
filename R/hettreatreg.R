#' OLS Weights on Heterogeneous Treatment Effects
#'
#' Computes diagnostics for linear regression when treatment effects are heterogeneous.
#'
#'@name hettreatreg
#'
#'@param outcome the outcome variable.
#'@param treatment the treatment variable. The variable must be binary and coded 0 for the untreated units and 1 for the treated units.
#'@param covariates the list of control variables. The list must not include the treatment variable.
#'@param verbose logical. If \code{TRUE} estimation output is displayed.
#'
#'@details \code{hettreatreg} represents ordinary least squares (OLS) estimates of the effect of a binary treatment as a weighted average of the average treatment effect on the treated (ATT) and the average treatment effect on the untreated (ATU). The program estimates the OLS weights on these parameters, computes the associated model diagnostics, and reports the implicit OLS estimate of the average treatment effect (ATE). See Sloczynski (2019) for the underlying theoretical results and further details.
#'
#'The arguments \code{outcome} and \code{treatment} are used to designate an outcome variable and a treatment variable, respectively. The treatment variable must be binary and coded 0 for the untreated units and 1 for the treated units. \code{covariates} is a list of control variables that must not include the treatment variable.
#'
#'\code{hettreatreg} displays a number of statistics. \code{OLS} is the estimated regression coefficient on the treatment variable. \code{P(d=1)} and \code{P(d=0)} are the sample proportions of treated and untreated units, respectively. \code{w1} and \code{w0} are the OLS weights on ATT and ATU, respectively. \code{delta} is a diagnostic for interpreting OLS as ATE. \code{ATE}, \code{ATT}, and \code{ATU} are the implicit OLS estimates of the corresponding parameters. See Sloczynski (2019) for further details.
#'
#'If you use this program in your work, please cite Sloczynski (2019).
#'
#'@return
#'  \item{\code{OLS}}{ OLS estimate of the treatment effect}
#'  \item{\code{P(d=1)}}{ proportion of treated units}
#'  \item{\code{P(d=0)}}{ proportion of untreated units}
#'  \item{\code{w1}}{ OLS weight on ATT}
#'  \item{\code{w0}}{ OLS weight on ATU}
#'  \item{\code{delta}}{ diagnostic for interpreting OLS as ATE}
#'  \item{\code{ATE}}{ implicit OLS estimate of ATE}
#'  \item{\code{ATT}}{ implicit OLS estimate of ATT}
#'  \item{\code{ATU}}{  implicit OLS estimate of ATU}
#' @export

#'@references Sloczynski, Tymon (2019). "Interpreting OLS Estimands When Treatment Effects Are Heterogeneous: Smaller Groups Get Larger Weights." Available at \url{http://people.brandeis.edu/~tslocz/Sloczynski_paper_regression.pdf}.
#'
#'@author Tymon Sloczynski, Brandeis University, \email{tslocz@@brandeis.edu}, \url{http://people.brandeis.edu/~tslocz/}
#'
#' Maintained by: Mark McAvoy, Brandeis University, \email{mcavoy@@brandeis.edu}
#'
#'Please feel free to report bugs and share your comments on this program.
#'
#'@examples
#' # load package
#' library(hettreatreg)
#'
#' # read in data
#' data("nswcps")
#'
#' # save the outcome variable
#' outcome <- nswcps$re78
#'
#' # save the treatment variable
#' treated <- nswcps$treated
#'
#' # select control variables
#' our_vars <- c("age", "age2", "educ", "black", "hispanic", "married", "nodegree")
#' covariates <- subset(nswcps, select = our_vars)
#'
#' # run function
#' results <- hettreatreg(outcome, treated, covariates)
#' print(results)
#'
#' @export
#'
#' @importFrom stats lm na.omit predict var
#-------------------------------------------------------------------------------

hettreatreg <- function(outcome, treatment, covariates, verbose = FALSE) {

  # is.binary will check if the argument is binary,
  # it has the added benefit of aggregating to true if any true is present in the list
  # while aggregating to false if all false in the list
  is.binary <- function(v) {
    x <- unique(v)
    length(x) - sum(is.na(x)) == 2L
  }

  # to check if treatment and outcome are the same, put into a data frame
  # and check if there are duplicates in this dataframe
  df.outcome.treatment <- data.frame(outcome, treatment)
  if (is.binary(duplicated(t(df.outcome.treatment))) == "TRUE") {
    stop("outcome and treatment must not be the same")
  }

  # similarly to check if outcome is in the set of control variables
  # place in the same data frame and apply the is.binary function
  df.outcome.covariates <- data.frame(covariates, outcome)
  if (is.binary(duplicated(t(df.outcome.covariates))) == "TRUE" ) {
    stop("outcome must not appear on the list of control variables")
  }

  df.treatment.covariates <- data.frame(covariates, treatment)
  if (is.binary(duplicated(t(df.treatment.covariates))) == "TRUE" ) {
    stop("treatment must not appear on the list of control variables")
  }

  if (is.binary(treatment) != TRUE) {
    stop("treatment must be a binary variable")
  }

  if (min(treatment, na.rm = TRUE) != 0 || max(treatment, na.rm = TRUE) != 1) {
    stop("treatment must only take on values zero or one")
  }

  treatment_name <- deparse(substitute(treatment)) # save name for printing

  # ---- Subset data to remove NAs
  df_first <- data.frame(outcome = outcome,
                         treatment = treatment,
                         covariates = covariates)
  df_first2 <- na.omit(df_first)

  covariates_names <- setdiff(names(df_first2), c("outcome", "treatment"))
  covariates <- subset(df_first2, select = covariates_names)

  outcome <- df_first2$outcome; treatment <- df_first2$treatment
  covariates <- as.matrix(covariates)

  # ---- Main regression
  m1 <- lm(outcome ~ treatment + covariates)
  beta_treated <- as.numeric(m1$coefficients['treatment'])

  # ---- Calculate treatment effects
  m_propensity <- lm(treatment ~ covariates)
  ps <- predict(m_propensity)

  df_combined <- data.frame(outcome, treatment, covariates, ps)
  df.1 <- df_combined[which(treatment == 1),]
  df.0 <- df_combined[which(treatment == 0),]

  m_ot <- lm(outcome ~ ps, data = df.1)
  ot <- suppressWarnings(predict(m_ot, newdata = as.data.frame(covariates))) # add newdata option to predict for full sample

  m_oc <- lm(outcome ~ ps, data = df.0)
  oc <- suppressWarnings(predict(m_oc, newdata = as.data.frame(covariates))) # add newdata option to predict for full sample

  te <- ot - oc
  ate <- mean(te)

  df_combined$te <- te

  att <- as.numeric(mean(df_combined[which(treatment == 1),'te']))
  atu <- as.numeric(mean(df_combined[which(treatment == 0),'te']))

  # ---- Calculate treatment weights and recalculate OLS
  N <- dim(df_combined[which(treatment == 1),])[1]
  var_ps <- as.numeric(var(df_combined[which(treatment == 1),'ps']))
  v1 <- var_ps*(N-1)/N

  N <- dim(df_combined[which(treatment == 0),])[1]
  var_ps <- as.numeric(var(df_combined[which(treatment == 0),'ps']))
  v0 <- var_ps*(N-1)/N

  p1 <- mean(treatment)
  p0 <- 1 - p1

  w1 <- as.numeric((p0*v0)/(p0*v0 + p1*v1))
  w0 <- as.numeric((p1*v1)/(p0*v0 + p1*v1))
  delta <- p1 - w1

  ols2 <- w1*att + w0*atu

  # ---- Save results in a list
  ret <- structure(list(OLS1 = beta_treated,
                        treatment_name = treatment_name,
                        'P(d=1)' = p1,
                        'P(d=0)' = p0,
                        w1 = w1,
                        w0 = w0,
                        delta = delta,
                        ATE = ate,
                        ATT = att,
                        ATU = atu,
                        OLS2 = ols2),


                   # ---- Define a new class
                   class = "hettreatreg")

  # ---- Display results
  if (verbose == TRUE) {
    print.hettreatreg(ret)
  }

  # ---- Return the list
  return(invisible(ret))
}
