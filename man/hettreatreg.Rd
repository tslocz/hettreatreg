% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hettreatreg.R
\name{hettreatreg}
\alias{hettreatreg}
\title{OLS Weights on Heterogeneous Treatment Effects}
\usage{
hettreatreg(outcome, treatment, covariates, verbose = FALSE)
}
\arguments{
\item{outcome}{the outcome variable.}

\item{treatment}{the treatment variable. The variable must be binary and coded 0 for the untreated units and 1 for the treated units.}

\item{covariates}{the list of control variables. The list must not include the treatment variable.}

\item{verbose}{logical. If \code{TRUE} estimation output is displayed.}
}
\value{
\item{\code{OLS}}{ OLS estimate of the treatment effect}
 \item{\code{P(d=1)}}{ proportion of treated units}
 \item{\code{P(d=0)}}{ proportion of untreated units}
 \item{\code{w1}}{ OLS weight on ATT}
 \item{\code{w0}}{ OLS weight on ATU}
 \item{\code{delta}}{ diagnostic for interpreting OLS as ATE}
 \item{\code{ATE}}{ implicit OLS estimate of ATE}
 \item{\code{ATT}}{ implicit OLS estimate of ATT}
 \item{\code{ATU}}{  implicit OLS estimate of ATU}
}
\description{
Computes diagnostics for linear regression when treatment effects are heterogeneous.
}
\details{
\code{hettreatreg} represents ordinary least squares (OLS) estimates of the effect of a binary treatment as a weighted average of the average treatment effect on the treated (ATT) and the average treatment effect on the untreated (ATU). The program estimates the OLS weights on these parameters, computes the associated model diagnostics, and reports the implicit OLS estimate of the average treatment effect (ATE). See Sloczynski (2019) for the underlying theoretical results and further details.

The arguments \code{outcome} and \code{treatment} are used to designate an outcome variable and a treatment variable, respectively. The treatment variable must be binary and coded 0 for the untreated units and 1 for the treated units. \code{covariates} is a list of control variables that must not include the treatment variable.

\code{hettreatreg} displays a number of statistics. \code{OLS} is the estimated regression coefficient on the treatment variable. \code{P(d=1)} and \code{P(d=0)} are the sample proportions of treated and untreated units, respectively. \code{w1} and \code{w0} are the OLS weights on ATT and ATU, respectively. \code{delta} is a diagnostic for interpreting OLS as ATE. \code{ATE}, \code{ATT}, and \code{ATU} are the implicit OLS estimates of the corresponding parameters. See Sloczynski (2019) for further details.

If you use this program in your work, please cite Sloczynski (2019).
}
\examples{
# load package
library(hettreatreg)

# read in data
data("nswcps")

# save the outcome variable
outcome <- nswcps$re78

# save the treatment variable
treated <- nswcps$treated

# select control variables
our_vars <- c("age", "age2", "educ", "black", "hispanic", "married", "nodegree")
covariates <- subset(nswcps, select = our_vars)

# run function
results <- hettreatreg(outcome, treated, covariates)
print(results)

}
\references{
Sloczynski, Tymon (2019). "Interpreting OLS Estimands When Treatment Effects Are Heterogeneous: Smaller Groups Get Larger Weights." Available at \url{http://people.brandeis.edu/~tslocz/Sloczynski_paper_regression.pdf}.
}
\author{
Tymon Sloczynski, Brandeis University, \email{tslocz@brandeis.edu}, \url{http://people.brandeis.edu/~tslocz/}

Maintained by: Mark McAvoy, Brandeis University, \email{mcavoy@brandeis.edu}

Please feel free to report bugs and share your comments on this program.
}
