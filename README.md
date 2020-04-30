# hettreatreg

## Overview

This package is used to compute diagnostics for linear regression when treatment effects are heterogeneous.

The package represents OLS estimates of the effect of a binary treatment as a weighted average of the average treatment effect on the treated (ATT) and the average treatment effect on the untreated (ATU). It follows Sloczynski (2019) ["Interpreting OLS Estimands When Treatment Effects Are Heterogeneous: Smaller Groups Get Larger Weights," available at http://people.brandeis.edu/~tslocz/Sloczynski_paper_regression.pdf] in estimating the OLS weights on these parameters, computing the associated model diagnostics, and reporting the implicit OLS estimate of the average treatment effect (ATE).

See Sloczynski (2019) for the underlying theoretical results and further details.

## Installing hettreatreg

This GitHub website hosts the source code. The package may be downloaded using:

```
library(devtools)
devtools::install_github("tslocz/hettreatreg")
```

## Authors

Tymon Sloczynski, Brandeis University, Waltham, MA. E-mail: tslocz [at] brandeis [dot] edu. Website: http://people.brandeis.edu/~tslocz/.

Mark McAvoy, Brandeis University, Waltham, MA. E-mail: mcavoy [at] brandeis [dot] edu.

Please do not hesitate to report bugs and share your comments on this program.
