<img src="sticker.svg" align="right" alt="" width="200" />

# remmm

`remmm` (resampling multivariate mixed models) provides resampling-based inference tools for multivariate mixed models and clustered data. 
The `remmm` package implements score-based sign-flipping procedures for fixed-effect inference when observations are not independent, including settings with multivariate responses.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("angeella/remmm")
library(remmm)
```

## Example

### Gaussian dependent variable


```r
remotes::install_github("livioivil/jointest")
library(remmm)

db <- simulateData(
  Sigma = "equicorrelation",
  rho = 0.5,
  beta = c(0.5, 1.2),
  gamma = 0.8,
  J = 8,
  nJ = rep(25, 8)
)

str(db)

fit <- clip(
  Y ~ X + Z,
  data = db,
  cluster = "id"
)

summary(jointest::p.adjust(fit))
summary(jointest::combine_tests(fit))
summary(jointest::combine_tests(fit, by = "model"))
summary(jointest::combine_tests(fit, by = "coefficient"))
```

### Non-Gaussian dependent variable

```r
remotes::install_github("livioivil/jointest")
library(remmm)

N <- 20
n <- rpois(N, 20)
reff <- rep(rnorm(N), n)

D <- data.frame(
  X1 = rnorm(length(reff)),
  X2 = rep(rnorm(N), n),
  Grp = factor(rep(rep(LETTERS[1:3], length.out = N), n)),
  Subj = rep(1:N, n)
)

D$Y <- rbinom(
  n = nrow(D),
  prob = plogis(2 * D$X1 * (D$Grp == "B") + 2 * D$X2 + reff),
  size = 1
)

res <- flip2sss(
  Y ~ Grp * X1 + X2,
  data = D,
  cluster = D$Subj,
  family = "binomial"
)

summary(res)
summary(jointest::combine_tests(res))
summary(jointest::combine_contrasts(res))
```

