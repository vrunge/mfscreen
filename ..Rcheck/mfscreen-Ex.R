pkgname <- "mfscreen"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mfscreen')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("screening_score")
### * screening_score

flush(stderr()); flush(stdout())

### Name: screening_score
### Title: Reciprocal model-free screening statistic for one predictor
### Aliases: screening_score

### ** Examples

set.seed(123)

x_signal <- rnorm(200)
x_null <- rnorm(200)
y <- 2 * x_signal + rnorm(200)

screening_score(x_signal, y)
screening_score(x_null, y)




cleanEx()
nameEx("screening_score_R")
### * screening_score_R

flush(stderr()); flush(stdout())

### Name: screening_score_R
### Title: Reciprocal model-free screening statistic for one predictor
### Aliases: screening_score_R

### ** Examples


n <- 2500
q <- 0.05

# One predictor with a marginal association with y
x_signal <- rnorm(n)

# One predictor independent of y
x_null <- rnorm(n)

# Response generated from x_signal only
y <- x_signal^2 + 0.5 * x_signal + rnorm(n, sd = 1)

# Reciprocal distance statistic
D2_signal <- screening_score_R(x_signal, y)
D2_null <- screening_score_R(x_null, y)

# screening decisions
gamma <- qnorm(1 - q / 2)

selected <- c(
  signal = D2_signal <= 1 / gamma^2,
  null = D2_null <= 1 / gamma^2
)

data.frame(
  variable = c("signal", "null"),
  D2 = c(D2_signal, D2_null),
  selected = selected
)




cleanEx()
nameEx("screening_score_old_R")
### * screening_score_old_R

flush(stderr()); flush(stdout())

### Name: screening_score_old_R
### Title: Model-free screening statistic for one predictor
### Aliases: screening_score_old_R

### ** Examples

set.seed(123)

x_signal <- rnorm(200)
x_null <- rnorm(200)
y <- 2 * x_signal + rnorm(200)

screening_score_old_R(x_signal, y)
screening_score_old_R(x_null, y)




cleanEx()
nameEx("screening_test_matrix")
### * screening_test_matrix

flush(stderr()); flush(stdout())

### Name: screening_test_matrix
### Title: Reciprocal model-free screening for a predictor matrix
### Aliases: screening_test_matrix

### ** Examples

set.seed(123)

n <- 300
p <- 500

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
y <- 2 * X[, 1] - 1.5 * X[, 2] + rnorm(n)

result <- screening_test_matrix(X, y, q = 0.10)

result$selected_indices
head(result$scores)
result$threshold




cleanEx()
nameEx("screening_test_matrix_R")
### * screening_test_matrix_R

flush(stderr()); flush(stdout())

### Name: screening_test_matrix_R
### Title: Reciprocal model-free screening for a predictor matrix
### Aliases: screening_test_matrix_R

### ** Examples

set.seed(123)

n <- 250
p <- 100

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))

y <- 2 * X[, 1] - 1.5 * X[, 2] + rnorm(n)

result <- screening_test_matrix_R(X, y, q = 0.10)

result$selected_variables
head(result$scores)
result$threshold




cleanEx()
nameEx("screening_test_matrix_old_R")
### * screening_test_matrix_old_R

flush(stderr()); flush(stdout())

### Name: screening_test_matrix_old_R
### Title: Model-free marginal screening for a predictor matrix
### Aliases: screening_test_matrix_old_R

### ** Examples

set.seed(123)

n <- 2000
p <- 100

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))

y <- 2 * X[, 1] - 1.5 * X[, 2] + rnorm(n)

result <- screening_test_matrix_old_R(X, y, q = 0.10)

result$selected_variables
head(result$scores)
result$threshold




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
