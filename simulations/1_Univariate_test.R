


################################################################################
####
#### FPR = q = mean(RES)
####

q <- 0.3

RES <- NULL
scores <- NULL
n <- 10000

for(i in 1:1000)
{
  x <- rnorm(n)
  y <- x^2 + rnorm(n)
  #y <- x + rnorm(n)

  res <- screening_score(x,y, tol = 0)

  RES <- c(RES, res < 1/(qnorm(1-q/2)^2))
  scores <- c(scores, res)
}

mean(RES)


################################################################################
####
#### PLOT scores
####

library(ggplot2)

score_data <- data.frame(
  score = scores,
  selected = factor(
    RES,
    levels = c(FALSE, TRUE),
    labels = c("Not selected", "Selected")
  )
)

ggplot(score_data, aes(x = score, fill = selected)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    position = "identity",
    alpha = 0.45
  ) +
  scale_x_log10() +
  labs(
    title = "Distribution of screening scores",
    subtitle = "Logarithmic x-axis; vertical line is the selection threshold",
    x = expression(log[10]("screening score")),
    y = "Density",
    fill = "Decision"
  ) +
  geom_vline(
    xintercept = 1 / qnorm(1 - q / 2)^2,
    linetype = "dashed"
  ) +
  theme_minimal()


################################################################################
####
#### Gaussian distribution comparison
####

z_abs <- 1 / sqrt(scores)

z_data <- data.frame(z_abs = z_abs)

ggplot(z_data, aes(x = z_abs)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    boundary = 0
  ) +
  stat_function(
    fun = function(z) 2 * dnorm(z),
    linewidth = 1
  ) +
  coord_cartesian(xlim = c(0, quantile(z_abs, 0.99))) +
  labs(
    title = "Gaussian comparison for the transformed screening score",
    subtitle = expression(
      "Histogram of " * abs(Z) * " with half-normal reference density"
    ),
    x = expression(abs(Z) == 1 / sqrt("score")),
    y = "Density"
  ) +
  theme_minimal()


#################################################################################
####
#### linear link with different strength
####


q <- 0.3


RES <- NULL
u <- NULL
n <- 100

for(i in 1:1000)
{
  x <- rnorm(n)
  u <- c(runif(1), u)
  #y <- u[1]*x + rnorm(n)^2 + rnorm(n)
  y <- u[1]*x + rnorm(n)

  res <- screening_score(x,y, tol = 0)
  q <- 0.3

  RES <- c(res < 1/(qnorm(1-q/2)^2), RES)

}

### should be close to n
sum(RES)
### should be close to 1
mean(RES)


################################################################################
####
#### density of the selected variable / density of the non-selected variable
####

library(ggplot2)

plot_data <- data.frame(
  u = u,
  RES = factor(
    RES,
    levels = c(FALSE, TRUE),
    labels = c("FALSE: not selected", "TRUE: selected")
  )
)

ggplot(plot_data, aes(x = u, fill = RES)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 0.025,
    position = "identity",
    alpha = 0.45,
    boundary = 0
  ) +
  coord_cartesian(xlim = c(0, 0.5)) +
  labs(
    x = "Signal coefficient u",
    y = "Density",
    fill = "Screening decision",
    title = "Distribution of u by screening decision"
  ) +
  theme_minimal()

