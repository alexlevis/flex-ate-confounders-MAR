library(SuperLearner)

set.seed(776)
n <- 50000
A <- rbinom(n, 1, 0.5)
Y <- rbeta(n, shape1 = ifelse(A == 0, 2, 4), shape2 = ifelse(A == 0, 4, 2)) 

L1 <- rbinom(n, 1, plogis(-0.6 + 0.5 * A + 0.25 * Y + 0.1 * A * Y))
L2 <- rnorm(n, 1.0 * A + Y + 2.5 * L1 * Y, 1.25)

mean.mod <- SuperLearner(Y, X = cbind.data.frame(L1, L2, A),
                       SL.library = c("SL.ranger", "SL.polymars", 
                                      "SL.glm", "SL.nnet"))
prop.mod <- SuperLearner(A, X = cbind.data.frame(L1, L2),
                         family = binomial(),
                         SL.library = c("SL.ranger", "SL.polymars", 
                                        "SL.glm", "SL.nnet"))

A <- rbinom(n, 1, 0.5)
Y <- rbeta(n, shape1 = ifelse(A == 0, 2, 4), shape2 = ifelse(A == 0, 4, 2)) 

L1 <- rbinom(n, 1, plogis(-0.6 + 0.5 * A + 0.25 * Y + 0.1 * A * Y))
L2 <- rnorm(n, 1.0 * A + Y + 2.5 * L1 * Y, 1.25)

dat.0 <- cbind.data.frame(L1 = L1, L2 = L2, A = 0)
dat.1 <- cbind.data.frame(L1 = L1, L2 = L2, A = 1)

mean.hat.0 <- predict(mean.mod, newdata = dat.0)$pred
mean.hat.1 <- predict(mean.mod, newdata = dat.1)$pred
prop.hat <- predict.SuperLearner(prop.mod, type = "response",
                                 newdata = cbind.data.frame(L1 = L1, L2 = L2))$pred

IF.0 <- mean.hat.0 + 1 * (A == 0) * (Y - mean.hat.0) / (1 - prop.hat)
IF.1 <- mean.hat.1 + 1 * (A == 1) * (Y - mean.hat.1) / prop.hat

ATE.0.true <- mean(IF.0)
ATE.1.true <- mean(IF.1)

# true ATE est
ATE.true <- mean(mean.hat.1 + 1 * (A == 1) * (Y - mean.hat.1) / prop.hat) -
  mean(mean.hat.0 + 1 * (A == 0) * (Y - mean.hat.0) / (1 - prop.hat))

ATE.true + c(-1,1) * qnorm(0.975) * sqrt(mean((IF.1 - ATE.1.true - IF.0 + ATE.0.true)^2) / n)

# rm(n, A, Y, L1, L2, mean.mod,prop.mod, dat.0, dat.1,
#    prop.hat, mean.hat.0, mean.hat.1, IF.0, IF.1)


## parametric truth computation
delta.y <- 0.0025
y.grid <- seq(0, 1, by = delta.y)
delta.l <- 0.013875
l.grid <- seq(-3.3, 7.8, by = delta.l)

## true nuisance models
eta <- 0.5
y.pdf <- function(a, y) {
  ifelse(a == 0, dbeta(y, 2, 4), dbeta(y, 4, 2))
}
y.pdf <- Vectorize(y.pdf)
l1.pmf <- function(a, y, l1) {
  p <- plogis(-0.6 + 0.5 * a + 0.25 * y + 0.1 * a * y)
  ifelse(l1 == 1, p, 1 - p)
}
l1.pmf <- Vectorize(l1.pmf)
l2.pdf <- function(a, y, l1, l2) {
  dnorm(l2, 1.0 * a + y + 2.5 * l1 * y, 1.25)
}
l2.pdf <- Vectorize(l2.pdf)

# L1 = 0, A = 0
gamma.0.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 0) *
        l2.pdf(a = 0, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid)) * delta.y
}, simplify = 0)
beta.0.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 0) *
        l2.pdf(a = 0, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.0.0 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 0, y = y.grid[j], l1 = 0, l2 = l.grid)
})

# L1 = 0, A = 1
gamma.0.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 0) *
        l2.pdf(a = 1, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid)) * delta.y
}, simplify = 0)
beta.0.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 0) *
        l2.pdf(a = 1, y = y.grid, l1 = 0, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.0.1 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 1, y = y.grid[j], l1 = 0, l2 = l.grid)
})

# L1 = 1, A = 0
gamma.1.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 1) *
        l2.pdf(a = 0, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid)) * delta.y
}, simplify = 0)
beta.1.0 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 0, y = y.grid, l1 = 1) *
        l2.pdf(a = 0, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 0, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.1.0 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 0, y = y.grid[j], l1 = 1, l2 = l.grid)
})


# L1 = 1, A = 1
gamma.1.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 1) *
        l2.pdf(a = 1, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid)) * delta.y
}, simplify = 0)
beta.1.1 <- sapply(1:length(l.grid), function(j) {
  sum(l1.pmf(a = 1, y = y.grid, l1 = 1) *
        l2.pdf(a = 1, y = y.grid, l1 = 1, l2 = l.grid[j]) *
        y.pdf(a = 1, y = y.grid) * y.grid) * delta.y
}, simplify = 0)
lambda.mat.1.1 <- sapply(1:length(y.grid), function(j) {
  l2.pdf(a = 1, y = y.grid[j], l1 = 1, l2 = l.grid)
})

ATE.true.0 <- 
  eta * sum(delta.y * y.pdf(a = 0, y = y.grid) *
              (l1.pmf(a = 0, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.0 / gamma.0.0 * 
                         lambda.mat.0.0[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 0, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.0 / gamma.1.0 * 
                         lambda.mat.1.0[,j] * delta.l)
                 }, simplify = 0))) +
  eta * sum(delta.y * y.pdf(a = 1, y = y.grid) *
              (l1.pmf(a = 1, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.0 / gamma.0.0 * 
                         lambda.mat.0.1[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 1, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.0 / gamma.1.0 * 
                         lambda.mat.1.1[,j] * delta.l)
                 }, simplify = 0)))

ATE.true.1 <- 
  eta * sum(delta.y * y.pdf(a = 0, y = y.grid) *
              (l1.pmf(a = 0, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.1 / gamma.0.1 * 
                         lambda.mat.0.0[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 0, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.1 / gamma.1.1 * 
                         lambda.mat.1.0[,j] * delta.l)
                 }, simplify = 0))) +
  eta * sum(delta.y * y.pdf(a = 1, y = y.grid) *
              (l1.pmf(a = 1, y = y.grid, l1 = 0) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.0.1 / gamma.0.1 * 
                         lambda.mat.0.1[,j] * delta.l)
                 }, simplify = 0) +
                 l1.pmf(a = 1, y = y.grid, l1 = 1) *
                 sapply(1:length(y.grid), function(j) {
                   sum(beta.1.1 / gamma.1.1 * 
                         lambda.mat.1.1[,j] * delta.l)
                 }, simplify = 0)))

ATE.true <- ATE.true.1 - ATE.true.0

