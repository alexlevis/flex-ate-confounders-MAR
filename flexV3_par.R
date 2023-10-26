library(EnvStats)

load("fold2.Rdata")
delta.y <- 0.02
y.grid <- seq(0, 1, by = delta.y)

### fit nuisance models

## A
eta.hat <- 0.5

## Y | A
## use MLE for two shape parameters, stratified by A
mu.pars.0 <- ebeta(Y[A == 0], method = "mle")$parameters
mu.pars.1 <- ebeta(Y[A == 1], method = "mle")$parameters

## S | A, Y
pi.mod <- glm(S ~ A * Y, family = "binomial",
              data = cbind.data.frame(S = S, A = A, Y = Y))

## L1 | A, Y, S = 1
lambda1.mod <- glm(L1 ~ A * Y,
                   family = "binomial",
                   data = cbind.data.frame(L1 = L1, A = A, Y = Y)[S == 1,])

## L2 | A, Y, S = 1, L1
lambda2.mod <- lm(L2 ~ A + Y + Y:L1,
                  data = cbind.data.frame(L1 = L1, L2 = L2, A = A, Y = Y)[S == 1,])


### estimate parameters

load("fold1.Rdata")

pi.hat.1 <- predict.glm(pi.mod, type = "response",
                        newdata = cbind.data.frame(A = A, Y = Y))

lambda1.pmf <- function(a, y, l1) {
  p <- predict.glm(lambda1.mod, type = "response",
                   newdata = cbind.data.frame(A = a, Y = y))
  ifelse(l1 == 1, p, 1 - p)
}
lambda1.pmf <- Vectorize(lambda1.pmf)
lambda1.hat <- predict.glm(lambda1.mod, type = "response",
                           newdata = cbind.data.frame(A = A, Y = Y))

mu.pdf <- function(a, y) {
  ifelse(a == 0, dbeta(y, mu.pars.0[1], mu.pars.0[2]),
         dbeta(y, mu.pars.1[1], mu.pars.1[2]))
}
mu.pdf <- Vectorize(mu.pdf)

mu.pred.0 <- mu.pdf(a = 0, y = y.grid)
mu.pred.1 <- mu.pdf(a = 1, y = y.grid)

lambda2.pdf <- function(a, y, l1, l2) {
  dnorm(l2, predict.lm(lambda2.mod, newdata = cbind.data.frame(A = a, Y = y, L1 = l1)), sigma(lambda2.mod))
}
lambda2.pdf <- Vectorize(lambda2.pdf)

lambda.obs.mat.0 <- sapply(1:n, function(i) {
  if (S[i] == 0) { rep(NA, length(y.grid)) }
  else { 
    ifelse(L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
      lambda2.pdf(a = 0, y = y.grid, l1 = L1[i], l2 = L2[i])
  }
})
lambda.obs.mat.1 <- sapply(1:n, function(i) {
  if (S[i] == 0) { rep(NA, length(y.grid)) }
  else { 
    ifelse(L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
      lambda2.pdf(a = 1, y = y.grid, l1 = L1[i], l2 = L2[i])
  }
})

# save(n, A, Y, S, L1, L2, delta.y, y.grid,
#      eta.hat, pi.mod, pi.hat.1,
#      lambda1.mod, lambda1.hat, lambda1.pmf,
#      lambda.obs.mat.0, lambda.obs.mat.1,
#      mu.pars.0, mu.pars.1, mu.pdf, mu.pred.0, mu.pred.1,
#      lambda2.mod, lambda2.pdf,
#      file = "flexV3_F2_par.Rdata")

## IWOR estimator

# load("flexV3_F1_par.Rdata")
# load("flexV3_F2_par.Rdata")

gamma.hat.1 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.1[,i] * mu.pred.1))
}, simplify = 0)
gamma.hat.0 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.0[,i] * mu.pred.0))
}, simplify = 0)

beta.hat.1 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(y.grid * lambda.obs.mat.1[,i] * mu.pred.1))
}, simplify = 0)
beta.hat.0 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(y.grid * lambda.obs.mat.0[,i] * mu.pred.0))
}, simplify = 0)

IOR.1 <- mean( ifelse(S == 1, ( beta.hat.1 / gamma.hat.1) / 
                        pi.hat.1, 0))

IOR.0 <- mean( ifelse(S == 1, ( beta.hat.0 / gamma.hat.0) / 
                        pi.hat.1, 0))
IOR.ATE2 <- IOR.1 - IOR.0

(IOR.ATE + IOR.ATE2)/2 # 0.2554671

## IF-based estimates
res.agg <- read.csv("res-int-par-F2.csv")

IF.0 <- res.agg$b.01 + ifelse(A == 0, res.agg$b.02 / (1 - eta.hat), 0) +
  ifelse(S == 1, (beta.hat.0 / gamma.hat.0 - res.agg$b.01 +
                    ifelse(A == 0, ((gamma.hat.0 * (1 - eta.hat) + gamma.hat.1 * eta.hat) * 
                                      (Y - beta.hat.0 / gamma.hat.0) / gamma.hat.0 - 
                                      res.agg$b.02) / (1 - eta.hat), 0)) / pi.hat.1, 0)
IF.1 <- res.agg$b.11 + ifelse(A == 1, res.agg$b.12 / eta.hat, 0) +
  ifelse(S == 1, (beta.hat.1 / gamma.hat.1 - res.agg$b.11 +
                    ifelse(A == 1, ((gamma.hat.0 * (1 - eta.hat) + gamma.hat.1 * eta.hat) * 
                                      (Y - beta.hat.1 / gamma.hat.1) / gamma.hat.1 - 
                                      res.agg$b.12) / eta.hat, 0)) / pi.hat.1, 0)

IF.0.est <- mean(IF.0)
IF.1.est <- mean(IF.1)
IF.ATE2 <- IF.1.est - IF.0.est
IF.ATE.var2 <- var(IF.1 - IF.0)/n
IF.ATE2 + c(-1,1) * qnorm(0.975) * sqrt(IF.ATE.var2)

(IF.ATE + IF.ATE2)/2 # 0.2512982
(IF.ATE + IF.ATE2)/2 + c(-1,1) * qnorm(0.975) * sqrt((IF.ATE.var + IF.ATE.var2)/2)
