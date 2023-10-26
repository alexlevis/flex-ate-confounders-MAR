library(haldensify)
library(data.table)
library(ggplot2)
library(SuperLearner)
library(statmod)

set.seed(158)

## training fold
n <- 2500
A <- rbinom(n, 1, 0.5)
Y <- rbeta(n, shape1 = ifelse(A == 0, 2, 4), shape2 = ifelse(A == 0, 4, 2)) # mean 1/3 when A = 0, 2/3 when A = 1

S <- rbinom(n, 1, plogis(-0.35 + 0.5 * A + 0.18 * Y + 0.05 * A * Y))
L1 <- rbinom(n, 1, plogis(-0.6 + 0.5 * A + 0.25 * Y + 0.1 * A * Y))
L2 <- rnorm(n, 1.0 * A + Y + 2.5 * L1 * Y, 1.25)

# HAL-based density estimate of Y|A
hal_mu <- haldensify(
  A = Y, W = A,
  n_bins = 25, grid_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 100)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 3,
  reduce_basis = 1 / sqrt(n)
)
hal_mu

delta.y <- 0.02
new_y <- seq(0, 1, by = delta.y)
new_dat <- as.data.table(list(
  y = new_y,
  a_zero = rep(0, length(new_y)),
  a_one = rep(1, length(new_y))
))
new_dat[, pred_a_zero := predict(hal_mu, trim = F,
                                 new_A = new_dat$y,
                                 new_W = new_dat$a_zero)]
new_dat[, pred_a_one := predict(hal_mu, trim = F,
                                new_A = new_dat$y,
                                new_W = new_dat$a_one)]

# visualize results
dens_dat <-  melt(new_dat, id = c("y"),
                  measure.vars = c("pred_a_zero", "pred_a_one"))
p_dens <- ggplot(dens_dat, aes(x = y, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  stat_function(fun = dbeta, args = list(shape1 = 2, shape2 = 4),
                colour = "#F8766D", linetype = "dashed") +
  stat_function(fun = dbeta, args = list(shape1 = 4, shape2 = 2),
                colour = "#00BFC4", linetype = "dashed") +
  labs(
    x = "",
    y = "Estimated conditional density",
    title = "Conditional density estimates p(Y|A)"
  ) +
  theme_bw() +
  theme(legend.position = "none")
p_dens
# ggplot_build(p_dens)$data


# HAL-based density estimate of L_p2 | L_p1, A, Y
hal_lambda <- haldensify(
  A = L2[S == 1], W = cbind(L1, A, Y)[S == 1, ],
  n_bins = 35, grid_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 50)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 3,
  reduce_basis = 1 / sqrt(n)
)
hal_lambda

eta.hat <- 0.5
# pi.mod <- glm(S ~ A * Y, family = "binomial", 
#               data = cbind.data.frame(S = S, A = A, Y = Y))
# lambda1.mod <- glm(L1 ~ A * Y,
#                    family = "binomial", 
#                    data = cbind.data.frame(L1 = L1, A = A, Y = Y)[S == 1,])
pi.mod <- SuperLearner(S, X = cbind.data.frame(A, Y),
                       family = "binomial",
                       SL.library = c("SL.glm.interaction",
                                      "SL.polymars"))
lambda1.mod <- SuperLearner(Y = L1[S == 1], X = cbind.data.frame(A, Y)[S == 1,],
                            family = "binomial",
                            SL.library = c("SL.glm.interaction",
                                           "SL.polymars"))

## testing fold
A <- rbinom(n, 1, 0.5)
Y <- rbeta(n, shape1 = ifelse(A == 0, 2, 4), shape2 = ifelse(A == 0, 4, 2)) # mean 1/3 when A = 0, 2/3 when A = 1

S <- rbinom(n, 1, plogis(-0.35 + 0.5 * A + 0.18 * Y + 0.05 * A * Y))
L1 <- rbinom(n, 1, plogis(-0.6 + 0.5 * A + 0.25 * Y + 0.1 * A * Y))
L2 <- rnorm(n, 1.0 * A + Y + 2.5 * L1 * Y, 1.25)

pi.hat.1 <- predict.SuperLearner(pi.mod, type = "response",
                                 newdata = cbind.data.frame(A = A, Y = Y))$pred
# pi.hat.1 <- predict(pi.mod, type = "response",
#                     newdata = cbind.data.frame(A = A, Y = Y))
lambda1.hat <- predict.SuperLearner(lambda1.mod, type = "response",
                                    newdata = cbind.data.frame(A = A, Y = Y))$pred
# lambda1.hat <- predict(lambda1.mod, type = "response", 
#                        newdata = cbind.data.frame(A = A, Y = Y))

lambda.obs.mat.1 <- sapply(1:n, function(i) {
  if (S[i] == 0) { rep(NA, length(new_y)) }
  else { 
    ifelse(L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
      predict(hal_lambda, trim = F,
              new_A = rep(L2[i], length(new_y)),
              new_W = cbind(L1 = rep(L1[i], length(new_y)), 
                                 A = rep(1, length(new_y)),
                                 Y = new_y))
  }
})
lambda.obs.mat.0 <- sapply(1:n, function(i) {
  if (S[i] == 0) { rep(NA, length(new_y)) }
  else { 
    ifelse(L1[i] == 1, lambda1.hat[i], 1 - lambda1.hat[i]) *
      predict(hal_lambda, trim = F,
              new_A = rep(L2[i], length(new_y)),
              new_W = cbind(L1 = rep(L1[i], length(new_y)), 
                            A = rep(0, length(new_y)),
                            Y = new_y))
  }
})

# save(n, A, Y, S, L1, L2, new_dat, delta.y, new_y,
#      eta.hat, pi.mod, pi.hat.1,
#      lambda1.mod, lambda1.hat,
#      lambda.obs.mat.0, lambda.obs.mat.1,
#      hal_mu, hal_lambda,
#      file = "flexV3_F2.Rdata")
# load("flexV3_F1.Rdata")
# load("flexV3_F2.Rdata")


## IWOR estimator

gamma.hat.1 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.1[,i] * new_dat$pred_a_one))
}, simplify = 0)
gamma.hat.0 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(lambda.obs.mat.0[,i] * new_dat$pred_a_zero))
}, simplify = 0)

beta.hat.1 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(new_y * lambda.obs.mat.1[,i] * new_dat$pred_a_one))
}, simplify = 0)
beta.hat.0 <- sapply(1:n, function(i) {
  ifelse(S[i] == 0, 0,
         delta.y * sum(new_y * lambda.obs.mat.0[,i] * new_dat$pred_a_zero))
}, simplify = 0)

IOR.1 <- mean( ifelse(S == 1, ( beta.hat.1 / gamma.hat.1) / 
                        pi.hat.1, 0))

IOR.0 <- mean( ifelse(S == 1, ( beta.hat.0 / gamma.hat.0) / 
                        pi.hat.1, 0))
IOR.ATE <- IOR.1 - IOR.0

(IOR.ATE + IOR.ATE2)/2


## IF-based estimates
res.agg <- read.csv("res-int-F2.csv")

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

(IF.ATE + IF.ATE2)/2
(IF.ATE + IF.ATE2)/2 + c(-1,1) * qnorm(0.975) * sqrt((IF.ATE.var + IF.ATE.var2)/2)

