library(haldensify)
library(SuperLearner)

load("flexV3_F2.Rdata")

args = commandArgs(TRUE)
i <- as.integer( Sys.getenv('SLURM_ARRAY_TASK_ID') )

delta.lambda <- 0.0555
lambda.grid <- seq(-3.3, 7.8, by = delta.lambda)

gamma.hyp.1.0 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.1.0 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.0.0 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.0.0 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.1.1 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.1.1 <- rep(NA, nrow = length(lambda.grid))
gamma.hyp.0.1 <- rep(NA, nrow = length(lambda.grid))
beta.hyp.0.1 <- rep(NA, nrow = length(lambda.grid))
for (j in 1:length(lambda.grid)) {
  lambda1.hat.1 <- predict.SuperLearner(lambda1.mod, type = "response",
                                        newdata = cbind.data.frame(A = 1, Y = new_y))$pred
  lambda1.hat.0 <- predict.SuperLearner(lambda1.mod, type = "response",
                                        newdata = cbind.data.frame(A = 0, Y = new_y))$pred
  
  # L1 = 1, A = 1
  lambda.hyp.hat <-
    lambda1.hat.1 * 
    predict(hal_lambda, trim = F,
            new_A = rep(lambda.grid[j], length(new_y)),
            new_W = cbind(L1 = rep(1, length(new_y)), 
                          A = rep(1, length(new_y)),
                          Y = new_y))
  gamma.hyp.1.1[j] <- delta.y * sum(lambda.hyp.hat * new_dat$pred_a_one)
  beta.hyp.1.1[j] <- delta.y * sum(new_y * lambda.hyp.hat * new_dat$pred_a_one)
  
  # L1 = 1, A = 0
  lambda.hyp.hat <-
    lambda1.hat.0 *
    predict(hal_lambda, trim = F,
            new_A = rep(lambda.grid[j], length(new_y)),
            new_W = cbind(L1 = rep(1, length(new_y)), 
                          A = rep(0, length(new_y)),
                          Y = new_y))
  gamma.hyp.1.0[j] <- delta.y * sum(lambda.hyp.hat * new_dat$pred_a_zero)
  beta.hyp.1.0[j] <- delta.y * sum(new_y * lambda.hyp.hat * new_dat$pred_a_zero)
  
  # L1 = 0, A = 1
  lambda.hyp.hat <-
    (1 - lambda1.hat.1) *
    predict(hal_lambda, trim = F,
            new_A = rep(lambda.grid[j], length(new_y)),
            new_W = cbind(L1 = rep(0, length(new_y)), 
                          A = rep(1, length(new_y)),
                          Y = new_y))
  gamma.hyp.0.1[j] <- delta.y * sum(lambda.hyp.hat * new_dat$pred_a_one)
  beta.hyp.0.1[j] <- delta.y * sum(new_y * lambda.hyp.hat * new_dat$pred_a_one)
  
  # L1 = 0, A = 0
  lambda.hyp.hat <-
    (1 - lambda1.hat.0) *
    predict(hal_lambda, trim = F,
            new_A = rep(lambda.grid[j], length(new_y)),
            new_W = cbind(L1 = rep(0, length(new_y)), 
                          A = rep(0, length(new_y)),
                          Y = new_y))
  gamma.hyp.0.0[j] <- delta.y * sum(lambda.hyp.hat * new_dat$pred_a_zero)
  beta.hyp.0.0[j] <- delta.y * sum(new_y * lambda.hyp.hat * new_dat$pred_a_zero)
}

b.12 <- b.02 <- 0
lambda.obs.hat.0 <- 
  predict(hal_lambda, trim = F,
          new_A = lambda.grid,
          new_W = cbind(L1 = rep(0, length(lambda.grid)), 
                        A = rep(A[i], length(lambda.grid)), 
                        Y = rep(Y[i], length(lambda.grid))))
# lambda.obs.hat.0 <- lambda.obs.hat.0 / (delta.lambda * sum(lambda.obs.hat.0))

lambda.obs.hat.1 <- 
  predict(hal_lambda, trim = F,
          new_A = lambda.grid,
          new_W = cbind(L1 = rep(1, length(lambda.grid)), 
                        A = rep(A[i], length(lambda.grid)), 
                        Y = rep(Y[i], length(lambda.grid))))
# lambda.obs.hat.1 <- lambda.obs.hat.1 / (delta.lambda * sum(lambda.obs.hat.1))

b.01 <- delta.lambda * 
  sum((1 - lambda1.hat[i]) * lambda.obs.hat.0 * beta.hyp.0.0 / 
        gamma.hyp.0.0 +
        lambda1.hat[i] * lambda.obs.hat.1 * beta.hyp.1.0 / 
        gamma.hyp.1.0)
b.11 <- delta.lambda * 
  sum((1 - lambda1.hat[i]) * lambda.obs.hat.0 * beta.hyp.0.1 / 
        gamma.hyp.0.1 +
        lambda1.hat[i] * lambda.obs.hat.1 * beta.hyp.1.1 / 
        gamma.hyp.1.1)

if (A[i] == 0) {
  b.02 <- delta.lambda * 
    sum((1 - lambda1.hat[i]) * lambda.obs.hat.0 * 
          ((1 - eta.hat) * (gamma.hyp.0.0) + 
             eta.hat * (gamma.hyp.0.1)) *
          (Y[i] - beta.hyp.0.0 / gamma.hyp.0.0) /
          gamma.hyp.0.0 +
          lambda1.hat[i] * lambda.obs.hat.1 *
          ((1 - eta.hat) * (gamma.hyp.1.0) + 
             eta.hat * (gamma.hyp.1.1)) *
          (Y[i] - beta.hyp.1.0 / gamma.hyp.1.0) /
          gamma.hyp.1.0)
} else {
  b.12 <- delta.lambda * 
    sum((1 - lambda1.hat[i]) * lambda.obs.hat.0 * 
          ((1 - eta.hat) * (gamma.hyp.0.0) + 
             eta.hat * (gamma.hyp.0.1)) *
          (Y[i] - beta.hyp.0.1 / gamma.hyp.0.1) /
          gamma.hyp.0.1 +
          lambda1.hat[i] * lambda.obs.hat.1 *
          ((1 - eta.hat) * (gamma.hyp.1.0) + 
             eta.hat * (gamma.hyp.1.1)) *
          (Y[i] - beta.hyp.1.1 / gamma.hyp.1.1) /
          gamma.hyp.1.1)
}



export <- data.frame(id = i, b.01 = b.01, b.11 = b.11,
                     b.02 = b.02, b.12 = b.12)

write.csv(export, file = paste0("res-int-", i, ".csv"), row.names = F)

