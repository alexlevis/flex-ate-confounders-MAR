load("flexV3_F1_par.Rdata")

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
  # L1 = 1, A = 1
  lambda.hyp.hat <-
    lambda1.pmf(a = 1, y = y.grid, l1 = 1) * 
    lambda2.pdf(a = 1, y = y.grid, l1 = 1, l2 = lambda.grid[j])
  gamma.hyp.1.1[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.1)
  beta.hyp.1.1[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.1)
  
  # L1 = 1, A = 0
  lambda.hyp.hat <-
    lambda1.pmf(a = 0, y = y.grid, l1 = 1) * 
    lambda2.pdf(a = 0, y = y.grid, l1 = 1, l2 = lambda.grid[j])
  gamma.hyp.1.0[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.0)
  beta.hyp.1.0[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.0)
  
  # L1 = 0, A = 1
  lambda.hyp.hat <-
    lambda1.pmf(a = 1, y = y.grid, l1 = 0) * 
    lambda2.pdf(a = 1, y = y.grid, l1 = 0, l2 = lambda.grid[j])
  gamma.hyp.0.1[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.1)
  beta.hyp.0.1[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.1)
  
  # L1 = 0, A = 0
  lambda.hyp.hat <-
    lambda1.pmf(a = 0, y = y.grid, l1 = 0) * 
    lambda2.pdf(a = 0, y = y.grid, l1 = 0, l2 = lambda.grid[j])
  gamma.hyp.0.0[j] <- delta.y * sum(lambda.hyp.hat * mu.pred.0)
  beta.hyp.0.0[j] <- delta.y * sum(y.grid * lambda.hyp.hat * mu.pred.0)
}

b.12 <- b.02 <- 0
lambda.obs.hat.0 <- 
  lambda2.pdf(a = A[i], y = Y[i], l1 = 0, l2 = lambda.grid)

lambda.obs.hat.1 <- 
  lambda2.pdf(a = A[i], y = Y[i], l1 = 1, l2 = lambda.grid)

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

