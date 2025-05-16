#Packages
if (!require("logspline")) install.packages("logspline")
if (!require("purrr")) install.packages("purrr")
library(logspline)
library(purrr)

# Logspline Density Estimator 
logspline_density <- function(x, eval_points) {
  fit <- logspline(x)
  return(dlogspline(eval_points, fit))
}

# SAST with logspline density
SAST_logspline <- function(z, alpha = 0.05, conservative = TRUE, init = 200, h = 50) {
  T <- length(z)
  decisions <- rep(0, T)
  rejs <- 0
  den_null <- numeric(T)

  for (t in 1:T) {
    if (t <= init) {
      decisions[t] <- 0
      next
    }
    null_data <- z[1:(t-1)]
    null_density <- logspline_density(null_data, z[t])
    if (conservative) {
      threshold <- quantile(null_data, probs = 1 - alpha)
      decisions[t] <- as.integer(z[t] > threshold)
    } else {
      decisions[t] <- as.integer(null_density < alpha)
    }
  }
  return(list(decision = decisions))
}

# Data Generator
simulate_data_block <- function(T = 5000, burnin = 200, mu = 2.5, pi1 = 0.01) {
  pi_t <- rep(0.01, T)
  block_indices <- c(1000:1200, 2000:2200, 3000:3200, 4000:4200)
  pi_t[block_indices] <- pi1
  true_labels <- rbinom(T, 1, pi_t)
  z <- rnorm(T, mean = true_labels * mu, sd = 1)
  z_burnin <- rnorm(burnin, mean = 0, sd = 1)
  return(list(z_burnin = z_burnin, z = z, labels = true_labels))
}

# Evaluate Performance
evaluate_performance <- function(decision, true_labels) {
  TP <- sum(decision == 1 & true_labels == 1)
  FP <- sum(decision == 1 & true_labels == 0)
  FN <- sum(decision == 0 & true_labels == 1)
  FDP <- ifelse((TP + FP) == 0, 0, FP / (TP + FP))
  Power <- ifelse(sum(true_labels) == 0, 0, TP / sum(true_labels))
  MDR <- ifelse(sum(true_labels) == 0, 0, FN / sum(true_labels))
  return(c(FDP = FDP, Power = Power, MDR = MDR))
}

# Run Test
run_logspline_test <- function(reps = 100, m = 5000, pi1 = 0.01, mu = 2.5, init = 200) {
  fdp_final <- numeric(reps)
  for (r in 1:reps) {
    sim <- simulate_data_block(T = m, burnin = init, mu = mu, pi1 = pi1)
    z_all <- c(sim$z_burnin, sim$z)
    res <- SAST_logspline(z_all, alpha = 0.05, conservative = TRUE, init = init)
    decision <- res$decision[(init+1):(init+m)]
    perf <- evaluate_performance(decision, sim$labels)
    fdp_final[r] <- perf["FDP"]
  }
  mean_fdr <- mean(fdp_final)
  png("FDR_logspline_pi1_0.01.png", width = 800, height = 600)
  plot(1:reps, fdp_final, type = "l", col = "blue", ylab = "FDR", xlab = "Repetition",
       main = "FDR Across Runs Using Logspline Density (pi1 = 0.01)")
  abline(h = 0.05, lty = 2)
  dev.off()
  cat("✅ Mean FDR using logspline =", mean_fdr, "\n")
  return(mean_fdr)
}


run_logspline_test()
