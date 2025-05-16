# 可选参数进行调参："h_sweep", "mu_sweep", "burnin_sweep", "m_sweep", "pi1_sweep", "conserv_sweep"
mode <- "pi1_sweep"  # 修改此处进行调参
pattern <- "block"    # 修改此处选择数据结构："block", "constant", "linear", "sine"

# 加载包
if (!require("kedd")) install.packages("kedd")
if (!require("purrr")) install.packages("purrr")
library(kedd)
library(purrr)

# 加载 SAST
source("/Users/xixinzhao/FDR/SAST-main/R/SAST.R")  # SAST函数的路径

# 通用数据生成器（支持四种结构："block", "constant", "linear", "sine"）
simulate_data <- function(pattern = "block", T = 5000, burnin = 500, mu = 2.5, pi1 = 0.2) {
  if (pattern == "block") {
    pi_t <- rep(0.01, T)
    block_indices <- c(1000:1200, 2000:2200, 3000:3200, 4000:4200)
    pi_t[block_indices] <- pi1
  } else if (pattern == "constant") {
    pi_t <- rep(pi1, T)
  } else if (pattern == "linear") {
    pi_t <- seq(0, pi1, length.out = T)
  } else if (pattern == "sine") {
    pi_t <- 0.005 + pi1 * sin(2 * pi * (1:T) / 1000)
    pi_t <- pmin(pmax(pi_t, 0), 1)
  } else {
    stop("Invalid pattern. Choose from 'block', 'constant', 'linear', 'sine'.")
  }
  labels <- rbinom(T, 1, pi_t)
  z <- rnorm(T, mean = labels * mu, sd = 1)
  z_burnin <- rnorm(burnin)
  return(list(z_burnin = z_burnin, z = z, labels = labels))
}

#  SAST运行
run_SAST <- function(z_burnin, z, init = 200, conservative = TRUE, h = 50) {
  z_total <- c(z_burnin, z)
  res <- SAST(z_total, alpha = 0.05, conservative = conservative, init = init, h = h)
  return(res$decision[(length(z_burnin)+1):(length(z_burnin)+length(z))])
}
#  SAST评估
evaluate_performance <- function(decision, true_labels) {
  TP <- sum(decision == 1 & true_labels == 1)
  FP <- sum(decision == 1 & true_labels == 0)
  FN <- sum(decision == 0 & true_labels == 1)
  FDP <- ifelse((TP + FP) == 0, 0, FP / (TP + FP))
  Power <- ifelse(sum(true_labels) == 0, 0, TP / sum(true_labels))
  MDR <- ifelse(sum(true_labels) == 0, 0, FN / sum(true_labels))
  c(FDP = FDP, Power = Power, MDR = MDR)
}

#参数调整 1：带宽h sweep 
run_h_sweep <- function() {
  h_values <- c(10, 30, 50, 100, 200) #调整不同的带宽h值
  reps <- 100
  fdrs <- numeric()
  for (h in h_values) {
    fdp_vec <- replicate(reps, {
      sim <- simulate_data()
      dec <- run_SAST(sim$z_burnin, sim$z, h = h)
      evaluate_performance(dec, sim$labels)["FDP"]
    })
    fdrs <- c(fdrs, mean(fdp_vec))
  }
  plot(h_values, fdrs, type = "b", pch = 19, col = "darkred", lwd = 2,
       ylim = c(0, 0.1), xlab = "Bandwidth h", ylab = "FDR",
       main = "FDR vs Bandwidth h")
  abline(h = 0.05, lty = 2)
}

# 参数调整 2：信号强度 mu sweep
run_mu_sweep <- function() {
  mu_values <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5) #调整不同的信号强度mu值
  reps <- 100
  results <- data.frame(mu = mu_values, FDR = NA, Power = NA, MDR = NA)
  
  for (i in seq_along(mu_values)) {
    metrics <- replicate(reps, {
      sim <- simulate_data(mu = mu_values[i])
      dec <- run_SAST(sim$z_burnin, sim$z)
      evaluate_performance(dec, sim$labels)
    })
    results[i, 2:4] <- rowMeans(metrics)
  }
  
  plot(results$mu, results$FDR, type = "b", col = "blue", lwd = 2, pch = 19,
       ylim = c(0, 0.1), xlab = "mu", ylab = "FDR", main = "FDR vs Signal Strength mu")
  abline(h = 0.05, lty = 2)
  print(results)
}

# 参数调整 3：burn-in sweep
run_burnin_sweep <- function() {
  burnin_values <- c(50, 100, 200, 400, 800) #调整不同的burn-in值
  reps <- 100
  fdrs <- numeric()
  for (b in burnin_values) {
    fdp_vec <- replicate(reps, {
      sim <- simulate_data(burnin = b)
      dec <- run_SAST(sim$z_burnin, sim$z, init = b)
      evaluate_performance(dec, sim$labels)["FDP"]
    })
    fdrs <- c(fdrs, mean(fdp_vec))
  }
  plot(burnin_values, fdrs, type = "b", col = "blue", lwd = 2, pch = 19,
       ylim = c(0, 0.1), xlab = "Burn-in (init)", ylab = "FDR",
       main = "FDR vs Burn-in")
  abline(h = 0.05, lty = 2)
}

# 参数调整 4：样本数m sweep
run_m_sweep <- function() {
  m_values <- c(300, 1000, 3000, 5000) #调整不同的样本数m值
  reps <- 100
  colors <- c("red", "blue", "green", "black")
  
  plot(NULL, xlim = c(200, 5000), ylim = c(0, 0.1), xlab = "Time", ylab = "FDR",
       main = "FDR vs Time (by sample size m)")
  abline(h = 0.05, lty = 2)
  
  for (i in seq_along(m_values)) {
    checkpoints <- seq(100, m_values[i], length.out = 5)
    fdp_mat <- matrix(0, reps, length(checkpoints))
    for (r in 1:reps) {
      sim <- simulate_data(T = m_values[i])
      dec <- run_SAST(sim$z_burnin, sim$z)
      for (j in seq_along(checkpoints)) {
        perf <- evaluate_performance(dec[1:checkpoints[j]], sim$labels[1:checkpoints[j]])
        fdp_mat[r, j] <- perf["FDP"]
      }
    }
    lines(checkpoints, colMeans(fdp_mat), col = colors[i], lwd = 2)
  }
  legend("topright", legend = paste0("m = ", m_values), col = colors, lwd = 2)
}

# 参数调整 5：pi1 sweep 
run_pi1_sweep <- function() {
  pi1_values <- c(0.01, 0.02, 0.05, 0.1, 0.2) #调整不同的pi1值
  reps <- 100
  fdrs <- numeric()
  for (pi1 in pi1_values) {
    fdp_vec <- replicate(reps, {
      sim <- simulate_data(pi1 = pi1)
      dec <- run_SAST(sim$z_burnin, sim$z)
      evaluate_performance(dec, sim$labels)["FDP"]
    })
    fdrs <- c(fdrs, mean(fdp_vec))
  }
  plot(pi1_values, fdrs, type = "b", pch = 19, col = "darkgreen", lwd = 2,
       ylim = c(0, 0.1), xlab = "pi1", ylab = "FDR",
       main = "FDR vs pi1")
  abline(h = 0.05, lty = 2)
}

# 参数调整 6：conservative vs non-conservative sweep
run_conserv_sweep <- function() {
  pi1_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4) #调整不同的pi1值，对比不同的pi1值的conservative vs non-conservative
  reps <- 100
  init <- 200
  results <- data.frame(pi1 = rep(pi1_vals, each = 2),
                        conservative = rep(c(TRUE, FALSE), times = length(pi1_vals)),
                        FDR = NA)
  for (i in seq_along(pi1_vals)) {
    for (cons in c(TRUE, FALSE)) {
      fdrs <- numeric(reps)
      for (r in 1:reps) {
        sim <- simulate_data(pattern = pattern, pi1 = pi1_vals[i], burnin = init)
        decision <- run_SAST(sim$z_burnin, sim$z, init = init, conservative = cons)
        perf <- evaluate_performance(decision, sim$labels)
        fdrs[r] <- perf["FDP"]
      }
      results$FDR[(i - 1) * 2 + ifelse(cons, 1, 2)] <- mean(fdrs)
    }
  }
  plot(results$pi1[results$conservative == TRUE],
       results$FDR[results$conservative == TRUE],
       type = "b", pch = 19, col = "blue", ylim = c(0, 0.1),
       xlab = "pi1", ylab = "FDR",
       main = paste("Conservative vs Non-Conservative (", pattern, ")"))
  lines(results$pi1[results$conservative == FALSE],
        results$FDR[results$conservative == FALSE],
        type = "b", pch = 17, col = "red")
  abline(h = 0.05, lty = 2)
  legend("topright", legend = c("Conservative", "Non-Conservative"),
         col = c("blue", "red"), pch = c(19, 17), lwd = 2)
  print(results)
}

# 控制入口
if (mode == "h_sweep") {
  run_h_sweep()
} else if (mode == "mu_sweep") {
  run_mu_sweep()
} else if (mode == "burnin_sweep") {
  run_burnin_sweep()
} else if (mode == "m_sweep") {
  run_m_sweep()
} else if (mode == "pi1_sweep") {
  run_pi1_sweep()
} else if (mode == "conserv_sweep") {
  run_conserv_sweep()
} else {
  stop("❌ Unknown mode")
}
