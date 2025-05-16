
mode <- "block"  # 修改这个值切换运行模式，选择运行模式：block / constant / linear / sine

# 加载必要包
if (!require("onlineFDR")) install.packages("onlineFDR")
if (!require("kedd")) install.packages("kedd")
if (!require("purrr")) install.packages("purrr")
library(onlineFDR)
library(kedd)
library(purrr)

# 加载 SAST 函数
source("/Users/xixinzhao/FDR/SAST-main/R/SAST.R")  # 路径

# 模式定义：block / constant / linear / sine
simulate_data <- switch(mode,
                        block = function(T = 5000, burnin = 500, mu = 2.5) {
                          pi_t <- rep(0.01, T)
                          block_indices <- c(1000:1200, 2000:2200, 3000:3200, 4000:4200)
                          pi_t[block_indices] <- ifelse(block_indices <= 3000, 0.2, 0.4)
                          true_labels <- rbinom(T, 1, pi_t)
                          z <- rnorm(T, mean = true_labels * mu, sd = 1)
                          z_burnin <- rnorm(burnin)
                          list(z_burnin = z_burnin, z = z, labels = true_labels)
                        },
                        constant = function(T = 5000, burnin = 500, mu = 2.5) {
                          pi_t <- rep(0.1, T)
                          true_labels <- rbinom(T, 1, pi_t)
                          z <- rnorm(T, mean = true_labels * mu, sd = 1)
                          z_burnin <- rnorm(burnin)
                          list(z_burnin = z_burnin, z = z, labels = true_labels)
                        },
                        linear = function(T = 5000, burnin = 500, mu = 2.5) {
                          pi_t <- seq(0.05, 0.2, length.out = T)
                          true_labels <- rbinom(T, 1, pi_t)
                          z <- rnorm(T, mean = true_labels * mu, sd = 1)
                          z_burnin <- rnorm(burnin)
                          list(z_burnin = z_burnin, z = z, labels = true_labels)
                        },
                        sine = function(T = 5000, burnin = 500, mu = 2.5) {
                          pi_t <- 0.1 + 0.05 * sin(2 * pi * (1:T) / 1000)
                          pi_t <- pmin(pmax(pi_t, 0), 1)
                          true_labels <- rbinom(T, 1, pi_t)
                          z <- rnorm(T, mean = true_labels * mu, sd = 1)
                          z_burnin <- rnorm(burnin)
                          list(z_burnin = z_burnin, z = z, labels = true_labels)
                        },
                        stop("Unknown mode.")
)

# SAST 和 onlineFDR的baseline 方法 
run_SAST <- function(z_burnin, z) {
  z_total <- c(z_burnin, z)
  res <- SAST(z_total, alpha = 0.05, conservative = TRUE, init = length(z_burnin), h = 50)
  decision <- res$decision[(length(z_burnin)+1):(length(z_burnin)+length(z))]
  return(decision)
}

run_onlineFDR <- function(z, method) {
  pvals <- 2 * pnorm(-abs(z))
  df <- data.frame(pval = pvals)
  res <- switch(method,
                LOND = LOND(df),
                LORD = LORD(df),
                ADDIS = ADDIS(df),
                stop("Unknown method.")
  )
  return(res$R)
}

evaluate_performance <- function(decision, true_labels) {
  TP <- sum(decision == 1 & true_labels == 1)
  FP <- sum(decision == 1 & true_labels == 0)
  FN <- sum(decision == 0 & true_labels == 1)
  FDP <- ifelse((TP + FP) == 0, 0, FP / (TP + FP))
  Power <- ifelse(sum(true_labels) == 0, 0, TP / sum(true_labels))
  MDR <- ifelse(sum(true_labels) == 0, 0, FN / sum(true_labels))
  c(FDP = FDP, Power = Power, MDR = MDR)
}

# 仿真主程序
T <- 5000
burnin <- 500
reps <- 100
checkpoint_times <- seq(1500, 5000, by = 500)
methods <- c("SAST", "LOND", "LORD", "ADDIS")

results <- list()
for (mtd in methods) {
  results[[mtd]] <- list(
    FDP = matrix(0, reps, length(checkpoint_times)),
    Power = matrix(0, reps, length(checkpoint_times)),
    MDR = matrix(0, reps, length(checkpoint_times))
  )
}

set.seed(123)

for (r in 1:reps) {
  sim <- simulate_data(T, burnin, mu = 2.5)
  sast_decision <- run_SAST(sim$z_burnin, sim$z)
  lond_decision <- run_onlineFDR(sim$z, "LOND")
  lord_decision <- run_onlineFDR(sim$z, "LORD")
  addis_decision <- run_onlineFDR(sim$z, "ADDIS")
  
  decisions <- list(SAST = sast_decision, LOND = lond_decision,
                    LORD = lord_decision, ADDIS = addis_decision)
  
  for (i in seq_along(checkpoint_times)) {
    t_point <- checkpoint_times[i]
    for (mtd in methods) {
      perf <- evaluate_performance(decisions[[mtd]][1:t_point], sim$labels[1:t_point])
      results[[mtd]]$FDP[r,i] <- perf["FDP"]
      results[[mtd]]$Power[r,i] <- perf["Power"]
      results[[mtd]]$MDR[r,i] <- perf["MDR"]
    }
  }
}

# 绘图
par(mfrow = c(1,3))
plot(checkpoint_times, colMeans(results$SAST$FDP), type = "l", ylim = c(0,0.1),
     main = paste("FDP vs Time (", mode, " Pattern)"), xlab = "Time", ylab = "FDP", col = "black", lwd = 2)
lines(checkpoint_times, colMeans(results$LOND$FDP), col = "red", lwd = 2)
lines(checkpoint_times, colMeans(results$LORD$FDP), col = "blue", lwd = 2)
lines(checkpoint_times, colMeans(results$ADDIS$FDP), col = "green", lwd = 2)
abline(h = 0.05, lty = 2)
legend("topright", legend = methods, col = c("black", "red", "blue", "green"), lwd = 2)

plot(checkpoint_times, colMeans(results$SAST$Power), type = "l", ylim = c(0,1),
     main = paste("Power vs Time (", mode, " Pattern)"), xlab = "Time", ylab = "Power", col = "black", lwd = 2)
lines(checkpoint_times, colMeans(results$LOND$Power), col = "red", lwd = 2)
lines(checkpoint_times, colMeans(results$LORD$Power), col = "blue", lwd = 2)
lines(checkpoint_times, colMeans(results$ADDIS$Power), col = "green", lwd = 2)
legend("topright", legend = methods, col = c("black", "red", "blue", "green"), lwd = 2)

plot(checkpoint_times, colMeans(results$SAST$MDR), type = "l", ylim = c(0,1),
     main = paste("MDR vs Time (", mode, " Pattern)"), xlab = "Time", ylab = "MDR", col = "black", lwd = 2)
lines(checkpoint_times, colMeans(results$LOND$MDR), col = "red", lwd = 2)
lines(checkpoint_times, colMeans(results$LORD$MDR), col = "blue", lwd = 2)
lines(checkpoint_times, colMeans(results$ADDIS$MDR), col = "green", lwd = 2)
legend("bottomright", legend = methods, col = c("black", "red", "blue", "green"), lwd = 2)

# 输出总结表格
summary_table <- data.frame(
  Method = methods,
  Mean_FDR = sapply(methods, function(m) mean(results[[m]]$FDP[, length(checkpoint_times)])),
  Mean_Power = sapply(methods, function(m) mean(results[[m]]$Power[, length(checkpoint_times)])),
  Mean_MDR = sapply(methods, function(m) mean(results[[m]]$MDR[, length(checkpoint_times)]))
)
print(summary_table)
