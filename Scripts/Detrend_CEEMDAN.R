## To detrend the raw TMM expression data using Complete Ensemble Empirical Mode Decomposition with Adaptive Noise (CEEMDAN)
## https://ieeexplore.ieee.org/document/5947265/
## https://github.com/helske/Rlibeemd
## We select the intrinsic mode functions (IMFs) to reconstruct the detrended data by calculating the correlation coefficients (cor) and variance contribution rates (vcr), and then follow the criterion suggested by ALBERT AYENU-PRAH and NII ATTOH-OKINE (2010); https://doi.org/10.1142/S1793536910000367

library(Rlibeemd)

raw_mean <- as.matrix(read.csv("23905_raw_mean.csv", header = TRUE, row.names = 1))
raw <- read.csv("23905_TMM_raw.txt", row.names = 1)

length_N <- length(colnames(raw_mean))
N <- length(rownames(raw_mean))

threshold_df <- data.frame()
cor_value_df <- data.frame()
vcr_value_df <- data.frame()
IMFs_df <- data.frame()
Residual_df <- data.frame()

for (i in 1:N){
  ceemdan_result <- as.data.frame(ceemdan(raw_mean[i,], num_imfs = 0, noise_strength = 0.4, S_number = 0, num_siftings = 50, ensemble_size = 250))
  cor_IMF1 <- cor.test(x = raw_mean[i,], y = ceemdan_result$`IMF 1`)$estimate
  cor_IMF2 <- cor.test(x = raw_mean[i,], y = ceemdan_result$`IMF 2`)$estimate
  cor_IMF3 <- cor.test(x = raw_mean[i,], y = ceemdan_result$`IMF 3`)$estimate
  threshold <- max(c(cor_IMF1, cor_IMF2, cor_IMF3))/(10*max(c(cor_IMF1, cor_IMF2, cor_IMF3))-3)
  threshold_df <- rbind(threshold_df, threshold)
  cor_value <- c(cor_IMF1, cor_IMF2, cor_IMF3)
  cor_value_df <- rbind(cor_value_df, cor_value)
  vcr <- var(ceemdan_result$`IMF 1`) + var(ceemdan_result$`IMF 2`) + var(ceemdan_result$`IMF 3`)
  vcr_IMF1 <- var(ceemdan_result$`IMF 1`)/vcr
  vcr_IMF2 <- var(ceemdan_result$`IMF 2`)/vcr
  vcr_IMF3 <- var(ceemdan_result$`IMF 3`)/vcr
  vcr_value <- c(vcr_IMF1, vcr_IMF2, vcr_IMF3)
  vcr_value_df <- rbind(vcr_value_df, vcr_value)
  order_vcr <- order(c(vcr_IMF1, vcr_IMF2, vcr_IMF3), decreasing = TRUE)
  components_value <- rbind(ceemdan_result$`IMF 1`, ceemdan_result$`IMF 2`, ceemdan_result$`IMF 3`)
  vcr_sum <- 0
  union_ceemdan <- 0
  max_vcr <- order_vcr[1]
  if (cor_value[max_vcr] > threshold & vcr_value[max_vcr] >= 2/3) {
    union_ceemdan <- components_value[max_vcr,]
  } else {
    for (n in order_vcr){
      if (cor_value[n] > threshold & vcr_sum < 2/3){
        union_ceemdan <- union_ceemdan + components_value[n,]
        vcr_sum <- vcr_sum + vcr_value[n]
      }
    }
  }
  if (union_ceemdan == 0){
    union_ceemdan <- ceemdan_result$`IMF 1` + ceemdan_result$`IMF 2` + ceemdan_result$`IMF 3`
  }
  IMFs_df <- rbind(IMFs_df, union_ceemdan)
  Residual_df <- rbind(Residual_df, ceemdan_result$Residual)
}

cor_value_df <- cbind(cor_value_df, threshold_df)
colnames(cor_value_df) <- c("cor_IMF1", "cor_IMF2", "cor_IMF3", "thershold")
rownames(cor_value_df) <- rownames(raw_mean)
write.csv(cor_value_df, "Ceemdan_cor.csv")

colnames(vcr_value_df) <- c("vcr_IMF1", "vcr_IMF2", "vcr_IMF3")
rownames(vcr_value_df) <- rownames(raw_mean)
write.csv(vcr_value_df, "Ceemdan_vcr.csv")

colnames(IMFs_df) <- colnames(raw_mean)
rownames(IMFs_df) <- rownames(raw_mean)
write.csv(IMFs_df, "Ceemdan_IMFs_mean.csv")

colnames(Residual_df) <- colnames(raw_mean)
rownames(Residual_df) <- rownames(raw_mean)
write.csv(Residual_df, "Ceemdan_residual.csv")

raw_mean <- as.data.frame(raw_mean)
ceemdan_IMFs_df <- data.frame()
for (i in 1:N){
  deviation_IMFs <- as.data.frame(rep(raw_mean[i,] - IMFs_df[i,], each = 3)) 
  ceemdan_IMFs_df <- rbind(ceemdan_IMFs_df, raw[i,] - deviation_IMFs) 
  } 
write.csv(ceemdan_IMFs_df, "Ceemdan_IMFs.csv")





