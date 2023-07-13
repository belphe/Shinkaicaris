## Determine whether or not there is a the monotonic trend in each raw transcript with the Mann-Kendall Test
## https://cran.r-project.org/web/packages/Kendall/index.html

library(Kendall)

raw_mean <- as.matrix(read.csv("Raw_TMM_mean.csv", header = TRUE, row.names = 1))

raw_mkt_df <- data.frame()

for (i in 1:nrow(raw_mean)) {
  raw_mkt_result <- MannKendall(raw_mean[i,])
  if (raw_mkt_result$S > 0) {
    standard_normal_variate = (raw_mkt_result$S - 1)/sqrt(raw_mkt_result$varS)
    } 
  else if (raw_mkt_result$S == 0) {
    standard_normal_variate = 0
    }
  else if (raw_mkt_result$S < 0) {
    standard_normal_variate = (raw_mkt_result$S + 1)/sqrt(raw_mkt_result$varS)
    }
  raw_mkt_pvalue <- raw_mkt_result$sl
  if (raw_mkt_pvalue >= 0.05) {
    trend <- "Constant"
  } else {
    if (standard_normal_variate < 0) {
      trend <- "Downward"
    } else if (standard_normal_variate > 0) {
      trend <- "Upward"
    }
  }
  raw_mkt_df <- rbind(raw_mkt_df, cbind(raw_mkt_pvalue, standard_normal_variate, trend, t(raw_mean[i,])))
}

raw_mkt_df <- cbind(rownames(raw_mean), raw_mkt_df)
colnames(raw_mkt_df)[1] <- "ID"
write.csv(raw_mkt_df, "Raw_TMM_MKT.csv")

#raw_mkt_df_significant <- raw_mkt_df[which(raw_mkt_df$raw_mkt_pvalue < 0.05),]
#write.csv(raw_mkt_df_significant, "Raw_TMM_MKT_significant.csv")
