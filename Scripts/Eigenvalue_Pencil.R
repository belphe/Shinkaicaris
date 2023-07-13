## Eigenvalue/Pencil
## https://www.sciencedirect.com/science/article/pii/S1550413117302905
## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0198503
## Customized Script by Hongyin Zhang

library(MASS)
library(dplyr)

data <- read.csv("Ceemdan_IMFs_mean.csv", header = TRUE, row.names = 1)

r <- 3 # the dimension of the reduced system, e.g r = 3, r = 5, r = 7 etc.
resolution <- 4 #sampling resolution/interval (hours)

df <- c()
for (i in 1:length(rownames(data))){
  matrix <- as.matrix(t(data[i,2:19])) # This method requires the length of a given time series to be even, so we remove the data of time-point T00
  k <- length(matrix)/2
  H <- embed(matrix,k+1)[1:k,(k+1):1]
  E <- H[1:k,1:k]
  A <- H[1:k,2:(k+1)]
  B <- H[1:k,1]
  C <- H[1,1:k]
  SVD1 <- svd(rbind(E,A))
  SVD2 <- svd(cbind(E,A))
  X <- SVD2$u[1:k,1:r]
  Y <- SVD1$v[1:k,1:r]
  Er <- t(X)%*%E%*%Y
  Ar <- t(X)%*%A%*%Y
  Br <- t(X)%*%B
  Cr <- C%*%Y
  Brr <- ginv(Er)%*%Br
  EVD <- eigen(ginv(Er)%*%Ar)
  L <- EVD$values
  V <- EVD$vectors
  V_inverse <- ginv(V)
  Period_col <- c()
  Amplitude_col <- c()
  Phase_circle_col <- c()
  Phase_col <- c()
  Decay_rate_col <- c()
  EP_result <- c()
  for (j in 1:r){
  P <- (Cr%*%V[,j])*(V_inverse[j,]%*%Brr)
  period <- abs(resolution*2*pi/Im(log(L[j])))
  amplitude <- 4*exp(Re(log(P))) # peak to trough amplitude
  phase_circle <- ((360-Im(log(P))*180/pi)%%360)
  phase <- ((1-(Im(log(P))/(2*pi)))%%1)*period
  decay_rate <- 1 + Re(log(L[j]))
  Period_col <- rbind(Period_col, period)
  Amplitude_col <- rbind(Amplitude_col, amplitude)
  Phase_circle_col <- rbind(Phase_circle_col, phase_circle)
  Phase_col <- rbind(Phase_col, phase)
  Decay_rate_col <- rbind(Decay_rate_col, decay_rate)
  }
  EP_result <- as.data.frame(cbind(Amplitude_col, Phase_circle_col, Phase_col, Decay_rate_col, Period_col))
  colnames(EP_result) <- c("Amplitude", "Phase_circle", "Phase", "Decay_rate", "Period")
  EP_result <- distinct(EP_result, Period, .keep_all = TRUE)
  EP_result <- EP_result[order(EP_result$r_Amplitude, decreasing = TRUE),]
  ID <- c()
  if (length(EP_result$Period) == 0) {
    next
  } else {
    for(n in 1:length(EP_result$Period)){
      ID <- rbind(ID, paste0(rownames(raw[i,]), sep = ".", n))
    }
  rownames(EP_result) <- ID
  df <- rbind(df, EP_result)
  }
}

#Decay_rate_percentage10 <- quantile(df[,5], 0.1)
#Decay_rate_percentage90 <- quantile(df[,5], 0.9)
Decay_rate_percentage_left <- 0.8
Decay_rate_percentage_right <- 1.2
Period_threshold_left <- 8 # twice the sampling resolution
Period_threshold_right <- 36 # half of the sampling duration

Rhythmicity <- c()
for (i in 1:length(rownames(df))){
  if (df[i,"Decay_rate"] < Decay_rate_percentage_left | df[i,"Decay_rate"] > Decay_rate_percentage_right | df[i,"Period"] <= Period_threshold_left | df[i,"Period"] > Period_threshold_right | df[i,"Period"] == "Inf") {
    rhythmicity <- "nonrhythmic"
  } else {
    rhythmicity <- "rhythmic"
  }
  Rhythmicity <- rbind(Rhythmicity, rhythmicity)
}
df_R <- cbind(df, Rhythmicity)

write.csv(df_R, "EP_T04-T72_r3.csv")


