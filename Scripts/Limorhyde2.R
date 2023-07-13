## Limorhyde2
## https://limorhyde2.hugheylab.org/

library(limorhyde2)
library(data.table)

data <- read.csv("Ceemdan_IMFs.csv", header = TRUE, row.names = 1)

for (i in 1:nrow(data)){data[i,] <- data[i,] + abs(min(data[i,]))}

sample <- colnames(data)
time <- rep(seq(0,72,4), each = 3)
design <- data.table(cbind(sample, time))
time <- apply(design[,2],2,as.numeric)
design <- data.table(design$sample, time)
colnames(design) <- c("sample","time")

## Fit linear models for rhythmicity with a period of 12 hours
fit <- getModelFit(data, design, period = 12, nKnots = 3, degree = 3, method = "voom")
fit <- getPosteriorFit(fit, covMethod = "data-driven")
rhyStats <- getRhythmStats(fit, rms = TRUE)
write.csv(rhyStats, "limorhyde2_p12.csv")

## Fit linear models for rhythmicity with a period of 24 hours
fit <- getModelFit(data, design, period = 24, nKnots = 3, method = "voom")
fit <- getPosteriorFit(fit, covMethod = "data-driven")
rhyStats <- getRhythmStats(fit, rms = TRUE)
write.csv(rhyStats, "limorhyde2_p24.csv")
