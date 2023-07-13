## Detect periodic singnals in the 351 days temperature records using Lomb-Scargle Periodogram analysis
## https://cran.r-project.org/package=lomb

library(lomb)

data <- read.csv("Temperature_peripheral.csv", header = TRUE)
temperature <- data$Temperature

my.date <- seq(0, (length(temperature)-1)*10/60, 10/60)

t.spec <- lsp(x = ts(temperature), times = my.date, ofac = 50, type = 'period', from = 2, to = 28, alpha = 1e-20)
sig.level <- t.spec$sig.level
plot.lsp(t.spec, main = NULL, xlabel = "Period (hours)", ylabel = "Normalized Power")
summary(t.spec)
peak_top15 <- getpeaks(t.spec, npeaks = 15)
peak_top15$data

peak_tide_constituent <- cbind(peak_top15$data[c(1,3,5,6,10,14),], c("M2", "S2", "O1", "M4","N2", "K1"))
colnames(peak_tide_constituent) <- c("time", "peak", "tide_constituent")
write.csv(peak_tide_constituent, "LS_peak.csv")

lsp_result <- as.data.frame(cbind(t.spec$scanned, t.spec$power))
colnames(lsp_result) <- c("scanned", "power")
write.csv(lsp_result, "LSP_results.csv")
