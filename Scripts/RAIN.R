## RAIN
## https://rdrr.io/bioc/rain/
## detect rhythmic transcripts with asymmetric waveform and a periodicity of 8, 12, 16, 20, 24, and 28 hours

library(rain)

data <- read.csv("Ceemdan_IMFs.csv", header = T, row.names = 1)

for (i in seq(8,28,4)){
  rain_result <- rain(t(data), 
                      period = i, 
                      deltat = 4, 
                      period.delta = 0, 
                      nr.series = 3, 
                      peak.border = c(0,1),
                      method = "independent",
                      adjp.method = "ABH",
                      verbose = FALSE
                      )
  rain_result_ordered <- rain_result$pVal[order(rain_result$pVal)]
  filter_num <- max(which(rain_result_ordered < 0.05))
  rain_result_filtered <- rain_result[which(rain_result$pVal %in% rain_result_ordered[1:filter_num]),]
  rain_result_filtered_ordered <- rain_result_filtered[order(rain_result_filtered$pVal),]
  write.csv(rain_result_filtered_ordered, file = paste0("rain_p",i,"_ABH005_",filter_num,".csv"))
  }


