## ECHO
## https://github.com/delosh653/ECHO

library(echo.find)

data <- read.csv("Ceemdan_IMFs.csv", header = TRUE)
colnames(data)[1] <- "ID"

begin <- 0 # first time point
end <- 72 # last time point
resol <- 4 # time point resolution
num_reps <- 3 # number of replicates
low <- 10 # low period seeking
#low <- 22 # low period seeking
high <- 14 # high period seeking
#high <- 26 # high period seeking
run_all_per <- FALSE # we are not looking for all periods
paired <- FALSE # these replicates are unrelated, that is, a replicate being 
# called "replicate 1" at time point 2 means nothing
rem_unexpr <- FALSE # do not remove unexpressed genes
# we do not assign rem_unexpr_amt, since we're not removing unexpressed genes
is_normal <- FALSE # do not normalize
is_de_linear_trend <- FALSE # do not remove linear trends
is_smooth <- FALSE # do not smooth the data
#run_conf <- TRUE # do run confidence intervals
#which_conf <- "Bootstrap" # string of which type of confidence interval to compute ("Bootstrap" or "Jackknife")

echo_result <- echo_find(genes = data, begin = begin, end = end, resol = resol, 
                         num_reps = num_reps, low = low, high = high, run_all_per = run_all_per,
                         paired = paired, rem_unexpr = rem_unexpr, is_normal = is_normal,
                         is_de_linear_trend = is_de_linear_trend, is_smooth = is_smooth)

echo_result$Initial.Amplitude <- abs(echo_result$Initial.Amplitude)
echo_result_ordered <- echo_result$`BH Adj P-Value`[order(echo_result$`BH Adj P-Value`)]
filter_num <- max(which(echo_result_ordered < 0.05))
echo_result_filtered <- echo_result[which(echo_result$`BH Adj P-Value` %in% echo_result_ordered[1:filter_num]),]
echo_result_filtered_ordered <- echo_result_filtered[order(echo_result_filtered$`BH Adj P-Value`),]

write.csv(echo_result_filtered_ordered[,c(1,4:16)], file = paste0("echo_p",low,"-",high,"_BH005_",filter_num,".csv"))

