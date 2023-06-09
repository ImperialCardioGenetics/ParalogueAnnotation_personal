rm(list = ls())
graphics.off()

# Set seed for reproducibility
set.seed(123)

### bootrstrapping =========

# create output vectors
vec_ratio_A <- rep(NA, 10000)
vec_ratio_B <- rep(NA, 10000)
vec_dif <- rep(NA, 10000)

# resampling 10000 times
for (i in 1:10000){
  
  # calculate ratios by sampling with the observed sample size and lamda values
  vec_ratio_A[i] <- sum(rpois(2508, 1650/2508)) / sum(rpois(1911, 1107/1911))
  vec_ratio_B[i] <- sum(rpois(2508, 9/2508)) / sum(rpois(1911, 2/1911))
  
  # calculate the difference between ratios
  vec_dif[i] <- vec_ratio_A[i] - vec_ratio_B[i]
  # print(i)
}

### calculate confidence intervals 

# Sort the differences
vec_dif <- sort(vec_dif)

# Find the 2.5th and 97.5th percentiles
lower_bound <- quantile(vec_dif, 0.025, names = FALSE)
upper_bound <- quantile(vec_dif, 0.975, names = FALSE)

# Display the 95% confidence interval
cat("95% Confidence Interval: [", lower_bound, ", ", upper_bound, "]")
# if confidence interval does not overlap 0 you can reject the null hypothesis. 
