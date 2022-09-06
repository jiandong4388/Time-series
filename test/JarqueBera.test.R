## Jarque–Bera test is a goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution
## 
##  Jarque, Carlos M.; Bera, Anil K. (1980). "Efficient tests for normality, homoscedasticity and serial independence of regression residuals". Economics Letters.  
##   

## Example: 
##   x <- rnorm(1000,2,3)
##   JarqueBera.test(x,alpha=0.05)


JarqueBera.test <- function(data, alpha = 0.05) {
  
  iT = length(data)
  
  m1 = mean(data)
  m2 = mean((data - m1)^2.0)
  m3 = mean((data - m1)^3.0)
  m4 = mean((data - m1)^4.0)
  
  dS = m3 / m2^(3.0 / 2.0)
  dK = m4 / m2^2.0
  
  dStat = iT / 6.0 * (dS^2.0 + (dK - 3.0)^2.0 / 4.0)
  
  dPval = 1.0 - pchisq(dStat, 2)
  
  dCritical = qchisq(1.0 - alpha, 2L)
  
  return(c("Statistic" = dStat, "critical" = dCritical, "p-Value" = dPval))
  
}
