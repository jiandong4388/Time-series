JarqueBera.test <- function(vRes, dAlpha = 0.05) {
  
  iT = length(vRes)
  
  m1 = mean(vRes)
  m2 = mean((vRes - m1)^2.0)
  m3 = mean((vRes - m1)^3.0)
  m4 = mean((vRes - m1)^4.0)
  
  dS = m3 / m2^(3.0 / 2.0)
  dK = m4 / m2^2.0
  
  dStat = iT / 6.0 * (dS^2.0 + (dK - 3.0)^2.0 / 4.0)
  
  dPval = 1.0 - pchisq(dStat, 2)
  
  dCritical = qchisq(1.0 - dAlpha, 2L)
  
  return(c("Statistic" = dStat, "p-Value" = dPval, "critical" = dCritical))
  
}
