LjungBox.test <- function(data, lag = c(10, 15, 20)) {
    
    iP    = length(lag)
    mTest = matrix(data = NA, iP, 2L, dimnames =
                       list(lag, c("Statistic", "p-Value")))
    
    for (p in 1:iP) {
        Test = Box.test(data, lag[p], type = "Ljung-Box")
        mTest[p, ] = c(Test$statistic, Test$p.value)
    }
    
    return(mTest)
    
}
