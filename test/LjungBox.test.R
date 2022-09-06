LjungBox.test <- function(vRes, vLag = c(10, 15, 20)) {
    
    iP    = length(vLag)
    mTest = matrix(data = NA, iP, 2L, dimnames =
                       list(vLag, c("Statistic", "p-Value")))
    
    for (p in 1:iP) {
        Test = Box.test(vRes, vLag[p], type = "Ljung-Box")
        mTest[p, ] = c(Test$statistic, Test$p.value)
    }
    
    return(mTest)
    
}
