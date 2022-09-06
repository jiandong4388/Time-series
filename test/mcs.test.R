### Model Confident Set



mcs.test <- function(losses, alpha = 0.05, nboot = 100, nblock = 1, boot = c("stationary", "block"))
{
    n = NROW(losses)
    m = NCOL(losses)
    bsdat = bootstrap(1:n, nboot, nblock, type = boot[1])$B
    dijbar = matrix(0, m, m)
    for(j in 1:m) dijbar[j, ] = colMeans(losses - losses[,j])
    dijbarstar = array(0, dim = c(m, m, nboot))
    for(i in 1:nboot){
        tmp = colMeans(losses[bsdat[,i],,drop = FALSE])
        for(j in 1:m){
            dijbarstar[j,,i] = tmp - tmp[j]
        }
    }
    tmp = array(NA, dim = c(m, m, nboot))
    for(i in 1:nboot) tmp[,,i] = dijbar
    xtmp = ((dijbarstar - tmp)^2)
    vardijbar = matrix(0, m, m)
    for(i in 1:nboot) vardijbar = vardijbar + xtmp[,,i]
    vardijbar = vardijbar/nboot
    vardijbar = vardijbar + diag(m)
    tmp2 = array(NA, dim = c(m, m, nboot))
    for(i in 1:nboot) tmp2[,,i] = sqrt(vardijbar)
    z0=(dijbarstar-tmp)/tmp2
    zdata0=dijbar/sqrt(vardijbar)
    excludedR = matrix(0, m,1)
    pvalsR = matrix(1,m,1)
    for(i in 1:(m-1)){
        included = setdiff(1:m,excludedR)
        mx = length(included)
        z = z0[included, included,]
        empdistTR = apply(abs(z), 3, function(x) max(x))
        zdata = zdata0[included,included]
        TR = max(.Call("colMaxRcpp", zdata, PACKAGE = "rugarch"))
        pvalsR[i] = mean(empdistTR>TR)
        dibar = colMeans(dijbar[included,included])*(mx/(mx-1))
        dibstar = apply(dijbarstar[included,included,,drop=FALSE], 3, function(x) colMeans(x))*(mx/(mx-1))
        vardi=colMeans((t(dibstar)-matrix(dibar, nrow = nboot, ncol = mx, byrow = TRUE))^2)
        tx=dibar/sqrt(vardi)
        temp = max(tx)
        modeltoremove = which(tx == max(tx))[1]
        excludedR[i]=included[modeltoremove]
    }
    maxpval=pvalsR[1]
    for(i in 2:m){
        if(pvalsR[i]<maxpval){
            pvalsR[i] = maxpval
        } else{
            maxpval = pvalsR[i]
        }
    }
    excludedR[NROW(excludedR)] = setdiff(1:m,excludedR)
    pl = which(pvalsR >= alpha)[1]
    includedR = excludedR[pl:m]
    if(pl==1) excludedR = NULL else excludedR = excludedR[1:(pl-1)]
    excludedSQ = matrix(0,m,1)
    pvalsSQ = matrix(1,m,1)
    for(i in 1:(m-1)){
        included = setdiff(1:m,excludedSQ)
        mx = length(included)
        z = z0[included,included,]
        empdistTSQ = apply(z^2, 3, function(x) sum(x))/2
        zdata = zdata0[included,included]
        TSQ = sum(colSums(zdata^2))/2
        pvalsSQ[i] = mean(empdistTSQ>TSQ)
        dibar = colMeans(dijbar[included,included])*(mx/(mx-1))
        dibstar = apply(dijbarstar[included,included,,drop=FALSE], 3, function(x) colMeans(x))*(mx/(mx-1))
        vardi = colMeans((t(dibstar)-matrix(dibar, nrow = nboot, ncol = mx, byrow = TRUE))^2)
        tx = dibar/sqrt(vardi)
        temp = max(tx)
        modeltoremove = which(tx == max(tx))[1]
        excludedSQ[i]=included[modeltoremove]
    }
    maxpval = pvalsSQ[1]
    for( i in 2:m){
        if(pvalsSQ[i]<maxpval){
            pvalsSQ[i] = maxpval
        } else{
            maxpval = pvalsSQ[i]
        }
    }
    excludedSQ[NROW(excludedSQ)] = setdiff(1:m,as.numeric(excludedSQ))
    pl=which(pvalsSQ>=alpha)[1]
    includedSQ = excludedSQ[pl:m]
    if(pl==1) excludedSQ = NULL else excludedSQ = excludedSQ[1:(pl-1)]
    return(list(includedR = includedR, pvalsR = pvalsR, excludedR = excludedR,
                includedSQ = includedSQ, pvalsSQ = pvalsSQ, excludedSQ = excludedSQ))
}
