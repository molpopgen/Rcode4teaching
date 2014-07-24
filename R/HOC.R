#' @export
HOC.VGVE = function(pop)
    {
        return( list(VG=var(pop[,1]+pop[,2]),VE=var(pop[,3])) )
    }

#' @export
HOC.VGVE.dist = function(N,mu,sigmaMU,sigmaE,sigmaS,ngens,nreps)
    {
        x=( replicate(nreps,HOC.VGVE(HOCsim(N,mu,sigmaMU,sigmaE,sigmaS,ngens))) )
        return( data.frame(VG=as.numeric(x['VG',]),VE=as.numeric(x['VE',])) )
    }
