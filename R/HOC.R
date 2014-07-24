#' @export
HOC.VGVE = function(pop)
    {
        return( list(VG=var(pop[,1]+pop[,2]),VE=var(pop[,3])) )
    }

#' Distribution of VG and VE under Gaussian HOC model
#' @param N The number of diploids
#' @param mu The mutation rate (per gamete, per generation)
#' @param sigmamu The standard deviation in the effect size of new mutations
#' @param sigmae The standard deviation in random effects added to phenotype
#' @param sigmas The standard deviation in the Gaussuian fitness function
#' @param ngens The number of generations to simulate
#' @param nreps The number of replicates to simulate
#' @return A data frame with the VG and VE components of phenotype from each replicate
#' @details uses replicate to call HOCsim repeatedly
#' @export
HOC.VGVE.dist = function(N,mu,sigmaMU,sigmaE,sigmaS,ngens,nreps)
    {
        x=( replicate(nreps,HOC.VGVE(HOCsim(N,mu,sigmaMU,sigmaE,sigmaS,ngens))) )
        return( data.frame(VG=as.numeric(x['VG',]),VE=as.numeric(x['VE',])) )
    }
