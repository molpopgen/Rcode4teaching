#' Trajectory of a new neutral mutation in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @return The allele frequency in each generation until fixation or loss
#' @examples
#' #Fate of a new mutant in a population of N = 100 diploids
#' n.traj = neutral.trajectory(100)
#' @export
neutral.trajectory <- function(N)
{
         g=1
         x=1
         trajectory = array()
         trajectory[g] = x/(2*N)
         while(  x> 0 && x < (2*N) )
         {
                  g=g+1
                  x=rbinom(1,(2*N),x/(2*N))
                  trajectory[g]=x/(2*N)
         }
         return(trajectory)
}

#' Trajectory of a new  mutation subject to genic/additive selection in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @param s The selection coefficient
#' @return The allele frequency in each generation until fixation or loss
#' @details Fitnesses are 1, 1+s, and 1+2s for genotypes AA, Aa, and aa, respectively.  Thus, a is the mutant allele.
#' @examples
#' g.traj = textbook.genic(100,0.01)
#' @export
textbook.genic <- function(N,s)
{
  return( rtraj(1/(2*N),N,1,1+s,1+2*s) )
}

#' Trajectory of a new  mutation subject to selection with arbitrary dominance in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The dominance of the mutant allele
#' @return The allele frequency in each generation until fixation or loss
#' @details Fitnesses are 1, 1+hs, and 1+2s for genotypes AA, Aa, and aa, respectively.  Thus, a is the mutant allele, and h = 1 is the genic/additive/codominant case.
#' @examples
#' d.traj = textbook.dominance(100,0.01,1)
#' @export
textbook.dominance <- function(N,s,h)
{
  return ( rtraj(1/(2*N),N,1,1+s*h,1+2*s) )
}

#' Trajectory of a new  mutation subject to overdomimant selection in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @param s1 The selection coefficient applied to AA
#' @param s2 The selection coefficient applied to aa
#' @param ifreq The initial frequency of the a allele.
#' @param ngen The number of generations to iterator the function.
#' @return The allele frequency in each generation until fixation or loss
#' @details Fitnesses are 1-s1, 1, and 1-s2 for genotypes AA, Aa, and aa, respectively.
#' @examples
#' od.traj = textbook.overdominance(100,0.1,0.1,0.2,100)
#' @export
textbook.overdominance <- function(N,s1,s2,ifreq,ngen)
{
    if(ngen <= 0)
        {
            stop("ngen must be > 0")
        }
    return ( rtraj(ifreq,N,1-s2,1,1-s1,ngen) )
}

#' Estimate fixation probability of a new mutation under genic/additive selection
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param howmany The number of simulation to run
#' @return An estimate of the fixation probability
#' @details Fitnesses are 1, 1+s, and 1+2s for genotypes AA, Aa, and aa, respectively.  Thus, a is the mutant allele.
#' @examples
#' pfix.genic = pfixgenic(100,0.1,10)
#' @export
pfixgenic <- function(N,s,howmany)
{
  return( pfix(1/(2*N),N,1,1+s,1+2*s,howmany) )
}

#' Estimate fixation probability of a new mutation under models with arbitrary dominance
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The domimance of the mutant allele
#' @param howmany The number of simulation to run
#' @return An estimate of the fixation probability
#' @details Fitnesses are 1, 1+hs, and 1+2s for genotypes AA, Aa, and aa, respectively.  Thus, a is the mutant allele.
#' @examples
#' pfix.dom = pfixdominance(100,0.1,1,10)
#' @export
pfixdominance <- function(N,s,h,howmany)
{
    return( pfix(1/(2*N),N,1,1+s*h,1+2*s,howmany) )
}

#' Returns a trajectory of a mutation that has fixed due to selection
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The domimance of the mutant allele
#' @return A list with the allele frequency trakectory plus the number of attempts needed before a mutation fixed
#' @details Fitnesses are 1, 1+hs, and 1+2s for genotypes AA, Aa, and aa, respectively.  Thus, a is the mutant allele.
#' @examples
#' ft = fixtrajdominance(100,0.1,0.05)
#' @export
fixtrajdominance <- function(N,s,h)
{
    ntries=0
    while( 1 )
        {
            ntries=ntries+1
            t<-textbook.dominance(N,s,h)
            if(t[length(t)]==1)
		{
                    return(list(traj=t,n=ntries))
		}
	}
}
