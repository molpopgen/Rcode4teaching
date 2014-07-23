#' Trajectory of a new neutral mutation in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @return The allele frequency in each generation until fixation or loss
#' @export
WFtrajectory <- function(N)
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
#' @export
WFgenic <- function(N,s)
{
  return( rtraj(1/(2*N),N,1,1+s,1+2*s) )
}

#' Trajectory of a new  mutation subject to selection with arbitrary dominance in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The dominance of the mutant allele
#' @return The allele frequency in each generation until fixation or loss
#' @export
WFdominance <- function(N,s,h)
{
  return ( rtraj(1/(2*N),N,1,1+s*h,1+2*s) )
}

#' Trajectory of a new  mutation subject to overdomimant selection in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @param s1 The selection coefficient applied to AA
#' @param s2 The selection coefficient applied to aa
#' @param h The dominance of the mutant allele
#' @return The allele frequency in each generation until fixation or loss
#' @export
WFoverdominance <- function(N,s1,s2,ngen,ifreq)
{
  return ( rtraj(1/(2*N),N,1-s2,1,1-s1) )
}

#' Estimate fixation probability of a new mutation under genic/additive selection
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param howmany The number of simulation to run
#' @return An estimate of the fixation probability
#' @export
pfixgenic <- function(N,s,howmany)
{
  return( pfix(1/(2*N),N,1,1+s,1+2*s,howmany) )
  ## i=0
  ## nfix = 0;
  ## while( i < howmany )
  ## {
  ## 	i=i+1
  ## 	t<-WFgenic(N,s)
  ## 	if ( t[length(t)] == 1 )
  ## 	{
  ## 		nfix=nfix+1
  ## 	}
  ## }
  ## return(nfix/howmany)
}

#' Estimate fixation probability of a new mutation under models with arbitrary dominance
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The domimance of the mutant allele
#' @param howmany The number of simulation to run
#' @return An estimate of the fixation probability
#' @export
pfixdominance <- function(N,s,h,howmany)
{
  return( pfix(1/(2*N),N,1,1+s*h,1+2*s,howmany) )
	## i=0
	## nfix = 0;
	## while( i < howmany )
	## {
	## 	i=i+1
        ##         t<-WFdominance(N,s,h)
	## 	if ( t[length(t)] == 1 )
	## 	{
	## 		nfix=nfix+1
	## 	}
	## }
	## return(nfix/howmany)
}

#' Returns a trajectory of a mutation that has fixed due to selection
#' @param N The diploid population size
#' @param s The selection coefficient
#' @param h The domimance of the mutant allele
#' @return A list with the allele frequency trakectory plus the number of attempts needed before a mutation fixed
#' @export
fixtrajdominance <- function(N,s,h)
{
	ntries=0
	while( 1 )
	{
		ntries=ntries+1
		t<-WFdominance(N,s,h)
		if(t[length(t)]==1)
		{
			return(list(traj=t,n=ntries))
		}
	}
}
