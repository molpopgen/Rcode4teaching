#Standard population-genetic models of selection in two-allele systems
#deltaq can be used to plot change in allele freq in infinite-N pops for autosomal loci
#rtraj and pfix can be used to relate simulation to results from stochastic theory

#' General formula for delta q for 2-allele model
#' @param q Current allele frequency
#' @param wAA Fitness of AA genotype
#' @param wAa Fitness of Aa genotype
#' @param waa Fitness of aa genotype
#' @return The change in allele frequency in an infinitely-large Wright-Fisher population
deltaq = function(q,wAA,wAa,waa)
{
  p = 1 - q
  wbar = wAA*p^2 + 2*p*q*wAa + waa*q^2
  return ( p*q*(p*(wAa-wAA)+q*(waa-wAa))/wbar )
}

#' Simulate random trajectory under selection starting from q0.
#' @param q0 initial mutant allele frequency
#' @param N Population size
#' @param wAA Fitness of AA genotype
#' @param wAa Fitness of Aa genotype
#' @param waa Fitness of aa genotype
#' @return A random allele frequency trajectory
rtraj = function(q0,N,wAA,wAa,waa)
  {
    q=q0
    traj=array()
    traj[1]=q
    generation=2
    while( q > 0 & q < 1 )
      {
        eqp = q + deltaq(q,wAA,wAa,waa)
        q = rbinom(1,2*N,eqp)/(2*N)
        traj[generation]=q
        generation=generation+1
      }
    return (traj)
  }

#' Estimate fixation prob and time to fixation by simulation
#' @param q0 initial mutant allele frequency
#' @param N Population size
#' @param wAA Fitness of AA genotype
#' @param wAa Fitness of Aa genotype
#' @param waa Fitness of aa genotype
#' @param nsims The number of replicates to run
#' @return A list.  prob = estimate of fixation probability.  times = distribution of fixation times
pfix = function(q0,N,wAA,wAa,waa,nsims)
  {
    times=array()
    dummy=0
    trials = 0
    while( trials < nsims )
      {
        t = rtraj(q0,N,wAA,wAa,waa)
        if( t[length(t)] == 1 ) #fixed
          {
            times[dummy]=length(t) #time in generations
            dummy = dummy + 1
          }
        trials = trials + 1
      }
    return(list(prob=dummy/nsims,
                times=times))
  }
