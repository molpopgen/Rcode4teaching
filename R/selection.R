#Standard population-genetic models of selection in two-allele systems
#deltaq can be used to plot change in allele freq in infinite-N pops for autosomal loci
#rtraj and pfix can be used to relate simulation to results from stochastic theory

#general formula for delta q for 2-allele model
deltaq = function(q,wAA,wAa,waa)
{
  p = 1 - q
  wbar = wAA*p^2 + 2*p*q*wAa + waa*q^2
  return ( p*q*(p*(wAa-wAA)+q*(waa-wAa))/wbar )
}

#random trajectory under selection starting from q0
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

#fixation prob and dist of times
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
