#' Trajectory of a new neutral mutation in a diploid Wright-Fisher population of size N
#' @param N The diploid population size
#' @return The allele frequency in each generation until fixation or loss
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

WFgenic <- function(N,s)
{
	wAA = 1+2*s
	wAa = 1+s
	waa = 1

	p = 1/(2*N)
	x = 1
	traj <- array();
	generation = 1
	traj[generation] = p
	while(p > 0 && p < 1)
	{
		q = 1-p
		wbar = p*p*wAA + 2*p*q*wAa + q*q*waa
		epprime = (p*p*wAA + p*q*wAa)/wbar
		x = rbinom(1,2*N,epprime)
		p = x/(2*N)
		generation = generation+1
		traj[generation] = p
	}
	return(traj)
}

WFdominance <- function(N,s,h)
{
	wAA = 1+s
	wAa = 1+s*h
	waa = 1

	p = 1/(2*N)
	x = 1
	traj <- array();
	generation = 1
	traj[generation] = p
	while(p > 0 && p < 1)
	{
		q = 1-p
		wbar = p*p*wAA + 2*p*q*wAa + q*q*waa
		epprime = (p*p*wAA + p*q*wAa)/wbar
		x = rbinom(1,2*N,epprime)
		p = x/(2*N)
		generation = generation+1
		traj[generation] = p
	}
	return(traj)
}

WFoverdominance <- function(N,s1,s2,ngen,ifreq)
{
	wAA = 1-s2
	wAa = 1
	waa = 1-s1

	p = ifreq
	x = as.integer(ifreq*2*N)
	traj <- array();
	generation = 1
	traj[generation] = p
	while(generation < ngen)
	{
		q = 1-p
		wbar = p*p*wAA + 2*p*q*wAa + q*q*waa
		epprime = (p*p*wAA + p*q*wAa)/wbar
		x = rbinom(1,2*N,epprime)
		p = x/(2*N)
		generation = generation+1
		traj[generation] = p
	}
	return(traj)
}

pfixgenic <- function(N,s,howmany)
{
	i=0
	nfix = 0;
	while( i < howmany )
	{
		i=i+1
		t<-WFgenic(N,s)
		if ( t[length(t)] == 1 )
		{
			nfix=nfix+1
		}
	}
	return(nfix/howmany)
}

pfixdominance <- function(N,s,h,howmany)
{
	i=0
	nfix = 0;
	while( i < howmany )
	{
		i=i+1
		t<-WFdominance(N,s,h)
		if ( t[length(t)] == 1 )
		{
			nfix=nfix+1
		}
	}
	return(nfix/howmany)
}

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
