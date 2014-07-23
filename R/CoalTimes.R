#' Waiting times from W-F coalescent
#' @param n The sample size
#' @return A list with TMRCA and the total time on the tree based on 1 simulated replicate. Units are 4N generations.
WFtimes = function(n)
{
	t=0
	tt=0
	k=n
	while(k>1)
	{
		ttemp = rexp(1,k*(k-1)/2)
		t = t + ttemp
		tt = tt + k*ttemp
		k = k - 1
	}
	return(list(tmrca=t,ttot=tt))
}

#' Expected total time on genealogy for W-F coalescent
#' @param n The sample size.
#' @return The expected total time on the tree (in units of 4N generations)
an = function(n)
{
	tt=0
	for(i in 1:(n-1))
	{
		tt=tt+1/i
	}
	return (2*tt)
}
