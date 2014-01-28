WFtimes = function(n)
#Simulate the waiting times for a sample of size n 
#from the WF coalescent.  
#Returns a list containing tmrca and ttot, 
#the total time on the tree
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

an = function(n)
#The expected value of the total time on the tree.
#Warning: this function returns it on the scale of 2N generations,
#so to use it as the denominator for Watterson's theta, etc.
#divide the return value by 2
{
	tt=0
	for(i in 1:(n-1))
	{
		tt=tt+1/i
	}
	return (2*tt)
}