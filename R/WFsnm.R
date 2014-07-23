#source("WFcoal.R")
#source("InfiniteSites.R")

#' Simulate a genealogy under the equilibrium WF model for a sample of size n.
#' @param n The sample size
#' @return A list containing:
#' nodes: An array of 2n-1 integers.  nodes[i] is the ancestor of the i-th node.  
#'	 The root is labeled with a -1
#' times: An array of 2n-1 times, corresponding to the values in nodes.
#' tmrca: The tmrca of the sample, equivalent to times[2n-1]
#' ttot:  The total time on the tree
simWFtree = function(n)

{
	nodes = array(-1,dim=2*n-1)
	times = array(0.,dim=2*n-1)
	sample = array(1:n)
	nextnode = n+1
	t=0
	ttot=0
	k=n
	while( k > 1 )
	{
		tcoal = rexp(1,k*(k-1)/2)
		t = t+tcoal
		ttot = ttot + k*tcoal
		two = pick2(k)
		nodes[sample[two$i]] = nextnode
		nodes[sample[two$j]] = nextnode
		times[nextnode]=t

		#we now need to update the sample array
		#We do this in the following steps:
		# 1. replace the i-th individual with the ancestor (nextnode)
		# 2. remove the j-th individual from the sample
		#Since R has no "pass-by-reference" capability, we do it by copying
		sample[two$i] = nextnode
		nextnode=nextnode+1
		temp=array(dim=k-1)
		dummy=1
		for(i in 1:k)
		{
			if(i != two$j)
			{
				temp[dummy]=sample[i]
				dummy=dummy+1
			}
		}
		sample=temp
		k=k-1
	}
	treedata = list(times=times,nodes=nodes,tmrca=t,ttot=ttot)
	class(treedata) <- 'coalTree'
	return(treedata)
}

#' Simulate a sample of size n under the idealized Wright-Fisher model
#' @param n The sample size
#' @param theta The mutation rate (4Nu)
#' @return An object of class coalSample
WFsample = function(n,theta)
{
	#generate a tree
	t=simWFtree(n)
	#apply mutations
	d=mutate(t,n,theta)
	s=list(times=t$times,nodes=t$nodes,tmrca=t$tmrca,ttot=t$ttot,pos=d$pos,types=d$types)
	class(s) = 'coalSample'
	return(s)
}

#' Simulate a sample of size n under the idealized Wright-Fisher model
#' @param n The sample size
#' @param S The number of mutations to place in the sample
#' @return An object of class coalSample
WFsampleS = function(n,S)
{
	#generate a tree
	t=simWFtree(n)
	#apply mutations
	d=mutateS(t,n,S)
	s=list(times=t$times,nodes=t$nodes,tmrca=t$tmrca,ttot=t$ttot,pos=d$pos,types=d$types)
	class(s) = 'coalSample'
	return(s)
}
