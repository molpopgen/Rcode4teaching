#the following functions are use to apply the
#infinitely-many sites model to a tree

#'pick a random branch from a tree
pickbranch=function(tree,n,r)
{	
	t = 0
	for(i in 1:(2*n-1))
	{
		t = t+(tree$times[tree$nodes[i]]-tree$times[i])
		if(t >= r)
		{
			return(i)
		}
	}
	return(i)
}

isdescendant=function(tree,ind,branch)
{
	i = ind
	while(i < branch)
	{
		i=tree$nodes[i]
	}
	return(i==branch)
}

mutateS = function( tree, n , S )
{
	positions=array(0.,dim=S);
	genotypes=matrix(data=0,nrow=n,ncol=S)
	if(S>0)
	{
	for(i in 1:S)
	{
		positions[i] = runif(1,0.,1.)
		branch = pickbranch(tree,n,runif(1,0.,tree$ttot))
		for(ind in 1:n)
		{
			if(isdescendant(tree,ind,branch)==TRUE)
			{
				genotypes[ind,i]=1
			}
		}
	}
	}
	return(list(pos=sort(positions),types=genotypes))
}

mutate=function( tree , n , theta )
{
	S = rpois(1,tree$ttot*theta/2)
	return( mutateS(tree,n,S) )
}
