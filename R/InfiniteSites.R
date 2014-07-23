#the following functions are use to apply the
#infinitely-many sites model to a tree

#' Pick a random branch from a tree
#' @param tree A simualted coalescent history
#' @param n The sample size corresponding to tree
#' @param r A random number corresponding to runif(1,0,tree$ttot)
#' @note This is an internal function
#' @return The index of the randomly-chosen branch
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

#' Asks if an individual on a tree is a descendant of a branch
#' @param tree A simulated coalescent history
#' @param ind The integer label of a sample
#' @param branch The index of the branch
#' @return TRUE or FALSE
isdescendant=function(tree,ind,branch)
{
	i = ind
	while(i < branch)
	{
		i=tree$nodes[i]
	}
	return(i==branch)
}

#' Add specific number of mutations to a tree
#' @param tree A simulated coalescent history
#' @param n The sample size in tree
#' @param S The number of mutations to add
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

#' Infinitely-many sites mutation
#' @param tree A simulated coalescent history
#' @param n The sample size in tree
#' @param theta 4Nu
mutate=function( tree , n , theta )
{
	S = rpois(1,tree$ttot*theta/2)
	return( mutateS(tree,n,S) )
}
