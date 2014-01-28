#the functions below are all for simulating trees under the WF coalescent

swap = function(i,j)
#return a list which swaps the values of i and j
{
	temp = j
	j=i
	i=temp
	return(list(i=i,j=j))
}

pick2 = function(n)
#randomly pick two integers, i and j,
#such that 1 \leq i \leq n and 1 \leq j \leq n and i < j
{
	i = as.integer(runif(1,0,n))+1
	j = as.integer(runif(1,0,n))+1
	while(i==j)
	{
		j = as.integer(runif(1,0,n))+1
	}
	if(j < i)
	{
		return(swap(i,j))
	}
	return(list(i=i,j=j))
}

pick2indeme = function(n1,n2,deme,config)
#randomly pick two integers, i and j, such that they are in the same deme, and
#such that 1 \leq i \leq n and 1 \leq j \leq n and i < j
{
	two = list()
	if( deme == 1)
	{
		two = pick2(n1)
	} else
	{
		two = pick2(n2)
	}
	c1=0
	c2=0
	i=1
	j=1
	while( i <= (n1+n2) )
	{
		if(config[i]==deme && j == two$i)
		{
			c1=i
		}
		else if(config[i]==deme && j == two$j)
		{	
			c2=i
		}
		if(config[i]==deme)
		{
			j=j+1
		}
		i=i+1
	}
	return(list(i=min(c1,c2),j=max(c1,c2)))
}

totaltime = function( tree, n )
{
	tt = 2*tree$times[2*n-1]
	for( i in 1:(2*n-2) )
	{
		tt = (tt+tree$times[i])
	}
	return(tt)
}


