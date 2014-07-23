simWFtreeBottle = function(n,tr,f,d)
{
	if( tr < 0 || f < 0 || d < 0 )
	{
		stop("bad parameters")
	}

	nodes = array(data=-1,dim=2*n-1)
	times = array(data=0,dim=2*n-1)
	sample = array(1:n)
	nextnode = n+1
	t=0
	ttot=0

	k=n

	tcoal=-1
	while(k>1)
	{
		if( t < tr || t >= (tr+d) ) #pop. size is large
		{
			tcoal = rexp(1,k*(k-1)/2)
		} else #during bottleneck
		{
			tcoal = rexp(1,k*(k-1)/(2*f))
		}
		
		if( t < tr && (t+tcoal) >= tr )
		{
			temp = tr-t
			t=tr
			ttot=ttot+(k*temp)
		} else if( t<(tr+d) && (t+tcoal) >= (tr+d) )
		{
			temp = (tr+d)-t
			t=(tr+d)
			ttot=ttot+(k*temp)
		} else #coalescent event
		{
			t=t+tcoal
			ttot=ttot+(k*tcoal)
			two = pick2(k)
			nodes[sample[two$i]] = nextnode
			nodes[sample[two$j]] = nextnode
			times[nextnode]=t
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
	}
	treedata = list(times=times,nodes=nodes,tmrca=t,ttot=ttot)
	class(treedata) <- 'coalTree'
	return(treedata)
}

WFsampleBottle = function(n,tr,d,f,theta)
{
	t=simWFtreeBottle(n,tr,d,f)
	d=mutate(t,n,theta)
	s=list(times=t$times,nodes=t$nodes,tmrca=t$tmrca,ttot=t$ttot,pos=d$pos,types=d$types)
	class(s) = 'coalSample'
	return(s)
}
