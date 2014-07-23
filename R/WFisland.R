##source("WFcoal.R")
##source("InfiniteSites.R")

simWFtreeMig = function(n1,n2,migrate)
#WF Island model w/symmetric migration
#migrate is per chromosome in the sample
{
	NSAM=n1+n2
	nodes = array(data=-1,dim=2*NSAM-1)	
	times = array(data=0,dim=2*NSAM-1)
	sample = array(1:NSAM)
	config = array(data=1,dim=NSAM)
	for(i in (n1+1):NSAM)
	{
		config[i]=2
	}

	nextnode=NSAM+1
	t=0
	ttot=0
	k=NSAM
	while(k>1)
	{
		if (k != n1+n2)
		{
			stop("sample sizes are funny")
		}
		#rates of events
		rcoal1 = n1*(n1-1)/2	#deme 1
		rcoal2 = n2*(n2-1)/2	#deme 2
		rmig = migrate*(n1+n2)

		tc1 = 0
		if(rcoal1 > 0)
		{
			tc1 = rexp(1,rcoal1)
		} else
		{
			tc1 = 9999999
		}
		tc2=0
		if(rcoal2 > 0)
		{
			tc2 = rexp(1,rcoal2)
		} else
		{
			tc2 = 9999999
		}
		tmig=0
		if(rmig > 0)
		{
			tmig = rexp(1,rmig)
		} else
		{
			tmig=999999
		}

		tmin = min(tc1,tc2,tmig)
		t = t+tmin
		ttot = ttot+k*tmin
		#print(paste(n1,n2,tmin,tc1,tc2,tmig))
		if( tmin == tc1 ) #coalescence in deme 1
		{
			two=pick2indeme(n1,n2,1,config)
			if(config[two$i] != 1 || config[two$j] != 1)
			{
				stop("chromosomes chosen incorrectly in CA1")
			}
			nodes[sample[two$i]] = nextnode
			nodes[sample[two$j]] = nextnode
			times[nextnode]=t
			sample[two$i] = nextnode
			nextnode=nextnode+1
			temp=array(dim=k-1)
			tempc=array(dim=k-1)
			dummy=1
			for(i in 1:k)
			{
				if(i != two$j)
				{
					temp[dummy]=sample[i]
					tempc[dummy]=config[i]
					dummy=dummy+1
				}
			}	
			sample=temp
			config=tempc
			n1=n1-1
			k=k-1
			if( length(config[config==1]) != n1 )
			{
				stop("error in CA1 at time ",t)
			}

			if( length(config) != k  )
			{
				stop("error2 in CA1 at time ",t)
			}
		} else if (tmin == tc2) #coalescence in deme 2
		{
			two=pick2indeme(n1,n2,2,config)
			if(config[two$i] != 2 || config[two$j] != 2)
			{
				stop("chromosomes chosen incorrectly in CA2")
			}
			nodes[sample[two$i]] = nextnode
			nodes[sample[two$j]] = nextnode
			times[nextnode]=t
			sample[two$i] = nextnode
			nextnode=nextnode+1
			temp=array(dim=k-1)
			tempc=array(dim=k-1)
			dummy=1		
			for(i in 1:k)
			{
				if(i != two$j)
				{
					temp[dummy]=sample[i]
					tempc[dummy]=config[i]
					dummy=dummy+1
				}
			}	
			sample=temp
			config=tempc
			n2=n2-1
			k=k-1
			if( length(config[config==2]) != n2 )
			{
				stop("error in CA2 at time ",t)
			}

			if( length(config) != k  )
			{
				stop("error2 in CA2 at time ",t)
			}
			
		} else #migration event
		{
			if( runif(1,0,1) < n1/(n1+n2) ) # mig from 1 to 2
			{
				#pick a chromosome from 1 to n1, and migrate it
				ch = as.integer(runif(1,0,n1))+1
				dummy=0
				for( temp in 1:length(config) )
				{
					if(config[temp]==1)
					{
						dummy=dummy+1
					}
					if(dummy==ch)
					{
						break
					}
				}
#				print(config)
				config[temp]=2
#				print(config)
				n1=n1-1
				n2=n2+1
				if( length(config[config==1]) != n1 )
				{
					stop("error in mig1")
				}
			} else #mig from 2 to 1
			{
				#pick a chromosome from 1 to n2, and migrate it
				ch = as.integer(runif(1,0,n2))+1
				dummy=0
				for( temp in 1:length(config) )
				{
					if(config[temp]==2)
					{
						dummy=dummy+1
					}
					if(dummy==ch)
					{
						break
					}
				}
				config[temp]=1
				n2=n2-1
				n1=n1+1
				if( length(config[config==2]) != n2 )
				{
					stop("error in mig2")
				}
			}
		}	
	}
	treedata = list(times=times,nodes=nodes,tmrca=t,ttot=ttot)
	class(treedata) <- 'coalTree'
	return(treedata)
}

WFsampleMig = function(n1,n2,migrate,theta)
{
	t=simWFtreeMig(n1,n2,migrate)
	d=mutate(t,(n1+n2),theta)
	s=list(times=t$times,nodes=t$nodes,tmrca=t$tmrca,ttot=t$ttot,pos=d$pos,types=d$types)
	class(s) = 'coalSample'
	return(s)
}
