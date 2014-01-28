
addLeftBranch = function(tree,ancestor,nextleft,nsam)
{
	for(i in 1:ancestor)
	{
		if(tree$nodes[i]==ancestor)
		{
			print(paste("left i = ",i))
			rheight=(tree$times[i]/tree$times[2*nsam-1])*(2*nsam-1)
			rheightAbv=(tree$times[tree$nodes[i]]/tree$times[2*nsam-1])*(2*nsam-1)
			segments(nsam/2-(nextleft)-1,rheight,nsam/2-(nextleft),rheightAbv)
			return(i)
			break;
		}
	}
	return(-1)
}

addRightBranch = function(tree,ancestor,nextright,nsam)
{
	for(i in ancestor:1)
	{
		if(tree$nodes[i]==ancestor)
		{
			print(paste("right i = ",i))
			rheight=(tree$times[i]/tree$times[2*nsam-1])*(2*nsam-1)
			rheightAbv=(tree$times[tree$nodes[i]]/tree$times[2*nsam-1])*(2*nsam-1)
			segments(nsam/2+(1)*nextright+1,rheight,nsam/2+(nextright),rheightAbv,col="red")
			return(i)
			break;
		}
	}
	return(-1)
}

plot.coalTree = function(tree,nsam)
{
	plot(c(0.5,nsam),c(0,2*nsam-1),pch="",xaxt="n",yaxt="n",xlab="",ylab="",main="An ugly tree!")
	xlabels=array(data=0,dim=nsam)
	nextleft = 0
	nextright= 0
	for(i in (2*nsam-1):1)
	{
		l=addLeftBranch(tree,i,0,nsam)
		while(l != -1)
		{
			l=addLeftBranch(tree,l,i-l,nsam)
		}
		r=addRightBranch(tree,i,0,nsam)
		while(r != -1)
		{
			r=addRightBranch(tree,r,i-r,nsam)
		}
	}
}
