
sfs = function(d)
#calculate the unfolded site-frequency spectrum
{
	s = array(data=0,dim=nrow(d$types)-1)
	S = length(d$pos)
	if(S==0)
	{
		return (s)
	}
	for(i in 1:S)
	{
		nmuts = length(d$types[,i][d$types[,i]==1])
		s[nmuts] = s[nmuts]+1
	}
	return(s)
}

thetaw = function(d)
{
	nsam = nrow(d$types)
	S = length(d$pos)
	if(S==0){ return (0) }
	#factor of 2 here is b/c our an function above is in units of 2N generations
	return( 2*S/(an(nsam)) )
}

pi = function(d)
#calculate pi as sum of site heterozygosity
{
	if(length(d$pos)==0) { return(0) }
	pq=0
	nsam = nrow(d$types)
	for(i in 1:ncol(d$types))
	{
		nmuts = length(d$types[,i][d$types[,i]==1])
		pq = pq + nmuts*(nsam-nmuts)
	}
	return ((2*pq)/(nsam*(nsam-1)))
}

thetah = function(d)
#Fay and Wu's ThetaH
{
	if(length(d$pos)==0) { return(0) }
	H=0
	for(i in 1:ncol(d$types))
	{
		nc = length(d$types[,i][d$types[,i]==1])
		H = H + (nc*nc)
	}
	nsam = nrow(d$types)
	return (2*H/(nsam*(nsam-1)))
}

#Tajima's D requires a lot of functions
bn = function(n)
{
	b=0
	for(i in 1:n)
	{
		b = b+1/(i^2)
	}
	return (b)
}

TajdDdenominator = function(S,n)
{
	a = an(n)/2;
	b = bn(n)

	b1 = (n+1)/(3*(n-1))
	b2 = (2*(n^2) + n + 3)/(9*n*(n-1))
	c1 = b1-1/a
	c2 = b2 - (n+2)/(a*n) + b/(a^2)
	e1 = c1/a
	e2 = c1/(a^2 + b)
	return ( sqrt(e1*S + e2*S*(S-1)) )
}

TajD = function(d)
{
	if( length(d$pos) == 0 )
	{
		#return an undefined value if there are no segregating sites
		return(0/0)
	}
	return ( (pi(d)-thetaw(d))/TajdDdenominator(length(d$pos),nrow(d$types)) )
}

Nhaps = function(d)
#calculate the number of haplotypes (unique sequences) in the sample
{
	if( length(d$pos) == 0 ){ return(1) }
	return( nrow(unique(as.data.frame(d$types))) )
}