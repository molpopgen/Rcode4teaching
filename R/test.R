## #source("Coal.R")

## n=20
## S=10
## theta=10
## nreps=1000
## p=0
## for(i in 1:nreps)
## {
## 	s=WFsample(n,theta)
## 	if( length(s$pos) <= S)
## 	{
## 		p=p+1
## 	}
## }
## print(p/nreps)
