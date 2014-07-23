## #load a file that contains functions we need
## #source("Coal.R")

## Dobs = -1.6

## p=0
## N=10000
## n=10
## S=11
## for(i in 1:N)
## {
## 	#generate a sample
## 	#under the WF model w/no recombination
## 	sim = WFsampleS(n,S)
## 	#calculate D for the sample
## 	sim.D = TajD(sim)
## 	if( sim.D <= Dobs )
## 	{
## 		p=p+1
## 	}
## }
## print(p/N)
