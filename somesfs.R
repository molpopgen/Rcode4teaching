source("WFsnm.R")
source("WFisland.R")
source("WFbottle.R")
source("SummaryStats.R")

#calculate and simulate the unfolded sfs for n=20, theta = 10

N=20
THETA=10

exact = array(data=0,dim=N-1)
for(i in 1:(N-1))
{
	exact[i] = THETA/i
}


sim = array(data=0,dim=N-1)
struct = array(data=0,dim=N-1)
bottle = array(data=0,dim=N-1)
NREPS=1000
for(i in 1:NREPS)
{
	s = WFsample(N,THETA);
	s.sfs = sfs(s)
	sim = sim + s.sfs

	s=WFsampleMig(N/2,N-(N/2),0.05,THETA)
	s.sfs=sfs(s)
	struct=struct+s.sfs


	s=WFsampleBottle(N,0.008,2*0.019,0.03,THETA)
	s.sfs=sfs(s)
	bottle=bottle+s.sfs
}
sim = sim/NREPS
struct = struct/NREPS
bottle = bottle/NREPS

#normalize all the SFS
exact = exact/sum(exact)
sim = sim/sum(sim)
struct = struct/sum(struct)
bottle = bottle/sum(bottle)
postscript("somesfs.ps",height=10,width=10,pointsize=18)
barplot(rbind(exact,sim,struct,bottle),beside=T,legend.text=c("neutral","neutral sim","N1=N/2,4Nm=0.05/chromo","Bottleneck"),main="Simulations based on 1000 replicates")
dev.off()