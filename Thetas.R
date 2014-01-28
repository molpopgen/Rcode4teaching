source("Coal.R")
#calculate the distribution of 3 estimators of theta
#for fun, let's do it under a bottleneck
watterson = array()
tajima = array()
faywu = array()
n=20
theta=10
tr=0.008
d=0.03
f=0.04
N=5000
for(i in 1:N)
{
	s=WFsampleBottle(n,tr,d,f,theta)
	watterson[i] = thetaw(s)
	tajima[i] = pi(s)
	faywu[i] = thetah(s)
}

plot(density(log(watterson)),type="l",xlab=expression(paste("log(",hat(theta),")")),main="")
abline(v=log(mean(watterson)))
lines(density(log(tajima)),col="red")
abline(v=log(mean(tajima)),col="red")
lines(density(log(faywu)),col="blue")
abline(v=log(mean(faywu)),col="blue")
abline(v=log(10),lty="dashed",lwd=3)