HOC.update = function(pop,sigmas)
    {
        ##Phenotypic values.  Additive under HOC, ergo rowSums
        phenos = rowSums(pop)
        ##Fitnesses are Gaussian
        w = dnorm(phenos,0,sigmas)

        rv = matrix(data=0,ncol=ncol(pop),nrow=nrow(pop))
        parents = 1:length(w)
        for(i in 1:nrow(rv))
            {
                ##Sample this individual's parents
                parents.i = sample(parents,2,replace=TRUE,prob=w)
                ##Mendel
                rv[i,1] = pop[parents.i[1],(rbinom(1,1,0.5)+1)]
                rv[i,2] = pop[parents.i[2],(rbinom(1,1,0.5)+1)]
            }
        rv
    }

#' House-of-cards simulation
#' @export
HOC.sim = function(N,mu,sigmamu,sigmae,sigmas,ngens)
    {
        warning("Not clear that this simulation is correct...")
        ##columns are: allele 1, allele 2, random effect
        pop = matrix(data=0,ncol=3,nrow=N)
        for( i in 1:ngens )
            {
                m=rbinom(N,1,mu)
                m.effects = rnorm(length(which(m!=0)),0,sigmamu)
                pop[which(m!=0),1]=m.effects
                m=rbinom(N,1,mu)
                m.effects = rnorm(length(which(m!=0)),0,sigmamu)
                pop[which(m!=0),2]=m.effects
                pop[,3] = rnorm(N,0,sigmae)
                pop=HOC.update(pop,sigmas)
            }
        pop[,3] = rnorm(N,0,sigmae)
        pop
    }

HOC.h2 = function(pop)
    {
        return( list(VG=var(pop[,1]+pop[,2]),VE=var(pop[,3])) )
    }
