##update the population frequencies
updateP <- function(data)
    {
        P <- array(1,c(2,length(data)-1));
        n1 <- length(data[,1][data[,1]==1]);
        n2 <- length(data[,1][data[,1]==2]);
        for (i in 2:length(data) )
            {
                x1 <- length(data[,i][data[,i]==1 & data[,1]==1]);# num of allele 1 in pop1
                x2 <- length(data[,i][data[,i]==1 & data[,1]==2]);
                
       P[1,i-1] <- rbeta(1,1+x1,1+n1-x1);
                P[2,i-1] <- rbeta(1,1+x2,1+n2-x2);
            }
        return(P);
    }

##update population assignments of individuals
updateZ <- function(data,P,N)
    {
        for (ind in 1:N)
            {
                ##eqn A8 from Pritchard et al. is pgeno/(denom1+denom2)
                ##for the case of 2 pops
                pgeno <- 1;
                denom1 <- 1; #component of denominator for A8 from pop1
                denom2 <- 1; #same, from pop2
                pop <- data[2*ind-1,1];
                for(m in 2:length(data))
                    {
                        A <- data[2*ind-1,m];
                        a <- data[2*ind,m];
                        if ( A == 1 )
                            {
                                pgeno <- pgeno * P[pop,m-1];
                                denom1 <- denom1 * P[1,m-1];
                                denom2 <- denom2 * P[2,m-1];
                            }
                        else
                            {	
                                pgeno <- pgeno * (1-P[pop,m-1]);
                                denom1 <- denom1 * (1-P[1,m-1]);
                                denom2 <- denom2 * (1-P[2,m-1]);
                            }
                        if ( a == 1 )
                            {
                                pgeno <- pgeno * P[pop,m-1];
                                denom1 <- denom1 * P[1,m-1];
                                denom2 <- denom2 * P[2,m-1];
                            }
                        else
                            {	
                                pgeno <- pgeno * (1-P[pop,m-1]);
                                denom1 <- denom1 * (1-P[1,m-1]);
                                denom2 <- denom2 * (1-P[2,m-1]);
                            }
                        pgeno <- pgeno * 2; #because it's a 2pq thing...
                        denom1 <- denom1 * 2;
                        denom2 <- denom2 * 2;
                    }
                pZ <- pgeno/(denom1+denom2);
                u <- runif(1);
                if ( u <= pZ )
                    {
                        ##do nothing...we accept this pop assignment
                    }
                else
                    {
                        if(pop == 1)
                            {
                                data[2*ind-1,1] <- 2;
                                data[2*ind,1] <- 2;
                            }
                        else
                            {
                                data[2*ind-1,1] <- 1;
                                data[2*ind,1] <- 1;
                            }
                    }
            }
        return(data);
    }

#' The "structure" algorithm for biallelic data and K = 2
#' @param data A data frame/matrix.  The first column must be sampling locations, and the remaining columns are genotypes
#' @param burn Number of MCMC iterations for burn in.
#' @param m Number of MCMC generations after the burn in step.
#' @return An array of assignment probabilities.  These are the probability that each individual belongs to population 1.
#' @note This is the most basic implementation of the structure algorithm possible.  The implementation here
#' was in response to homework assignment number 4 in Carlos Bustamante's Biometry 694 course at Cornell.
#' @details See help(genoData) for more detail on the format that data must take.
#' @references Pritchard, J. K., M. Stephens, and P. Donnelly (2000) Inference of Population Structure Using Multilocus Genotype Data. Genetics 155:945-959
#' @examples
#' data(genoData)
#' genoData.assigned = bistruct(genoData)
#' @export
bistruct = function(data,burn=2000,m=2000)
    {
        N = nrow(data)/2
        pops=unique(data[,1])
        randomPops = sample(pops,N,replace=T)
        for( i in 1:N )
            {
                data[2*i-1,1] = randomPops[i]
                data[2*i,1] = randomPops[i]
            }
        assignments = array(0,10)
         for( i in 1:burn)
             {
                 P <- updateP(data);
                 data <- updateZ(data,P,N);
             }
        for (i in 1:m)
            {
                ##step 1 of the algorithm
                P <- updateP(data);
                ##step 2
                data <- updateZ(data,P,N);
                ##in assignments array, keep track
                ##of num runs for which individual j
                ## is in pop 1
                for (j in 1:N)
                    {
                        if (data[2*j-1,1] == 1)
                            {
                                assignments[j] <- assignments[j] + 1;
                            }
                    }
            }
        return(assignments/m)
    }
