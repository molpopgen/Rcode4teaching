#include <Rcpp.h>
#include <random>
#include <limits>
#include <vector>
#include <functional>

using namespace Rcpp;

NumericMatrix updatePop( const NumericMatrix & pop,
			 const std::vector<double> & fitnesses,
			 std::mt19937 & generator )
{
  RNGScope scope;
  NumericMatrix rv(pop.nrow(),pop.ncol());
  std::discrete_distribution<double> d(fitnesses.begin(),fitnesses.end());
  for( unsigned i = 0 ; i < pop.nrow() ; ++i)
    {
      size_t p1 = d(generator),p2=d(generator);
      int c1 = rbinom( 1, 1, 0.5 )[0],
	c2 = rbinom( 1, 1, 0.5 )[0];
      rv(i,0) = pop(p1, c1);
      rv(i,1) = pop(p2, c2);
    }
  return rv;
}

//Generic HOC simulation machine.
//The emodel parameter must be able to 
//take sigmamu (std. dev of effect sizes),
//and convert into a deviate from the desired distribution.
template<typename effectsizes>
NumericMatrix HOCsim_generic(const unsigned & N,
			     const double & mu,
			     const double & sigmamu,
			     const double & sigmae,
			     const double & sigmas,
			     const unsigned & ngens,
			     const effectsizes & emodel)
{
  NumericMatrix pop(N,3);
  std::vector<double> fitnesses(N);
  RNGScope scope;
  std::mt19937 generator;
  generator.seed( unsigned(runif(1,0,std::numeric_limits<unsigned>::max())[0]) );
  for( unsigned gen = 0 ; gen < ngens ; ++gen )
    {
      NumericVector mutants1 = rbinom(N,1,mu),mutants2=rbinom(N,1,mu);
      //Mutate and add Gaussian noise
      for( unsigned i = 0 ; i < N ; ++i )
	{
	  if( mutants1[i] )
	    {
	      pop(i,0) = emodel(sigmamu);
	    }
	  if( mutants2[i])
	    {
	      pop(i,1) = emodel(sigmamu);
	    }
	  pop(i,2) = as<double>(rnorm(1,0,sigmae));
	  fitnesses[i] = dnorm(NumericVector(1,sum( pop(i,_))),0.,sigmas)[0];
	}
      pop = updatePop(pop,fitnesses,generator); 
    }
  pop(_,2) = rnorm(N,0,sigmae);
  return pop;
}

//' Simulate a quantitative trait under a Gaussian House-of-Cards mutation model
//' @param N The number of diploids
//' @param mu The mutation rate (per gamete, per generation)
//' @param sigmamu The standard deviation in the effect size of new mutations.  Mutation effect sizes ~ N(0,sigmamu^2).
//' @param sigmae The standard deviation in random effects added to phenotype
//' @param sigmas The standard deviation in the Gaussuian fitness function
//' @param ngens The number of generations to simulate
//' @return A matrix.  Columns 1 and 2 are the alleles for each diploid.  Column 3 is the random effect.  Phenotype is therefore given by rowSums(result)
//' @details  Tersely: the mutation model is Kingman's house-of-cards with Gaussian effects.  The inspiration for this function is Turelli.
//' @references Kingman, J. (1978). A simple model for the balance between selection and mutation. Journal of Applied Probability, 15, 1-12.
//' @references Turelli, M. (1984). Heritable genetic variation via mutation-selection balance: Lerch's zeta meets the abdominal bristle. Theoretical Population Biology, 25(2), 138-193.
//' @export
//' @examples
//' \dontrun{
//' N=5e2
//' mu=0.001
//' sige=0.15
//' sigmu=0.5
//' sigs=1
//' pop = HOCsim(N,mu,sigmu,sige,sigs,8*N)
//' }
// [[Rcpp::export]] 
NumericMatrix HOCsim(const unsigned & N,
		     const double & mu,
		     const double & sigmamu,
		     const double & sigmae,
		     const double & sigmas,
		     const unsigned & ngens)
{
  return HOCsim_generic(N,mu,sigmamu,sigmae,sigmas,ngens, [](const double & x){ return rnorm(1,0.,x)[0]; } );
}

//' Simulate a quantitative trait under an exponential House-of-Cards mutation model
//' @param N The number of diploids
//' @param mu The mutation rate (per gamete, per generation)
//' @param sigmamu The standard deviation in the effect size of new mutations.  Mutation effect sizes are exponentially-distributed with rate lambda = 1/sigmamu.
//' @param sigmae The standard deviation in random effects added to phenotype
//' @param sigmas The standard deviation in the Gaussuian fitness function
//' @param ngens The number of generations to simulate
//' @return A matrix.  Columns 1 and 2 are the alleles for each diploid.  Column 3 is the random effect.  Phenotype is therefore given by rowSums(result)
//' @references Kingman, J. (1978). A simple model for the balance between selection and mutation. Journal of Applied Probability, 15, 1-12.
//' @export
//' @examples
//' \dontrun{
//' N=5e2
//' mu=0.001
//' sige=0.15
//' sigmu=0.5
//' sigs=1
//' pop = HOCsim_exp(N,mu,sigmu,sige,sigs,8*N)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix HOCsim_exp(const unsigned & N,
			 const double & mu,
			 const double & sigmamu,
			 const double & sigmae,
			 const double & sigmas,
			 const unsigned & ngens)
{
  return HOCsim_generic(N,mu,sigmamu,sigmae,sigmas,ngens, [](const double & x){ return rexp(1,x)[0]; } );
}

//Returns distribution of VG and VE
//Internal function
// [[Rcpp::export]] 
DataFrame VGVEdistHOC( const unsigned & N,
		       const double & mu,
		       const double & sigmamu,
		       const double & sigmae,
		       const double & sigmas,
		       const unsigned & ngens,
		       const unsigned & nreps, 
		       const bool & gaussian = true)
{
  NumericVector VG(nreps),VE(nreps);

  for(unsigned rep=0;rep<nreps;++rep)
    {
      NumericMatrix pop = (gaussian) ? HOCsim(N,mu,sigmamu,sigmae,sigmas,ngens) : HOCsim_exp(N,mu,sigmamu,sigmae,sigmas,ngens);
      NumericVector G = pop(_,0)+pop(_,1);
      VG[rep] = var(G);
      VE[rep] = var(pop(_,2));
    }
  return DataFrame::create( Named("VG")=VG,
			    Named("VE")=VE );
}
