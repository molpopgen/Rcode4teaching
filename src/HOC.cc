#include <Rcpp.h>
#include <random>
#include <limits>
#include <vector>
using namespace Rcpp;

NumericMatrix updatePop( const NumericMatrix & pop,
			 const std::vector<double> & fitnesses,
			 std::mt19937 & generator )
{
  RNGScope scope;
  NumericMatrix rv(pop.nrow(),pop.ncol());
  std::discrete_distribution<double> d(fitnesses.begin(),fitnesses.end());
  //stop("No more crashes!");
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

//' @export
// [[Rcpp::export]] 
NumericMatrix HOCsim(const unsigned & N,
		     const double & mu,
		     const double & sigmamu,
		     const double & sigmae,
		     const double & sigmas,
		     const unsigned & ngens)
{
  NumericMatrix pop(N,3);
  //std::vector< std::vector<double> > pop(3,std::vector<double
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
	      pop(i,0) = as<double>(rnorm(1,0,sigmamu));
	    }
	  if( mutants2[i])
	    {
	      pop(i,1) = as<double>(rnorm(1,0,sigmamu));
	    }
	  pop(i,2) = as<double>(rnorm(1,0,sigmae));
	  NumericVector pheno(1,pop(i,0)+pop(i,1)+pop(i,2));
	  fitnesses[i] = dnorm( pheno,0.,sigmas )[0];
	}
      pop = updatePop(pop,fitnesses,generator); 
    }
  pop(_,2) = rnorm(N,0,sigmae);
  return pop;
}
