
// von Bertalanffy growth function
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data
DATA_VECTOR(Age);
DATA_VECTOR(Length);
int n = Age.size(); // get number of data points to loop over

// parameters:
PARAMETER(Linf); // asymptoptic length
PARAMETER(K); // growth rate
PARAMETER(T); // Age at length 0
PARAMETER(LogSigma); // log(residual SD)
// Fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type sigma = exp(LogSigma);
vector<Type> LengthPred(n);

// Report transformed variables
ADREPORT(sigma);

Type nll = 0.0; // Initialize negative log-likelihood

for(int i = 0; i < n; i++){ // C++ starts loops at 0!
// get negative log likelihood (last argument is log = TRUE)
LengthPred(i) = Linf * (1.0-exp(-K*(Age(i) - T))) ;
}

nll = -sum(dnorm(Length,LengthPred , sigma,true)); // Calculate the negative log-likelihood
// NOTE: the negative log-likelihood can also be included in the for loop using the subtraction assignment '-=' where:
// nll -= dnorm(Length,LengthPred , sigma,true)
// This will subtract the each log-likelihood from the initial 0 value
return nll;
}
