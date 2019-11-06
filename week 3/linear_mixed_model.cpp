#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( Ngroups );  //number of groups
  DATA_IVECTOR( g_i );      //group indicator
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 );
  PARAMETER( log_SD );      //likelihood noise term
  PARAMETER( log_SDG );     //random effect
  PARAMETER_VECTOR( eps_g );//vector for random effects
  
  // Objective funcction
  Type jnll = 0;
  int n_i = y_i.size();
  
  // Probability of random coefficients
  for( int g=0; g<Ngroups; g++){
    jnll -= dnorm( eps_g(g), Type(0.0), exp(log_SDG), true );
  }
  
  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
    jnll -= dnorm( y_i(i), beta0 + eps_g(g_i(i)), exp(log_SD), true );
  }
  
  // Reporting
  Type SDG = exp(log_SDG);
  Type SD = exp(log_SD);
  
  REPORT( SDG );
  REPORT( SD );
  ADREPORT( SDG );
  ADREPORT( SD );
  
  return jnll;
}
