#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt );
  DATA_VECTOR( y_t );       
  
  // Parameters
  PARAMETER( x0 );
  PARAMETER( log_sigmaP );  
  PARAMETER( log_sigmaO );  
  PARAMETER( alpha ); //drift
  PARAMETER_VECTOR( x_t );
  
  // Objective funcction
  Type jnll = 0;
  
  //random effect
  jnll -= dnorm( x_t(0), x0, exp(log_sigmaP), true ); 
  
  for( int t=1; t<nt; t++){
    jnll -= dnorm( x_t(t), x_t(t-1) + alpha, exp(log_sigmaP), true );
  }
  
  // Probability of observations conditional on fixed and random effect values
  for( int t=0; t<nt; t++){
    jnll -= dnorm( y_t(t), x_t(t), exp(log_sigmaO), true );
  }
  
  // Reporting
  Type sigmaP = exp(log_sigmaP);
  Type sigmaO = exp(log_sigmaO);
  
  REPORT( sigmaP );
  REPORT( sigmaO );
  REPORT( x_t );
  
  ADREPORT( sigmaP );
  ADREPORT( sigmaO );
  ADREPORT( x_t );
  
  return jnll;
}
