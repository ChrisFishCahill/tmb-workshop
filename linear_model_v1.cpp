// Linear regression
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( Y_i );  //observation for critter i
  DATA_VECTOR( x1_i ); //covariate critter i 
  
  // Parameters
  PARAMETER( Beta0 ); 
  PARAMETER( Beta1 ); 
  PARAMETER( logSD ); 
  
  Type SD = exp( logSD );
  int Nobs = Y_i.size();
  vector<Type> Y_pred_i( Nobs ); 
  
  // Objective function
  Type jnll = 0;
  
  for( int i=0; i<Nobs; i++){
    Y_pred_i(i) = Beta0 + Beta1*x1_i(i); 
    jnll -= dnorm( Y_i(i), Y_pred_i(i), SD, true ); //works like dnorm in R
  }
  
  // Reporting
  ADREPORT(Y_pred_i); //SEs -- can be slow for big models 
  //for predictions--see also REPORT()
  return jnll;
}
