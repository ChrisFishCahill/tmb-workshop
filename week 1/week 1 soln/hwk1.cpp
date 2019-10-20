#include <TMB.hpp>
/* Nonlinear regression example
*  Juveniles = (alpha*Adults) / (1 + beta*Adults) + c*temp + eps
*  where eps ~ N(0, SD)
*  
*  Cahill Oct 2019
*/
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( Y_i );  //measurement of juveniles i
  DATA_VECTOR( A_i );  //measurement of adults i
  DATA_VECTOR( T_i );  //measurement of temp i   
  
  // Parameters
  PARAMETER_VECTOR( parm_vec ); 
  
  Type alpha = parm_vec(0); 
  Type beta = parm_vec(1); 
  Type c_temp = parm_vec(2); 
  Type sigma = exp(parm_vec(3)); //sigma cannot be negative

  int Nobs = Y_i.size();
  vector<Type> Y_pred_i( Nobs ); 
  
  // Objective function
  Type nll = 0;
  
  for( int i=0; i<Nobs; i++ ){
    Y_pred_i(i) = (alpha*A_i(i)) /( (Type(1.0) + beta*A_i(i)) ) + c_temp*T_i(i); 
    nll -= dnorm( Y_i(i), Y_pred_i(i), sigma, true ); 
  }
  
  SIMULATE{
   Y_i = rnorm(Y_pred_i, sigma); 
   REPORT(Y_i); 
  }
  return nll;
}

