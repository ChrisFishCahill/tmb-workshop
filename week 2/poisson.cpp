/*    Generalized Linear Model for Gorbachev's Holy Goose  
*     y_i ~ Poisson(lambda_i)
*     E(y_i) = lambda_i
*     Var(y_i) = lambda_i 
*     log(lambda_i) = Beta0 + Beta1*x1 + Beta2*x2 
*     Cahill Oct 23 2019
*/
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );  //observation for critter i
  DATA_MATRIX( X_ij ); //design matrix

  // Parameters
  PARAMETER_VECTOR( b_j );
  vector<Type> eta_fixed_i = X_ij * b_j; //calculate the systematic component

  // Objective function
  int Nobs = y_i.size();
  Type nll = 0;
  
  for( int i=0; i<Nobs; i++ ){
    nll -= dpois( y_i(i), /*inverse link*/ exp(eta_fixed_i(i)), true ); //Random component
  }

  return nll;
}

