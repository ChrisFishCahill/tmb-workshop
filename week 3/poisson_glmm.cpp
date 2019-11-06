#include <TMB.hpp>
//  y_s_i     ~  Poisson(exp(eta_s_i))
//  eta_s_i<- exp(beta0 + delta_i + eps_s)
//  delta_i ~  N(0, sig_i)  //obs dispersion from site-level prediction
//  eps_s   ~  N(0, sig_s)  //site deviations away from global mu
//  Cahill 6 November 2019 

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( Nobs );
  DATA_INTEGER( Nsite );
  DATA_IVECTOR( s_i );
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 ); 
  PARAMETER( log_SD_s ); //site 
  PARAMETER( log_SD_i ); //obs
  
  //Random Effects
  PARAMETER_VECTOR( eps_s );  //site
  PARAMETER_VECTOR( delta_i );//nobs
    
  // Objective function
  Type jnll = 0;
  
  // Probability of random coefficients for site level eps_s
  // Among-site data-generation process
  for(int s=0; s<Nsite; s++){
    jnll -= dnorm( eps_s(s), Type(0.0), exp(log_SD_s), true );
  }
  
  // Probability of random coefficients for obs level delta_i
  // Observation-level data generation process
  for(int i=0; i<Nobs; i++){
    jnll -= dnorm( delta_i(i), Type(0.0), exp(log_SD_i), true );
  }
  
  // Probability of data conditional on fixed and random effect values
  // Ecological process of interest 
  vector<Type> eta_s_i(Nobs);
  for(int i=0; i<Nobs; i++){
    eta_s_i(i) =  beta0  + eps_s( s_i(i) ) + delta_i(i); //linear predictor
    if( !isNA(y_i(i)) ) jnll -= dpois( y_i(i), exp(eta_s_i(i)), true );
  }

  // Reporting
  Type SD_s = exp(log_SD_s);
  Type SD_i = exp(log_SD_i); 

  ADREPORT( SD_s );
  ADREPORT( SD_i ); 

  return jnll;
}
