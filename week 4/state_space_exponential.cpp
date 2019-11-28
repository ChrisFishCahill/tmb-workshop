#include <TMB.hpp>
/* 
*y_obs_i ~  Normal(biomass_t, sigmaO)
*biomass_t = lambda(t-1)*biomass(t-1)
*lambda_t==0 = mu_lambda
*biomass_t==0 = B0 
*with the exception of b0, biomass_t is derived (or reconstructed if it helps)
*lambda_t is a latent state where lambda_t ~ N(lambda_t, lambda_mu, sigP)
*Cahill 26 November 2019 
*/ 
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( Nyears );   
  DATA_VECTOR( Y_obs_t );   //observed data
  
  // Parameters
  PARAMETER( logB0 );          //initial biomass measurement
  PARAMETER( log_sigmaP );     //process sd
  PARAMETER( log_sigmaO );     //measurement sd 
  PARAMETER( mu_lambda); 
  PARAMETER_VECTOR( lambda_t ); //population growth rate 

  // Objective function
  Type jnll = 0;

  vector<Type> biomass_t(Nyears); //create a vector to store biomass values
  biomass_t(0) = exp(logB0); //condition on B0
  // Probability of random coefficients--lambda_t --> vector of latent states
  //sweep downstream through time-series | B0, mu_lambda
  //jnll = dnorm(lambda_t(0), mu_lambda, exp(log_sigmaP), true); 
  
  for( int t=0; t<(Nyears-1); t++ ){
    jnll -= dnorm( lambda_t(t), mu_lambda, exp(log_sigmaP), true );
    //jnll -= dnorm( lambda_t(t), lambda_t(t-1), exp(log_sigmaP), true ); 
    biomass_t(t+1) = lambda_t(t)*biomass_t(t); 
  }
  
  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<Nyears; t++){
    if(!isNA(Y_obs_t(t))) jnll -= dnorm( Y_obs_t(t), biomass_t(t), exp(log_sigmaO), true ); 
  }
  
  // Reporting
  Type sigmaP = exp(log_sigmaP);
  Type sigmaO = exp(log_sigmaO);
  
  REPORT( sigmaP );
  REPORT( sigmaO );
  REPORT( lambda_t );
  REPORT( biomass_t); 
  
  ADREPORT( sigmaP );
  ADREPORT( sigmaO );
  ADREPORT( lambda_t );
  ADREPORT( biomass_t ); 
  
  return jnll;
}
