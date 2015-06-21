// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_t );                
  
  // Parameters
  PARAMETER( log_r );
  PARAMETER( log_k );
  PARAMETER( log_sigmap );
  PARAMETER( log_sigmam );
  PARAMETER_VECTOR( log_x_t )
  
  // 
  int nobs = y_t.size();
  vector<Type> xpred_t(nobs);
  
  // Objective funcction
  Type jnll = 0;
  
  // Reconstruct time series
  for( int t=1; t<nobs; t++){
    xpred_t(t) = exp(log_x_t(t-1)) * exp( exp(log_r) * (1-exp(log_x_t(t-1))/exp(log_k)));
    jnll -= dnorm( log_x_t(t), log(xpred_t(t)), exp(log_sigmap), true );
  }
  
  // Probability of data
  for( int t=0; t<nobs; t++){
    jnll -= dnorm( log(y_t(t)), log_x_t(t), exp(log_sigmam), true );
  }
  
  // Reporting
  ADREPORT( log_x_t );
  REPORT( log_x_t );

  return jnll;
}
