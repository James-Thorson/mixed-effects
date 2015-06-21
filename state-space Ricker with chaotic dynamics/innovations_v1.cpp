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
  PARAMETER_VECTOR( e_t )
  
  // 
  int nobs = y_t.size();
  vector<Type> x_t(nobs);
  
  // Objective funcction
  Type jnll = 0;
  
  // Reconstruct time series
  x_t(0) = exp(log_k) * exp( e_t(0) );
  for( int t=1; t<nobs; t++){
    x_t(t) = x_t(t-1) * exp( exp(log_r) * (1-x_t(t-1)/exp(log_k))) * exp(e_t(t));
    jnll -= dnorm( e_t(t), Type(0), exp(log_sigmap), true );
  }
  
  // Probability of data
  for( int t=0; t<nobs; t++){
    jnll -= dnorm( log(y_t(t)), log(x_t(t)), exp(log_sigmam), true );
  }
  
  // Reporting
  ADREPORT( x_t );
  REPORT( x_t );

  return jnll;
}
