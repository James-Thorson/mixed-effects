// Space time 
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(n_data);         // Total number of observations
  DATA_VECTOR(Y);       	// Count data
  DATA_FACTOR(NAind);		// 1 = Y is NA, 0 = is not NA
  DATA_INTEGER(n_knots);
  DATA_INTEGER(n_stations)	// Number of stations 
  DATA_FACTOR(meshidxloc);	// Pointers into random effects vector x
  DATA_INTEGER(n_years)          // Number of years  
  DATA_INTEGER(n_p)          	// number of columns in covariate matrix X
  DATA_MATRIX(X);		// Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(alpha);   // Mean of Gompertz-drift field
  PARAMETER(phi);            // Offset of beginning from equilibrium
  PARAMETER(log_tau_E);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
  //PARAMETER_MATRIX(Epsilon_input);  // Spatial process variation
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in carrying capacity

  // objective function -- joint negative log-likelihood 
  using namespace density;
  Type jnll = 0;
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();
  
  // Spatial parameters
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Objects for derived values
  vector<Type> eta(n_data); 
  vector<Type> nu(n_data);
  vector<Type> mean_abundance(n_years);
  matrix<Type> log_Dji(n_stations,n_years);
  matrix<Type> Epsilon(n_knots,n_years);
  vector<Type> Omega(n_knots);
  vector<Type> Equil(n_knots);
 
  // Probability of Gaussian-Markov random fields (GMRFs)
  jnll_comp(0) += GMRF(Q)(Omega_input);
  jnll_comp(1) = SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);
  //jnll += SCALE(GMRF(Q),exp(-log_tau_E))(Omega_input);
  //for(int t=0; t<n_years; t++){ 
  //  if(t==0) jnll_comp(1) += GMRF(Q)(Epsilon_input.col(t));
  //  if(t>=1) jnll_comp(1) += GMRF(Q)(Epsilon_input.col(t)-rho*Epsilon_input.col(t-1));
  //}
  
  // Transform GMRFs
  eta = X*alpha.matrix();
  int ii = 0;
  for(int j=0; j<n_knots; j++){
    Omega(j) = Omega_input(j) / exp(log_tau_O);
    Equil(j) = eta(ii) + Omega(j) / (1-rho);
    ii++;
    for(int i=0; i<n_years; i++){ 
      Epsilon(j,i) = Epsilon_input(j,i) / exp(log_tau_E);
    }
  }
  
  // Likelihood contribution from observations
  ii = 0;
  for (int i=0;i<n_years;i++){
    mean_abundance(i) = 0;
    for (int j=0;j<n_stations;j++){ 
      nu[ii] = phi * pow(rho, i);
      log_Dji(j,i) = nu[ii] + Epsilon(meshidxloc(j),i) / exp(log_tau_O) + ( eta(ii) + Omega(meshidxloc(j)) )/(1-rho);
      mean_abundance[i] = mean_abundance[i] + exp( log_Dji(j,i) );      
      if(!NAind(ii)){                
        jnll_comp(2) -= dpois( Y[ii], exp( log_Dji(j,i) ), true );
      }
      ii++;      
    }
    mean_abundance[i] = mean_abundance[i] / n_stations;      
  }
  jnll = jnll_comp.sum();

  // Diagnostics
  REPORT( jnll_comp );
  REPORT( jnll );
  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  ADREPORT( Range );
  ADREPORT( SigmaE );
  ADREPORT( SigmaO );
  // Fields
  REPORT( log_Dji );
  REPORT( Epsilon );
  REPORT( Omega );
  REPORT( Equil );
  // Total abundance
  ADREPORT( log(mean_abundance) ); // standard errors in log-space
  ADREPORT( mean_abundance );      // standard errors in nominal-space
  
  return jnll;
}
