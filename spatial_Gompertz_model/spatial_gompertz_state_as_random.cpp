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
  PARAMETER(log_tau_U);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)

  // Random effects
  PARAMETER_ARRAY(log_Dji);  // Spatial process variation
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
  Type SigmaU = 1 / sqrt(4*pi*exp(2*log_tau_U)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Objects for derived values
  vector<Type> eta(n_data); 
  vector<Type> nu(n_data);
  vector<Type> mean_abundance(n_years);
  array<Type> log_Dji_hat(n_knots,n_years);
  vector<Type> Omega(n_knots);
  vector<Type> Equil(n_knots);
 
  // Transform GMRFs
  eta = X*alpha.matrix();
  int ii = 0;
  for(int j=0; j<n_knots; j++){
    //Omega(j) = Omega_input(j) * exp(-log_tau_O);
    Omega(j) = Omega_input(j);
    Equil(j) = ( eta(ii) + Omega(j) ) / (1-rho);
    ii++;
  }
  
  // Calculate expectation for state-vector
  ii = 0;
  for(int j=0; j<n_knots; j++){
  for(int i=0; i<n_years; i++){ 
    if(i==0) log_Dji_hat(j,i) = phi + Equil(j);
    if(i>=1) log_Dji_hat(j,i) = rho*log_Dji(j,i-1) + eta(ii) + Omega(j);
    ii++;      
  }}
  
  // Probability of Gaussian-Markov random fields (GMRFs)
  // Omega
  // Default -- DOES WORK
    //jnll_comp(0) = GMRF(Q)(Omega_input);
  // SCALE function -- DOES WORK
    // If doing it this way, Omega=Omega_input (i.e., Omega_input not rescaled prior to being used in the likelihood)
    jnll_comp(0) = SCALE(GMRF(Q), exp(-log_tau_O))(Omega);    
  // Simple rescaling -- DOES WORK
    //jnll_comp(0) = GMRF(Q)( Omega * exp(log_tau_O));
  // Experimental -- DOESN'T WORK!
    //Eigen::SparseMatrix<Type> Qtmp = Q * exp(2*log_tau_O);
    //jnll_comp(0) = GMRF(Qtmp)(Omega);
  
  // Epsilon
  //jnll += SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);
  for(int i=0; i<n_years; i++){ 
    // SCALE(f,b) calculates density g(y) = f(x), where f is a known density and y=x*b
    jnll_comp(1) += SCALE( GMRF(Q), exp(-log_tau_U) )(log_Dji.col(i)-log_Dji_hat.col(i));
    //jnll_comp(1) += GMRF(Q)( (log_Dji.col(i)-log_Dji_hat.col(i)) * exp(-log_tau_U));       // DOESN'T WORK!
  }
  
  // Likelihood contribution from observations
  mean_abundance.setZero();
  ii = 0;
  for (int i=0;i<n_years;i++){
  for (int j=0;j<n_stations;j++){ 
    mean_abundance(i) += exp( log_Dji(j,i) ) / n_stations;      
    if( !NAind(ii) ){                
      jnll_comp(2) -= dpois( Y(ii), exp( log_Dji(meshidxloc(j),i) ), true );
    }
    ii++;
  }}
  jnll = jnll_comp.sum();

  // Diagnostics
  REPORT( jnll_comp );
  REPORT( jnll );
  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaU );
  REPORT( SigmaO );
  ADREPORT( Range );
  ADREPORT( SigmaU );
  ADREPORT( SigmaO );
  // Fields
  REPORT( log_Dji );
  REPORT( log_Dji_hat );
  REPORT( Omega );
  REPORT( Equil );
  // Total abundance
  ADREPORT( log(mean_abundance) ); // standard errors in log-space
  ADREPORT( mean_abundance );      // standard errors in nominal-space
  
  return jnll;
}
