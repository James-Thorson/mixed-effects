////////////////////////////////////////////////////////
// ADMB-RE code for a separable state-space Schaefer surplus production model
//
// This model should be compiled with -r flag for random effects,
// I.e. admb -r admbmodel
//
// Code modified from original by Casper W. Berg, 13-07-2010
// Available here: https://groups.nceas.ucsb.edu/non-linear-modeling/projects
// Published: Bolker et al. (2013) Methods in Ecology and Evolution, 4: 501-512.
//
// Data compiled in Meyer and Millar (1999) CJFAS, 56: 1078-1087.
//
// Code modified by James T. Thorson, 2014-07-29
// 
////////////////////////////////////////////////////////
DATA_SECTION
  !!USER_CODE ad_comm::change_datafile_name("tuna.dat");
  init_number N;
  init_vector Catch(1,N);
  init_vector BiomassIndex(1,N);
PARAMETER_SECTION
  init_bounded_number logSdU(-4.6,0);
  init_bounded_number logr0(-4.6,0.183);  //# Bounds: log(0.01) - log(1.2) from Meyer and Millar 1999
  init_bounded_number logQ(-4.6,4.6);
  init_bounded_number K(10.0,2000.0);
  init_bounded_number logU1(0,6.9);
  random_effects_bounded_vector logU(2,N,0,6.9);
  objective_function_value jnll;
  
  // derived quantities
  number logSdy;

  // reporting quantities
  sdreport_number r0;
  sdreport_number SdU;
  sdreport_number Sdy;
  sdreport_number Q;
  sdreport_vector U(1,N);
  sdreport_vector qU(1,N);
  
PRELIMINARY_CALCS_SECTION
  logSdU = log(0.1);
  logr0 = log(0.2);
  logQ = log(1);
  K = 500;
  logU(1) = logU1;
  for(int i=2; i<=N; i++){
    logU(i) = log(80);
  }

PROCEDURE_SECTION
  r0 = mfexp(logr0);
  logSdy = logSdU;
  SdU = mfexp(logSdU);
  Sdy = mfexp(logSdy);
  Q = mfexp(logQ);
  U(1) = mfexp(logU1);
  qU(1) = Q * mfexp( logU1 );
  for(int i=2; i<=N; ++i){
    qU(i) = Q*U(i);
    U(i) = mfexp(logU(i));
  }
  
  jnll = 0.0; // joint negative log-likelihood to be minimized.

  // transition equation (Eq. 3)
  step(logU1,logU(2),logSdU,logr0,K,2);
  for(int i=3; i<=N; i++){
    step(logU(i-1),logU(i),logSdU,logr0,K,i-1);
  }
  // observation equation (Eq. 4)
  obs( logQ, logSdU, logU1, 1);
  for(int i=2; i<=N; i++){
    obs( logQ, logSdU, logU(i), i);
  }

  prior( logr0 );

  
SEPARABLE_FUNCTION void step(const dvariable& logU1, const dvariable& logU2, const dvariable& logSdU, const dvariable& logr0, const dvariable& K, int i)
  dvariable Penalty = 0.0;
  dvariable pred = posfun( exp(logU1) + exp(logU1)*exp(logr0)*(1.0 - exp(logU1)/K) - Catch(i-1), 0.001, Penalty );
  jnll -= ( -log(2*M_PI)/2 - log(pred) - logSdU - square(logU2-log(pred))/(2.0*mfexp(2.0*logSdU)) );
  jnll += square(Penalty);

SEPARABLE_FUNCTION void obs(const dvariable& logQ, const dvariable& logSdU, const dvariable& logU, int i)
  jnll -= ( -log(2*M_PI)/2 - (logQ+logU) - logSdU - square((logQ+logU)-log(BiomassIndex(i)))/(2.0*mfexp(2.0*logSdU)) );

SEPARABLE_FUNCTION void prior(const dvariable& logr0)
  jnll -= ( -log(2*M_PI)/2 - log(0.252) - log(0.51) - square(logr0-log(0.252))/(2.0*square(0.51)) );

TOP_OF_MAIN_SECTION
  // Set maximum number of independent variables to 1000 and increase memory.
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  arrmblsize=2000000;
