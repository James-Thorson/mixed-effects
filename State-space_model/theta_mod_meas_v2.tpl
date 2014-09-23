////////////////////////////////////////////////////////
// ADMB code for a measurement-error Schaefer surplus production model
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
  init_bounded_number logSdy(-4.6,0);
  init_bounded_number logr0(-4.6,0.183);  //# Bounds: log(0.01) - log(1.2) from Meyer and Millar 1999
  init_bounded_number logQ(-4.6,4.6);
  init_bounded_number K(10.0,2000.0);
  init_bounded_number logU1(0,6.9);
  objective_function_value jnll;
  
  // derived quantities
  vector Penalty(1,N); 

  // reporting quantities
  sdreport_number r0;
  sdreport_number Sdy;
  sdreport_number Q;
  sdreport_vector U(1,N);
  sdreport_vector qU(1,N);
  
PRELIMINARY_CALCS_SECTION
  logSdy = log(0.1);
  logr0 = log(0.2);
  logQ = log(1);
  K = 500;
  logU1 = log(200);

PROCEDURE_SECTION
  jnll = 0.0; // joint negative log-likelihood to be minimized.
  r0 = mfexp(logr0);
  Sdy = mfexp(logSdy);
  Q = mfexp(logQ);
  qU = Q*U;
  
  // transition equation (Eq. 3)
  U(1) = mfexp( logU1 );
  Penalty(1) = 0.0;
  for(int i=2; i<=N; i++){
    Penalty(i) = 0.0;
    U(i) = posfun( U(i-1) + U(i-1)*r0*(1.0-U(i-1)/K) - Catch(i-1), 0.001, Penalty(i) );
    jnll += square(Penalty(i));
  }
  // observation equation (Eq. 4)
  for(int i=1; i<=N; i++){
    jnll -= ( -log(2*M_PI)/2 - (logQ+log(U(i))) - logSdy - square((logQ+log(U(i)))-log(BiomassIndex(i)))/(2.0*mfexp(2.0*logSdy)) );
  }
  // prior on r
  jnll -= ( -log(2*M_PI)/2 - log(0.252) - log(0.51) - square(logr0-log(0.252))/(2.0*square(0.51)) );
  
TOP_OF_MAIN_SECTION
  // Set maximum number of independent variables to 1000 and increase memory.
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  arrmblsize=2000000;
