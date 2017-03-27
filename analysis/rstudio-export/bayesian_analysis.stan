data {
  #nTrials: num trials
  int<lower=1> nTrials ;
  #power: power outcomes
  vector[nTrials] power ;
  #nS: num subjects
  int<lower=1> nS ;
  #S: trial-by-trial subject labels
  int<lower=1,upper=nS> S[nTrials] ;
  #nWmpower: num within predictors on mean power
  int<lower=1> nWmpower ;
  #nBmpower: num group predictors on mean power
  int<lower=1> nBmpower ;
  #nWspower: num within predictors on sd power
  int<lower=1> nWspower ;
  #nBspower: num group predictors on sd power
  int<lower=1> nBspower ;
  #Wmpower: within predictors for mean power
  matrix[nTrials, nWmpower] Wmpower ;
  #Bmpower: between predictors for mean power
  matrix[nS,nBmpower] Bmpower ;
  #Wspower: within predictors for sd power
  matrix[nTrials, nWspower] Wspower ;
  #Bspower: between predictors for sd power
  matrix[nS,nBspower] Bspower ;
}
transformed data{
  #nWtot: total number of within predictors
  int nWtot ;
  #compute nWtot
  nWtot = nWmpower + nWspower ;
}
parameters {
  #normal01: a helper variable
  matrix[nWtot, nS] normal01;
  #ZsdW: population-level sds (on the z-scale) for each within-subject predictor
  vector<lower=0>[nWtot] sdsW ;
  #ZcorW: population-level correlations (on cholesky factor scale) amongst within-subject predictors
  cholesky_factor_corr[nWtot] cors_helper ;
  #coefMpower: coefficientsfor between and within subject predictors on mpower
  matrix[nBmpower, nWmpower] coefMpower ;
  #coefSpower: coefficients for between and within subject predictors on spower
  matrix[nBspower, nWspower] coefSpower ;
}
model {
  #normal01 must have normal(0,1) prior for multivariate trick
  to_vector(normal01) ~ normal(0, 1) ;
  #flat prior on correlations
  cors_helper ~ lkj_corr_cholesky(1) ;
  #weibull prior on subject deviations
  sdsW ~ weibull(2, 1) ;
  #normal(0,1) priors on all mpower intercept & coefs
  to_vector(coefMpower) ~ normal(0, 1) ;
  #normal(0,1) priors on all spower intercept & coefs
  to_vector(coefSpower) ~ normal(0, 1) ;
  {
    #Sdevs: subject-by-subject deviations
    matrix[nS, nWtot] Sdevs ;
    #SvalsMpower: subject-by-subject mpower values
    matrix[nS, nWmpower] SvalsMpower ;
    #SvalsSpower: subject-by-subject spower values
    matrix[nS, nWspower] SvalsSpower ;
    #mpower: trial-by-trial mpower values
    vector[nTrials] mpower ;
    #spower: trial-by-trial spower values
    vector[nTrials] spower ;
    #compute
    Sdevs = transpose(diag_pre_multiply(sdsW,cors_helper) * normal01) ;
    SvalsMpower = Bmpower * coefMpower + Sdevs[,(1):(nWmpower)] ;
    SvalsSpower = Bspower * coefSpower + Sdevs[,(nWmpower+1):nWtot] ;
    mpower = rows_dot_product( SvalsMpower[S] , Wmpower ) ;
    spower = rows_dot_product( SvalsSpower[S] , Wspower ) ;
    target += normal_lpdf( power| mpower, exp(spower) );
  }
}
generated quantities {
  corr_matrix[nWtot] corsW;
  corsW = multiply_lower_tri_self_transpose(cors_helper);
}
