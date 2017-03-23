data {
	# nSubj: number of subjects
	int nSubj ;
	# nY: total number of observations
	int nY ;
	# nW: number of within-subject predictors
	int nW ;
	# nB: number of between-subject predictors
	int nB ;
	# Y: vector of outcomes
	vector[nY] Y ;
	# W: predictor matrix
	matrix[nY,nW] W ;
	# B: predictor matrix
	matrix[nSubj,nB] B ;
	# Subj: integer label of subject associated with each entry in Y
	int Subj[nY] ;
}
transformed data{
	# Y_scaled: z-transform of Y
	vector[nY] Y_scaled ;
	Y_scaled = (Y-mean(Y))/sd(Y) ;
}
parameters {
	# noise: measurement noise
	real<lower=0> noise ;
	# scaled_coef_means: coefficient means
	matrix[nB,nW] scaled_coef_means ;
	# scaled_coef_sds: coefficient sds
	vector<lower=0>[nW] scaled_coef_sds ;
	# cor_mat: correlations amongst coefficients
	corr_matrix[nW] cor_mat ;
	# scaled_subj_coefs: coefficients for each subject
	matrix[nSubj,nW] scaled_subj_coefs ;
}
model {
	# weakly informed priors
	to_vector(scaled_coef_means) ~ normal(0,1) ; #z-scale prior
	scaled_coef_sds ~ weibull(2,1) ; #z-scale prior
	cor_mat ~ lkj_corr(1) ; #flat prior on correlations
	noise ~ weibull(2,1) ; #prior on measuement noise
	#coefficients as multivariate normal
	for(s in 1:nSubj){
		scaled_subj_coefs[s,] ~ multi_normal(
			  B[s,] * scaled_coef_means
			, quad_form_diag(cor_mat, scaled_coef_sds)
		) ;
	}
	#loop over entries in the outcome vector Y
	for(i in 1:nY){
		Y_scaled[i] ~ normal(
			sum( scaled_subj_coefs[Subj[i]] .* W[i,] ) 
			, noise
		) ;
	}
}
generated quantities{
	# coef_means: coefficient means on the original scale of Y
	matrix[nB,nW] coef_means ;
	# coef_sds: coefficient sds on the original scale of Y
	vector[nW] coef_sds ;
	coef_sds = scaled_coef_sds * sd(Y) ;
	coef_means = scaled_coef_means * sd(Y) ;
	#adding back the mean to the intercept only
	coef_means[1,1] = coef_means[1,1] + mean(Y) ;
}