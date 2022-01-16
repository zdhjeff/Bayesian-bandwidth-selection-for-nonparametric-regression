data{
  int	n_sample;														//	number	of	patients
  int p;                          // number of dimensions
  matrix[n_sample, p] X;              // covaraite matrix 
  real	l0;													//	for prior,	 lowbound for bandwidth
  //real	u0;													//	for prior,	 upbound for bandwidth
  real	Y[n_sample];			//	vector	of	responses
  //real Yf[n_sample];   // same as y but for the loop 
}

parameters	{
  vector <lower=0> [p]	bandwidth;										//	bandwidth
  real<lower=0> ysd;         //  the sd for y , assuming the sd is the same bewtween each treatment, i.e. homogeneity of varaince 
}

// transformed parameters {
//   matrix[p, p] D; 
//   															
//   D = diag_matrix(bandwidth);
// }

model	{
  matrix[n_sample,n_sample] log_density;
  matrix[n_sample, n_sample] reg; 
  real denominator[n_sample];
  real numerator[n_sample];
  vector [n_sample]	m_hat;		
  
  for	(j	in	1:p)	bandwidth[j]~	exponential(l0);	//	prior	on	bandwidth 
  1/ysd ~ inv_gamma(1, 1);  // prior for the variance 
  
  for	(i	in	1:n_sample)	{

      for (j in 1: n_sample) {

            log_density[i,j] = multi_normal_lpdf(X[i,] | X[j,], diag_matrix(bandwidth));
            log_density[i,i]=-20;
            
         //    (2 * 3.14)^(-p/2) *  sqrt(1/prod(bandwidth)) 
         //   * exp( (-1/2) 
         //   // *  ( (X[i,]-X[j,])' * ((X[i,]-X[j,]) .* (1 ./ bandwidth) )  ) ;
         // // * quad_form( inverse(D),  (X[i,]-X[j,] ) ) );
         //    
         //   * ( (X[i,]-X[j,])  * diag_matrix(inv(bandwidth)) * (X[i,]-X[j,])'   )  );   

            reg[i,j] = exp (log_density[i,j]) * Y[j];
           // print("reg[i,j]",reg[i,j]);
      }

      denominator[i]= log_sum_exp( log_density[i,]);
      
      numerator[i] = sum(reg[i,]);
      m_hat[i]= numerator[i]/ exp(denominator[i]) ;
     // print("denominator[i]",denominator[i]);
      // print("numerator[i]",numerator[i]);
      Y[i] ~ normal (m_hat[i], ysd);  //	likelihood
    
  }
  
}

