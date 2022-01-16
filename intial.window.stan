data{
  int	n_sample;														//	number	of	patients
  int p;                          // number of dimensions
  matrix[n_sample, p] X;              // covaraite matrix 
  vector [p] x;                     // the x to be predicted
  real	alpha;													//	for prior,	 lowbound for bandwidth
  real	beta;													//	for prior,	 upbound for bandwidth
}

parameters	{
  vector <lower=0> [p]	bandwidth;										//	bandwidth
}

// transformed parameters {
//   matrix[p, p] D; 
//   															
//   D = diag_matrix(bandwidth);
// }

model	{
  vector [n_sample] log_density;
  real lh ;

  for	(j	in	1:p)	bandwidth[j]~ inv_gamma(alpha, beta);	//	prior	on	bandwidth 
 
  for	(i	in	1:n_sample)	{

            log_density[i] = multi_normal_lpdf(x | X[i,], diag_matrix(bandwidth));
  }
  
  lh =  exp (log_sum_exp (log_density))/ n_sample ;
        
  target += log (lh)  ;
      
    
}
  

