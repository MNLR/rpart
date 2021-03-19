/*
  * The four routines for gamma Log Likelihood Moments Matching splitting routine
*/

#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
  
static int *countn;
static int *tsplit;

int gammaLLMME_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaLLMME_init: Categorical predictors detected. This is not yet implemented for gammaLLMME");
  }
  
  *size = 3;  // First value is alpha/beta. 2nd: alpha; 3rd: beta
  return 0;
}


double gammaLogLikelihood(double y, double shape, double rate){
  double density_gamma = shape*log(rate) + (shape-1)*log(y) - rate*y - lgamma(shape) ;
  return density_gamma;
}

  /*
   * The gamma Log Likelihood Moments Matching evaluation function.
   * Negative Log Likelihood -> Lower is better (can be negative and positive)
   * assigns:
   * *risk = Negative Log Likelihood
   * value[0] = mean = shape/rate 
   * value[1] = shape 
   * value[2] = rate
   */
  
  void
  gammaLLMME_eval(int n, double *y[], double *value, double *risk, double *wt)
{
    // weights disabled
  int i;
  double meany, sdy, shape, rate, gamma_negloglikelihood;

  meany = 0.;
  for (i = 0; i < n; i++) {
    meany += *y[i];
  }
  meany = meany / n;
  
  
    sdy = 0.;
    for (i = 0; i < n; i++) {
      sdy += (*y[i] - meany) * (*y[i] - meany);
    }
    
    if (sdy == 0) sdy = meany; // degenerated caseS: 1 element - all elements equal
    else sdy /= (n-1);
  
    shape = (meany*meany)/sdy;
    rate = shape/meany;
  
  
  gamma_negloglikelihood = 0;
  for (i = 0; i < n; i++)
    gamma_negloglikelihood -= gammaLogLikelihood(*y[i], shape, rate);
  
  // 
  
  value[0] = meany;
  value[1] = shape;
  value[2] = rate;
  *risk = gamma_negloglikelihood;
  
}

/*
Gamma spliting functon
*/
  void
  gammaLLMME_split(int n, double *y[], double *x, int nclass,
      int edge, double *improve, double *split, int *csplit,
      double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum, right_sum;
  int left_n, right_n;
  int direction = LEFT;
  int where = 0;
  
  
  double meany0, sdy0, shape0, rate0, gamma_loglikelihood0;
  double left_mean, left_sd, right_mean, right_sd, left_gamma_logLikelihood, right_gamma_logLikelihood, aux;
  
  right_sum = 0;
  for (i = 0; i < n; i++) right_sum += *y[i];
  meany0 = right_sum / n;
  
  sdy0 = 0;
  for (i = 0; i < n; i++) sdy0 += (*y[i] - meany0) * (*y[i] - meany0);
  
  if (sdy0 == 0){ // Either one element or all elements are the same value. No point in spliting further
    *improve = 0;
  } else {
    sdy0 /= n-1; // sample variance
    gamma_loglikelihood0 = 0;
    for (i = 0; i < n; i++) gamma_loglikelihood0 += gammaLogLikelihood(*y[i], (meany0*meany0)/sdy0, meany0/sdy0);

    right_n = n;
  
    if (nclass == 0) {/* continuous predictor */
      temp = 0;
      left_sum = 0;           /* No data in left branch, to start, right_sum is sum(y) */
      left_n = 0;
	    best = 0;
	  
	    for (i = 0;  right_n > edge; i++) {
	      left_n++;
	      right_n--;
	      left_sum += *y[i];
	 	    right_sum -= *y[i];
	 	  
	      if (x[i + 1] != x[i] && left_n >= edge) {
	      
	        left_mean = left_sum/left_n;
	        right_mean = right_sum/right_n;
	      

	      left_sd = 0;
	      for (k = 0; k <= i; k++) left_sd += (*y[k] - left_mean) * (*y[k] - left_mean);
	      if (left_sd != 0) left_sd = left_sd/(left_n - 1); // sample variance
	      else left_sd = left_mean;

        right_sd = 0;
        for (k = n-1; k > i; k--) right_sd += (*y[k] - right_mean) * (*y[k] - right_mean);
        if (right_sd != 0) right_sd = right_sd/(right_n - 1); // sample variance
        else right_sd = right_mean;

	      left_gamma_logLikelihood = 0;
	      for (k = 0; k <= i; k++)
	        left_gamma_logLikelihood += gammaLogLikelihood(*y[k], (left_mean*left_mean)/left_sd, left_mean/left_sd);
	      
	      right_gamma_logLikelihood = 0;
	      for (k = n-1; k > i; k--) 
	        right_gamma_logLikelihood += gammaLogLikelihood(*y[k], (right_mean*right_mean)/right_sd, right_mean/right_sd);

	      
	      temp = left_gamma_logLikelihood + right_gamma_logLikelihood - gamma_loglikelihood0;
		    if (temp > best) {
		      best = temp;
		      where = i;
		       if (left_sum < right_sum)	direction = LEFT; 
		       else direction = RIGHT;
		    }
	    }
	    }
	  

	    *improve = best;
	    if (best > 0) {         
	      csplit[0] = direction;
	      *split = (x[where] + x[where + 1]) / 2;
	    }
	  
    }
  
    /*
     * Categorical predictor
     */
    else {
      error("In gammaLLMME_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
    }
  }
}

double
  gammaLLMME_pred(double *y, double *yhat)
  {
    Rprintf("Predicting out-of-sample... \n");
    double meanpred = fabs(y[0] - *yhat);  // y[0] is obs, *yhat is first value of prediction

    return meanpred;
  }
