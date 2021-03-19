/*
  * The four routines for gamma Log Likelihood Moments Matching splitting routine
*/
  
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
  
static int *countn;
static int *tsplit;

int gammaLLmean_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaLLmean_init: Categorical predictors detected. This is not yet implemented for gammaLLMME");
  }
  
  *size = 1;  // 
  return 0;
}



double gammaLogLikelihoodFixedRate(double y, double shape){
  double density_gamma =  (shape-1)*log(y) - y - lgamma(shape) ;
  return density_gamma;
}

/*
  * The gamma Log Likelihood - Mean evaluation function. 
  * Assumes rate = 1 ~ exponential distribution 
  * Negative Log Likelihood -> Lower is better (can be negative and positive)
  * assigns:
  * *risk = Negative Log Likelihood
  * *value = mean
*/
  
  void
  gammaLLmean_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i;
  double meany, gamma_negloglikelihood;
  
  meany = 0.;
  for (i = 0; i < n; i++) {
    meany += *y[i];
  }
  meany = meany / n;
  
  gamma_negloglikelihood = 0;
  for (i = 0; i < n; i++)
    gamma_negloglikelihood -= gammaLogLikelihoodFixedRate(*y[i], meany);
  
  *value = meany;
  *risk = gamma_negloglikelihood;
  
}

/*
  Gamma spliting functon
*/
  void
  gammaLLmean_split(int n, double *y[], double *x, int nclass,
                 int edge, double *improve, double *split, int *csplit,
                 double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum, right_sum;
  int left_n, right_n;
  int direction = LEFT;
  int where = 0;
  
  
  double meany0, gamma_loglikelihood0;
  double left_mean, right_mean, left_gamma_logLikelihood, right_gamma_logLikelihood;
  
  right_sum = 0;
  for (i = 0; i < n; i++) right_sum += *y[i];
  meany0 = right_sum / n;
  
  gamma_loglikelihood0 = 0;
  for (i = 0; i < n; i++)
    gamma_loglikelihood0 += gammaLogLikelihoodFixedRate(*y[i], meany0);
  
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
            
          left_gamma_logLikelihood = 0;
          for (k = 0; k <= i; k++)
            left_gamma_logLikelihood += gammaLogLikelihoodFixedRate(*y[k], left_mean);
              
          right_gamma_logLikelihood = 0;
          for (k = n-1; k > i; k--) 
            right_gamma_logLikelihood += gammaLogLikelihoodFixedRate(*y[k], right_mean);

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
      if (best > 0) {         /* found something */
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

double
  gammaLLmean_pred(double *y, double *yhat)
{
  Rprintf("Predicting out-of-sample... \n");
  double meanpred = fabs(y[0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}
