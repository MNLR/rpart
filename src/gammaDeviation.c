/*
  * The four routines for gamma distribution deviation splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  
  static int *countn;
static int *tsplit;

int gammaDeviation_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaDeviation_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
  }
  
  *size = 1;  // mean value
  return 0;
}



  void
gammaDeviation_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i;
  double meany, gamma_dev;
  
  meany = 0.;
  for (i = 0; i < n; i++) {
    meany += *y[i];
  }
  meany = meany / n;

  gamma_dev = 0;  
  for (i = 0; i < n; i++) {
    gamma_dev += - log(*y[i]/meany) + (*y[i] - meany)/meany;
  }
  
  *value = meany;
  *risk = gamma_dev;
  
}

/*
  Gamma spliting functon
*/
  void
gammaDeviation_split(int n, double *y[], double *x, int nclass,
                 int edge, double *improve, double *split, int *csplit,
                 double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum, right_sum;
  int left_n, right_n;
  int direction = LEFT;
  int where = 0;
  
  
  double meany0, gamma_dev0 = 0;
  double left_mean, right_mean, left_gamma_deviance, right_gamma_deviance;
  
  right_sum = 0;
  for (i = 0; i < n; i++) right_sum += *y[i];
  meany0 = right_sum / n;
  
  gamma_dev0 = 0;  
  for (i = 0; i < n; i++) {
    gamma_dev0 += - log(*y[i]/meany0) + (*y[i] - meany0)/meany0;
  }

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

              left_gamma_deviance = 0;
              for (k = 0; k <= i; k++)
                left_gamma_deviance += - log(*y[k]/left_mean) + (*y[k] - left_mean)/left_mean;

              right_gamma_deviance = 0;
              for (k = n-1; k > i; k--) 
                right_gamma_deviance +=  - log(*y[k]/right_mean) + (*y[k] - right_mean)/right_mean;
                
              temp = gamma_dev0 - left_gamma_deviance - right_gamma_deviance;  // >0 means improvement, i.e. LESS deviance
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
        error("In gammaDeviation_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
      }
  
}

double
gammaDeviation_pred(double *y, double *yhat)
{
  double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}
