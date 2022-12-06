/*
 * The four routines for gamma distribution deviation splitting routine
 */

#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
#include <stdbool.h>

static bool log_debug = false;

static int *countn;
static int *tsplit;

int mae_init(int n, double *y[], int maxcat, char **error,
                        double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In mae_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
  }
  
  
  *size = 1;  // mean value
  if (log_debug) Rprintf("Training model with Absolute Error (%d variables) \n", *size);
  
  return 0;
}



void
  mae_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double meany, mae;
    
    meany = 0.;
    for (i = 0; i < n; i++) {
      meany += *y[i];
    }
    meany = meany / n;
    
    mae = 0;  
    for (i = 0; i < n; i++) {
      mae += fabs(*y[i] - meany);
    }
    
    if (log_debug){
      Rprintf("###\n");
      Rprintf("Node with %d elements.\n", n);
      Rprintf("Mean = %lf. AbsError = %lf \n", 
              meany, 
              mae);
    }
    
    *value = meany;
    *risk = mae;
    
  }

/*
 Gamma spliting functon
 */
void
  mae_split(int n, double *y[], double *x, int nclass,
                       int edge, double *improve, double *split, int *csplit,
                       double myrisk, double *wt)
  {
    int i, j, k;
    double temp, best;
    double left_sum, right_sum;
    int left_n, right_n;
    int direction = LEFT;
    int where = 0;
    
    
    double meany0, mae0 = 0;
    double left_mean, right_mean, left_mae, right_mae;
    
    right_sum = 0;
    for (i = 0; i < n; i++) right_sum += *y[i];
    meany0 = right_sum / n;
    
    mae0 = 0;  
    for (i = 0; i < n; i++) mae0 += fabs(*y[i] - meany0);

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
    
    left_mae = 0;
    for (k = 0; k <= i; k++) left_mae += fabs(*y[k] - left_mean);
      
    right_mae = 0;
    for (k = n-1; k > i; k--) right_mae += fabs(*y[k] - right_mean);

    temp = mae0 - left_mae - right_mae;  // >0 means improvement, i.e. LESS deviance
    
    
    if (log_debug){
      Rprintf("Candidate split with number of elements: (LEFT,RIGHT) = (%d,%d)\n", 
              left_n , right_n);
      
      Rprintf("0 - (left - right) = %lf - (%lf + %lf). AE = %lf\n", 
              mae0,
              left_mae,
              right_mae,
              temp);
    }
    
    
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
    if (log_debug) Rprintf("BEST Improvement %lf ... in %d \n \n \n", 
        best, where + 1); // I add 1 for R
    
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
  mae_pred(double *y, double *yhat)
  {
    Rprintf("XPRED: \n");
    double meanpred =  fabs(y[0]-(*yhat)) ;  // y[0] is obs, *yhat is first value of prediction
    return meanpred;
  }
