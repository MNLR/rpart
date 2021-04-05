/*
  * The four routines for gamma distribution deviation splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  
  static int *countn;
  static int *tsplit;

int binaryCrossEntropyGammaDeviation_init(int n, double *y[], int maxcat, char **error,
                                          double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  *size = 2;  // mean value
  return 0;
}


double gammaDeviation(double y, double mean){
  return - log(y/mean) + (y - mean)/mean;
}

double negBinaryCrossEntropySum(int n0, int n1, double p){
  if (p == 0 || p == 1) return 0;
    return -n1*log(p) - n0*log(1-p);
}

void
  binaryCrossEntropyGammaDeviation_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i;
  double meany, err, p;
  int n1;
  
  meany = 0.;
  n1 = 0;
  for (i = 0; i < n; i++) {
    if (*y[i] > 0){
      n1++;
      meany += *y[i];
    }
  }
  
  err = 0;  
  p = n1/((double)n);
  
  if (n1 > 0){
    meany = meany / n1; // continuous mean
    for (i = 0; i < n; i++){
      if (*y[i] > 0) err += gammaDeviation(*y[i], meany);
    }
    
    err += negBinaryCrossEntropySum(n-n1, n1, p); // =0 if n1==0
  } else meany = 0;
  
  value[0] = p;
  value[1] = meany; // todo code to NaN
  *risk = err;
  
  /*
  Rprintf("p = %lf", p);
  Rprintf(", meany = %lf", meany);
  Rprintf(",err = %lf", err);
  Rprintf(", nbces = %lf", negBinaryCrossEntropySum(n-n1, n1, p));
  */
  
}

/*
  Split function
*/
  void
  binaryCrossEntropyGammaDeviation_split(int n, double *y[], double *x, int nclass,
                     int edge, double *improve, double *split, int *csplit,
                     double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum, right_sum;
  int direction = LEFT;
  int where = 0;
  
  double meany0, err0, p0;
  double left_mean, right_mean, left_err, right_err;
  int left_n, right_n, right_n1, left_n1;

  
  right_sum = 0;
  right_n1 = 0;
  for (i = 0; i < n; i++) {
    if (*y[i] > 0){
      right_sum += *y[i]; 
      right_n1++;
    }
  }
  
  err0 = 0; //shouldn't enter here with meany0 = 0, as err=0
  if (right_n1 > 0){
    meany0 = right_sum / right_n1;
    for (i = 0; i < n; i++) {
      if (*y[i]>0) err0 += gammaDeviation(*y[i], meany0);
    }
  }
  right_n = n;

  err0 += negBinaryCrossEntropySum(right_n-right_n1, 
                                   right_n1, 
                                   right_n1/( (double)right_n )
                                   );
  if (nclass == 0) {/* continuous predictor */
      temp = 0;
      left_sum = 0;           /* No data in left branch, to start, right_sum is sum(y) */
      left_n = 0;
      left_n1 = 0;
      best = 0;
        
      for (i = 0;  right_n > edge; i++) {
          left_n++;
          right_n--;
          
          if (*y[i] > 0){
            left_n1++;
            right_n1--;
            left_sum += *y[i];
            right_sum -= *y[i];
          }
          
          if (x[i + 1] != x[i] && left_n >= edge) {
            left_err = 0;
            if (left_n1 > 0){
              left_mean = left_sum/left_n1; //positive
              for (k = 0; k <= i; k++){
                if (*y[k]>0) left_err += gammaDeviation(*y[k], left_mean);
              }
            }
            left_err += negBinaryCrossEntropySum(left_n-left_n1, 
                                                 left_n1, 
                                                 left_n1/( (double)left_n ));
              
            right_err = 0;
            if (right_n1 > 0){
              right_mean = right_sum/right_n1; //positive
              for (k = n-1; k > i; k--) {
                if (*y[k]>0) right_err += gammaDeviation(*y[k], right_mean);
              }
            }
            right_err += negBinaryCrossEntropySum(right_n-right_n1,
                                                  right_n1,
                                                  right_n1/( (double)right_n ));
                
            temp = err0 - left_err - right_err;  // >0 means improvement, i.e. LESS deviance
                if (temp > best) {
                  best = temp;
                  where = i;
                  if (left_sum < right_sum)	direction = LEFT; 
                  else direction = RIGHT;
                }
          }
          
          /*
          Rprintf(", left_sum=%lf; ", left_sum);
          Rprintf(", right_sum=%lf; ", right_sum);
          
          Rprintf(", left_mean=%lf; ", left_mean);
          Rprintf(", right_mean=%lf; ", right_mean);
          Rprintf(", left_n1=%d; ", left_n1);
          Rprintf(", right_n1=%d; ", right_n1);

          Rprintf(", nbces_right = %lf", negBinaryCrossEntropySum(right_n-right_n1,
                                               right_n1,
                                               right_n1/( (double)right_n )));
          
          Rprintf(", nbces_left = %lf", negBinaryCrossEntropySum(left_n-left_n1, 
                                               left_n1, 
                                               left_n1/( (double)left_n )));
          

          Rprintf(";%lf", err0);
          Rprintf(",%lf", left_err);
          Rprintf(",%lf; ", right_err);
          Rprintf("\n");
           */
        }
        
        *improve = best;
        if (best > 0) {         /* found something */
            csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2;
        }
        
  } else { // categorical
    error("Categorical predictors detected in binaryCrossEntropyGammaDeviation");
  }
  
}

double
  binaryCrossEntropyGammaDeviation_pred(double *y, double *yhat)
{
    Rprintf("Detected xpred");
  double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}