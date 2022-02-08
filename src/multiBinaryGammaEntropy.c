/*
  * The four routines for gamma distribution deviation splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  #include  "stdbool.h"
  static int *countn;
static int *tsplit;
static int number_of_responses;
static int number_of_ys;


static int right_n1[256];
static int left_n1[256];
static double right_sum[256];
static double left_sum[256];

bool debug_log = false;


int multiBinaryGammaEntropy_init(int n, double *y[], int maxcat, char **error,
                                 double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  number_of_ys = (int)(*parm);
  number_of_responses = 2*number_of_ys;
  *size = number_of_responses;
  
  if (debug_log)
    Rprintf("Training model with %d responses (%d outputs)", number_of_ys, number_of_responses);
  
  return 0;
}



  double gammaDeviationM(double y, double mean){
  return -log(y/mean) + (y - mean)/mean;
}

double binaryCrossEntropy(int n0, int n1, double p){
  if (p == 0 || p == 1) return 0;
  return (-n1*log(p) - n0*log(1-p)) ;
}



void
  multiBinaryGammaEntropy_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i, varindex;
  double meany, err, aux_err0, p;
  int n1;
  
  err = 0;  
  
  if (debug_log){
    Rprintf("\n \n _____________________ ", n);
    Rprintf("\n NODE with %d elements... \n", n);
  }
  
  for (varindex = 0; varindex < number_of_ys; varindex++){
    aux_err0 = 0;
    meany = 0.;
    n1 = 0;
    for (i = 0; i < n; i++) {
      if (y[i][varindex] > 0){
        n1++;
        meany += y[i][varindex];
      }
    }
    
    p = n1/((double)n);
    
    if (n1 > 0){
      meany = meany / n1; // continuous mean
      for (i = 0; i < n; i++){
        if (y[i][varindex] > 0) aux_err0 += gammaDeviationM(y[i][varindex], meany);
      }
      
      err += ( aux_err0 + binaryCrossEntropy(n-n1, n1, p) ); // =0 if n1==0

      
      if (debug_log) Rprintf("\n gammaDev = %lf; ",  aux_err0);
    } else {
      meany = 0;
      if (debug_log) Rprintf("\n gammaDev = %lf; ",  0.);
    }
    
    if (debug_log) Rprintf("binCE = %lf \n",  binaryCrossEntropy(n-n1, n1, p));
    
    
    value[2*varindex] = p;
    value[(2*varindex) + 1] = meany; 
    
    if (debug_log){
      Rprintf("   p = %lf", p);
      Rprintf(", meany = %lf", meany);
      Rprintf(", err = %lf \n", err);
    }

  }
  
  if (debug_log) Rprintf("\n TOTAL err = %lf \n", err);
  
  *risk = err;

}

/*
  Split function
*/

  void
  multiBinaryGammaEntropy_split(int n, double *y[], double *x, int nclass,
                                int edge, double *improve, double *split, int *csplit,
                                double myrisk, double *wt)
{
  int i, j, k, varindex;
  double temp, best;
  int direction = LEFT;
  int where = 0;
  
  double meany0, err0, aux_err0, p0;
  double left_mean, right_mean, left_err, right_err, aux_left_err, aux_right_err;
  int left_n, right_n;
  
  left_n = 0;
  right_n = n;
  err0 = 0; 
  
  for (varindex = 0; varindex < number_of_ys; varindex++){
    aux_err0 = 0;
    right_n1[varindex] = 0;
    left_n1[varindex] = 0;
    
    right_sum[varindex] = 0;
    left_sum[varindex] = 0;
    
    for (i = 0; i < n; i++) {
      if (y[i][varindex] > 0){
        right_sum[varindex] += y[i][varindex]; 
        right_n1[varindex]++;
      }
    }
    
    if (right_n1[varindex] > 0){
      meany0 = right_sum[varindex] / (double)right_n1[varindex];
      for (i = 0; i < n; i++) {
        if (y[i][varindex] > 0) aux_err0 += gammaDeviationM(y[i][varindex], meany0);
      }
    }
    

    err0 += (aux_err0 + binaryCrossEntropy(right_n-right_n1[varindex],
                                           right_n1[varindex], 
                                           right_n1[varindex]/( (double)right_n ))
            );
  }
  
  
  if (nclass == 0) {/* continuous predictor */
    temp = 0;
    best = 0;
      
for (i = 0;  right_n > edge; i++) {
  left_n++;
  right_n--;
  right_err = 0;
  left_err = 0;

  for (varindex = 0; varindex < number_of_ys; varindex++){
    if (y[i][varindex] > 0){
      left_n1[varindex]++;
      right_n1[varindex]--;
      left_sum[varindex] += y[i][varindex];
      right_sum[varindex] -= y[i][varindex];
    }
    
    if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
      aux_left_err = 0;
      if (left_n1[varindex] > 0){
        left_mean = left_sum[varindex]/( (double)left_n1[varindex] ); //positive
        for (k = 0; k <= i; k++){
          if (y[k][varindex]>0) aux_left_err += gammaDeviationM(y[k][varindex], left_mean);
        }
      }
      aux_left_err += binaryCrossEntropy(left_n-left_n1[varindex], 
                                         left_n1[varindex], 
                                         left_n1[varindex]/( (double)left_n ));
        
      aux_right_err = 0;
      if (right_n1[varindex] > 0){
        right_mean = right_sum[varindex]/( (double)right_n1[varindex] ); //positive
        for (k = n-1; k > i; k--) {
          if (y[k][varindex]>0) aux_right_err += gammaDeviationM(y[k][varindex], right_mean);
        }
      }
      aux_right_err += binaryCrossEntropy(right_n-right_n1[varindex],
                                          right_n1[varindex],
                                          right_n1[varindex]/( (double)right_n ));
    
          
      left_err += aux_left_err;
      right_err += aux_right_err;
    }
  }
  
  if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
    //temp = err0 - left_err*(left_n/(double)n) - right_err*(right_n/(double)n);  // >0 means improvement, i.e. LESS error
    //temp = err0 - (left_err + right_err)/2;  // >0 means improvement, i.e. LESS error
    temp = err0 - (left_err + right_err);  // >0 means improvement, i.e. LESS error
    
    if (debug_log){
      Rprintf("left_n = %d", left_n);
      Rprintf("; right_n = %d", right_n);
      Rprintf("; improvement = %lf", temp);
      Rprintf("; left_err = %lf", left_err);
      Rprintf("; right_err = %lf", right_err);
      Rprintf("\n");
    }
  
    if (temp > best) {
      best = temp;
      where = i;
      if (left_err < right_err)	direction = LEFT; 
      else direction = RIGHT;    
    }
  }
  
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
  multiBinaryGammaEntropy_pred(double *y, double *yhat)
{
    Rprintf("Detected xpred");
  //double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
  
  return 0;
}