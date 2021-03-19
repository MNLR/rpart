/*
 * The four routines for gamma distribution deviation splitting routine
 */

#include <math.h>
#include "rpart.h"
#include "rpartproto.h"

static int *countn;
static int *tsplit;

int bernoulliLL_init(int n, double *y[], int maxcat, char **error,
                     double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In bernoulliLL_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
  }
  
  *size = 2;  // First value: p, second: n
  return 0;
}



void
  bernoulliLL_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double sumy, bernoulli_nll;
    double ndoub;
    
    ndoub = (double) n;
      
    sumy = 0;
    for (i = 0; i < n; i++) {
      sumy += *y[i];
    }
    
    // NEGATIVE Log Likelihood (always positive)
    bernoulli_nll = 0;  
    if (sumy != 0 && sumy != n ){ // bernoulli_nll = 0 means perfect node
      for (i = 0; i < n; i++) {
        bernoulli_nll -= sumy*log(sumy/ndoub) + (ndoub-sumy)*( log(1-(sumy/ndoub)) );
      }
    }
    

    /*
    Rprintf("\n ________________________ ");
    
    Rprintf("Node bernoulli_nll: ");
    Rprintf("%lf", bernoulli_nll);
    Rprintf("\n");  
     */
    
    value[0] = sumy/ndoub;
    value[1] = ndoub;
    *risk = bernoulli_nll;
  }


/*
 Bernoulli Log Likelihood spliting functon
 */
void
  bernoulliLL_split(int n, double *y[], double *x, int nclass,
                       int edge, double *improve, double *split, int *csplit,
                       double myrisk, double *wt)
  {
    int i, j, k;
    double temp, best;
    double left_sum, right_sum;
    double left_n, right_n;
    int direction = LEFT;
    int where = 0;
    
    double bernoulli_ll0, left_bernoulli_ll, right_bernoulli_ll;;

    if (nclass != 0) {
      error("In bernoulliLL_split: Categorical predictors detected. This is not yet implemented");
    } else { // continuous predictor
      right_sum = 0;
      for (i = 0; i < n; i++) right_sum += *y[i];
      
      if (right_sum != 0 && right_sum != n ) { // otherwise no point in spliting further
        bernoulli_ll0 = 0;
        for (i = 0; i < n; i++) {
bernoulli_ll0 += right_sum*log(right_sum/n) + (n-right_sum)*( log(1-(right_sum/n)) );
      }
        
      bernoulli_ll0 /= n;
    
    
/*
     Rprintf("bernoulli_ll0: ");
     Rprintf("%lf", bernoulli_ll0);
     Rprintf("\n");  
*/
      right_n = (double) n;
    
      temp = 0;
      left_sum = 0;           /* No data in left branch, to start, right_sum is sum(y) */
      left_n = 0;
      best = 0;
    
    for (i = 0;  right_n > edge; i++) {
      left_n += 1;
      right_n--;
      left_sum += *y[i];
      right_sum -= *y[i];
      
      if (x[i + 1] != x[i] && left_n >= edge) {
        
        
  /*      
         Rprintf("left_sum: ");
         Rprintf("%lf", left_sum);
         Rprintf("\n");  
         Rprintf("right_sum: ");
         Rprintf("%lf", right_sum);
         Rprintf("\n");  
    */    
         
        left_bernoulli_ll = 0;
        if (left_sum != 0 && (left_sum != left_n)){ // else left_bernoulli_ll = 0
          for (k = 0; k <= i; k++) {
left_bernoulli_ll += left_sum*log(left_sum/left_n) + (left_n-left_sum)*( log(1-(left_sum/left_n)) );
          }
          left_bernoulli_ll /= left_n;
        }
        
        /*
         Rprintf("left_bernoulli_ll: ");
         Rprintf("%lf", left_bernoulli_ll);
         Rprintf("\n");  
         */
        
        right_bernoulli_ll = 0;
        if (right_sum != 0 && (right_sum != right_n)){
          for (k = n-1; k > i; k--) {
right_bernoulli_ll += right_sum*log(right_sum/right_n) + (right_n-right_sum)*( log(1-(right_sum/right_n)) );
          }
          right_bernoulli_ll /= right_n;
        }
        
        /*
         Rprintf("right_bernoulli_ll: ");
         Rprintf("%lf", right_bernoulli_ll);
         Rprintf("\n");  
         */
        temp = left_bernoulli_ll + right_bernoulli_ll - bernoulli_ll0;  // >0 means improvement, i.e. LESS deviance
        
        /*
         Rprintf("temp: ");
         Rprintf("%lf", temp);
         Rprintf("\n");  
         */
         
        if (temp > best) {
          best = temp;
          where = i;
          if (left_sum < right_sum)	direction = LEFT; 
          else direction = RIGHT;
        }
      } //if (x[i + 1] != x[i] && left_n >= edge)
    } //for (i = 0;  right_n > edge; i++)
    
    
      *improve = best;
      if (best > 0) {         /* found something */
        csplit[0] = direction;
        *split = (x[where] + x[where + 1]) / 2;
      }
    
    } else { // right_sum == 0 || right_sum == n 
      *improve = 0;
    }
    }
    

  }

double // the NEGATIVE LogLikelihood (lower is better)
  bernoulliLL_pred(double *y, double *yhat)
  {
    Rprintf("Predicting out-of-sample... \n");
    double meanpred = abs(*yhat - y[0]); 

    return meanpred;
  }
