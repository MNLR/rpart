
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"

static int *countn;
static int *tsplit;

int binaryCrossEntropy_init(int n, double *y[], int maxcat, char **error,
                            double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  *size = 1;  // p
  return 0;
}



double negBinaryCrossEntropy(double p){
  if (p == 0 || p == 1) return 0;
    return -p*log(p) - (1-p)*log(1-p);
}



void
  binaryCrossEntropy_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double err = 0;
    double p;
    int n1;
    
    n1 = 0;
    for (i = 0; i < n; i++)  if (*y[i] > 0) n1++; 
    
    p = n1/((double)n);
    
    err = negBinaryCrossEntropy(p); // =0 if n1==0 | n1=n

    *value = p;
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
  binaryCrossEntropy_split(int n, double *y[], double *x, int nclass,
                                         int edge, double *improve, double *split,
                                         int *csplit,
                                         double myrisk, double *wt)
  {
    int i, j, k;
    double temp, best;
    int direction = LEFT;
    int where = 0;
    
    double err0;
    double left_err = 0;
    double right_err = 0;
    int left_n, right_n, right_n1, left_n1;
    
    right_n1 = 0;
    for (i = 0; i < n; i++) if (*y[i] > 0) right_n1++;
    
    
    err0 = n * negBinaryCrossEntropy( right_n1/((double)n) );
    
    if (nclass == 0) {/* continuous predictor */
      temp = 0;
      right_n = n;
      left_n = 0;
      left_n1 = 0;
      best = 0;

      for (i = 0;  right_n > edge; i++) {
        left_n++;
        right_n--;
        
        if (*y[i] > 0){
          left_n1++;
          right_n1--;
        }
        
        if (x[i + 1] != x[i] && left_n >= edge) {
          left_err = left_n * negBinaryCrossEntropy( left_n1/((double)left_n) );
          
          right_err = right_n * negBinaryCrossEntropy( right_n1/((double)right_n) );
          
          temp = err0 - left_err - right_err;  // >0 means improvement, i.e. LESS deviance
          if (temp > best) {
            best = temp;
            where = i;
            if (left_n1 < right_n1)	direction = LEFT; 
            else direction = RIGHT;
          }
        }
        
        /*

         Rprintf(", left_n1=%d; ", left_n1);
         Rprintf(", right_n1=%d; ", right_n1);
         
         Rprintf(", nbces_right = %lf", 
         negBinaryCrossEntropy( right_n1/( (double)right_n )));
         
         Rprintf(", nbces_left = %lf", negBinaryCrossEntropy(
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
            error("Categorical predictors detected in binaryCrossEntropy");
    }
          
  }

double
  binaryCrossEntropy_pred(double *y, double *yhat)
  {
    Rprintf("Detected xpred");
    double misc = abs(y[0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
    
    return misc;
  }