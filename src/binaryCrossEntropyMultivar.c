
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"

static int *countn;
static int *tsplit;
static int number_of_responses;
  
int binaryCrossEntropyMultivar_init(int n, double *y[], int maxcat, char **error,
                            double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  
  number_of_responses = (int)(*parm);
  *size = number_of_responses;
  
  //Rprintf("Training model with %d responses...\n", number_of_responses);
  return 0;
  
}



double negBinaryCrossEntropyMultivar(double p){
  if (p == 0 || p == 1) return 0;
  return -p*log(p) - (1-p)*log(1-p);
}



void
binaryCrossEntropyMultivar_eval(int n, 
                                double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i, varindex;
  double err = 0;
  double p;
  int n1;
  
  for (varindex = 0; varindex < number_of_responses; varindex++){
    n1 = 0;
    
    for (i = 0; i < n; i++) if (y[i][varindex] > 0) n1++; 
    p = n1/((double)n);
    
    err += n*negBinaryCrossEntropyMultivar(p); // =0 if n1==0 | n1=n
    
    value[varindex] = p;
    //Rprintf("p = %lf, ", value[varindex]);
  }
  //Rprintf("\n");
  
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
binaryCrossEntropyMultivar_split(int n, double *y[], double *x, int nclass,
                         int edge, double *improve, double *split,
                         int *csplit,
                         double myrisk, double *wt)
{
  int i, j, k, varindex;
  double temp, best;
  int direction = LEFT;
  int where = 0;
  
  double err0 = 0;
  double left_err = 0;
  double right_err = 0;
  int left_n = 0;
  int right_n;
  
  int *right_n1 = malloc(number_of_responses*sizeof(int));
  int *left_n1 = malloc(number_of_responses*sizeof(int));
  
  for (varindex = 0; varindex < number_of_responses; varindex++){
    right_n1[varindex] = 0;
    left_n1[varindex] = 0;
    
    for (i = 0; i < n; i++){
      if (y[i][varindex] > 0){
        right_n1[varindex]++;
      }
    }
    
    /*
    Rprintf("Var: %d = ", varindex);
    Rprintf("right_n1: %d; ", right_n1[varindex]);
    Rprintf("left_n1: %d \n", left_n1[varindex]);
    */
    
    err0 += n*negBinaryCrossEntropyMultivar(right_n1[varindex]/((double)n)); // =0 if n1==0 | n1=n
  }
  
  //Rprintf("err0=%lf \n" , err0);
  
  if (nclass == 0) {/* continuous predictor */
    temp = 0;
    best = 0;
    
    right_n = n;
    
      for (i = 0;  right_n > edge; i++) {
        right_err = 0; 
        left_err = 0;
        left_n++;
        right_n--;
        
        //Rprintf("i: %d = \n", i);
        
        for (varindex = 0; varindex < number_of_responses; varindex++){
            
          if (y[i][varindex] > 0){
            left_n1[varindex]++;
            right_n1[varindex]--;
          }
          
          /*
          Rprintf("Var: %d = ", varindex);
          Rprintf("right_n1: %d; ", right_n1[varindex]);
          Rprintf("left_n1: %d \n", left_n1[varindex]);
          */
          if (x[i + 1] != x[i] && left_n >= edge) {
            left_err += (left_n * negBinaryCrossEntropyMultivar( left_n1[varindex]/((double)left_n) ));
            right_err += (right_n * negBinaryCrossEntropyMultivar( right_n1[varindex]/((double)right_n) ));
          }
        }
        
        /*
        Rprintf("left_err: %lf ; ", left_err);
        Rprintf("right_err: %lf \n", right_err);
        */
        
        if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
          temp = err0 - left_err - right_err;  // >0 means improvement, i.e. LESS deviance
        
          if (temp > best) {
            //Rprintf("Improvement %lf ...\n", temp);
            best = temp;
            where = i;
            if (left_n < right_n)	direction = LEFT; 
            else direction = RIGHT;
          }
        }
    }
      
      free(right_n1);
      free(left_n1);  
      
      *improve = best;
      if (best > 0) {         /* found something */

          //Rprintf("BEST Improvement %lf ... in %d \n", best, where);
        
          csplit[0] = direction;
          *split = (x[where] + x[where + 1]) / 2;
      }
      
  } else { // categorical
    error("Categorical predictors detected in binaryCrossEntropy");
  }
  
}

double
binaryCrossEntropyMultivar_pred(double *y, double *yhat)
{
  Rprintf("Detected xpred");
  double misc = 0; //abs(y[i][0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return misc;
}








