
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
#include <stdbool.h>
static int *countn;
static int *tsplit;
static int number_of_responses;
static int n_correlations;
static double cor_penalization;
static double pred_penalization;
static double cor_obs[4950]; // just enough for 100 responses
static double default_cor_penalization = 0; // .5*.5
static bool print_log = true;


int binaryCrossEntropyMultivarCorPenalization_init(int n, double *y[], int maxcat, 
                                                   char **error,
                                                   double *parm, int *size, int who, 
                                                   double *wt)
{
  /*
   * parm[0] number of responses
   * parm[1] alpha 
   * parm[2] beta
   * parm[3:(2 + ( ((parm[0]-1)*parm[0])/2 )] cor values
   * 
   */
  
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  
  number_of_responses = (int)(parm[0]);
  *size = number_of_responses;
  
  pred_penalization = (double)parm[1];
  cor_penalization = (double)parm[2];
  
  n_correlations = ((number_of_responses-1)*number_of_responses)/2 ;
  if (n_correlations > 4950) *error = _("Maximum number of responses is 100");
  
  for (int i = 0; i < n_correlations; i++){
    cor_obs[i] = (double)parm[3+i];
  }
  
  if (print_log){
    Rprintf("Training model with %d responses...\n", number_of_responses);
    Rprintf("Penalizations (a = %lf, b = %lf) \n", pred_penalization, cor_penalization);
    Rprintf("Correlation values: \n");
    for (int i = 0; i < n_correlations; i++){
      Rprintf("%lf, ", cor_obs[i]);
    }
  }

  return 0;
  
} 



double negBinaryCrossEntropyMultivarCorPenalization(double p){
  if (p == 0 || p == 1) return 0;
  return -p*log(p) - (1-p)*log(1-p);
}



void
  binaryCrossEntropyMultivarCorPenalization_eval(int n, 
                                double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i, varindex;
  double err = 0;
  double cor_err = 0;
  double p;
  int corindex;
  static double right_sumprodxy[4950]; // just enough for 100 responses
  int *n1 = malloc(number_of_responses*sizeof(int));
  
  if (print_log) Rprintf("\n --- NODE ---\n");
  
  if (print_log)  Rprintf("p = ");
  for (varindex = 0; varindex < number_of_responses; varindex++){
    n1[varindex] = 0;
    
    for (i = 0; i < n; i++) if (y[i][varindex] > 0) n1[varindex]++; 
    p = n1[varindex]/((double)n);
    
    err += negBinaryCrossEntropyMultivarCorPenalization(p); // =0 if n1==0 | n1=n
    
    value[varindex] = p;
    
    if (print_log)    Rprintf("%lf, ", value[varindex]);
  }
  if (print_log) Rprintf("\n");
  
 
 
 
 for (i = 0; i < n; i++){
   // restarts, for each i sums to each position
   if (print_log) Rprintf("\n Yprod: \n");
   
   corindex = 0;
   for (varindex = 0; varindex < number_of_responses; varindex++){
     for (int k = (1 + varindex) ; k < number_of_responses; k++){
       if (i == 0 && k == (1 + varindex)) right_sumprodxy[corindex] = 0; // initialization
       right_sumprodxy[corindex] += (y[i][varindex]*y[i][k]);
       if (print_log) Rprintf("%lf, " , y[i][varindex]*y[i][k] );
       corindex++;
     }
   }
 }
 
 corindex = 0;
  for (varindex = 0; varindex < number_of_responses; varindex++){
    for (int k = (1 + varindex) ; k < number_of_responses; k++){
      double cov0 = sqrt( (n*n1[varindex] - n1[varindex]*n1[varindex])*(n*n1[k] - n1[k]*n1[k]) );
      if (abs(cov0) > 0.0000001)
        cor_err += pow(cor_obs[corindex] - ( n*right_sumprodxy[corindex] - n1[varindex]*n1[k] )/cov0, 2);
      else cor_err += default_cor_penalization;
      
      corindex++;
    }
  }
  
  err /= number_of_responses;
  cor_err /= n_correlations;
 
 
 
  *risk = pred_penalization*err + cor_penalization*cor_err;
  
  
  if (print_log) Rprintf("\n err = %lf \n", pred_penalization*err + cor_penalization*cor_err);

}

/*
  Split function
*/
  void
  binaryCrossEntropyMultivarCorPenalization_split(int n, double *y[], double *x, int nclass,
                                 int edge, double *improve, double *split,
                                 int *csplit,
                                 double myrisk, double *wt)
{
  int i, j, k, varindex;
  double temp, best;
  int direction = LEFT;
  int where = 0;
  
  
  double err0 = 0;
  double cor_err0 = 0;
  double left_err = 0;
  double right_err = 0;
  int left_n = 0;
  int right_n;
  int corindex;
  
  double left_split_cor_vector[4950];
  double right_split_cor_vector[4950];
  double left_sumprodxy[4950];
  double right_sumprodxy[4950];

  double left_cor_error;
  double right_cor_error;
  double cov0;
    
  int *right_n1 = malloc(number_of_responses*sizeof(long));
  int *left_n1 = malloc(number_of_responses*sizeof(long));
  
  Rprintf("\n \n %d \n \n", n);
  
  for (varindex = 0; varindex < number_of_responses; varindex++){
    right_n1[varindex] = 0;
    left_n1[varindex] = 0;
    
    for (i = 0; i < n; i++){
      if (y[i][varindex] > 0){
        right_n1[varindex]++;
      }
    }
    
    if (print_log){
      Rprintf("Var: %d = ", varindex);
      Rprintf("right_n1: %d; ", right_n1[varindex]);
      Rprintf("left_n1: %d \n", left_n1[varindex]);
    }
      err0 += negBinaryCrossEntropyMultivarCorPenalization(right_n1[varindex]/((double)n)); // =0 if n1==0 | n1=n
  }
  
  for (i = 0; i < n; i++){
    Rprintf("\n");
    // restarts, for each i sums to each position
    corindex = 0;
    for (varindex = 0; varindex < number_of_responses; varindex++){
      for (int k = (1 + varindex) ; k < number_of_responses; k++){
        if (i == 0 && k == (1 + varindex)) right_sumprodxy[corindex] = 0; // initialization
        //if (print_log) Rprintf("%lf, ", y[i][varindex]*y[i][k]);
        right_sumprodxy[corindex] += y[i][varindex]*y[i][k];
        corindex++;
      }
    }
  }
  
  corindex = 0;
  cor_err0 = 0;
  for (varindex = 0; varindex < number_of_responses; varindex++){
    for (int k = (1 + varindex) ; k < number_of_responses; k++){
      cov0 = sqrt( ( (((double)n)*right_n1[varindex]) - (right_n1[varindex]*right_n1[varindex]))*((((double)n)*right_n1[k]) - (right_n1[k]*right_n1[k])) );
      if (abs(cov0) > 0.0000001){
        double aux_err0 = cor_obs[corindex] - (( ((double)n)*right_sumprodxy[corindex] - (double)(right_n1[varindex]*right_n1[k]) )/cov0);
        cor_err0 += aux_err0*aux_err0;
      }
      else cor_err0 += default_cor_penalization;
      
      
      if (print_log){
        Rprintf("\n Cor pair=%d \n" , corindex);
        Rprintf("\n cov0=%lf \n" , cov0);
        Rprintf("cor_err0=%lf \n" , cor_err0);
        Rprintf("cor_obs=%lf \n" , cor_obs[corindex]);
        Rprintf("right_sumprodxy=%lf \n" , right_sumprodxy[corindex]);
        Rprintf("rn1*rn1=%d",(right_n1[varindex]*right_n1[k]));
      }
      
      corindex++;
    }
  }
  
  err0 /= number_of_responses;
  cor_err0 /= n_correlations;
  
  if (print_log){
  Rprintf("\n err0=%lf \n" , err0);
  Rprintf("cor_err0=%lf \n" , cor_err0);
}
  if (nclass == 0) {/* continuous predictor */
      temp = 0;
      best = 0;
      
      right_n = n;
      
      for (i = 0;  right_n > edge; i++) {
        right_err = 0; 
        left_err = 0;
        left_n++;
        right_n--;
        
        if (print_log){
          
           Rprintf("i: %d = \n", i);
          
          /*
          Rprintf("\n \n Data_Left: \n");
          for (int bb = 0; bb <= i; bb++){
            for (int aa = 0;  aa< number_of_responses; aa++) {
              Rprintf("%lf," ,y[bb][aa]);
            }
            Rprintf("\n");
          }
           */
        }
           
        

        for (varindex = 0; varindex < number_of_responses; varindex++){
          
          if (y[i][varindex] > 0){
            left_n1[varindex]++;
            right_n1[varindex]--;
          }
          
          if (print_log){
            Rprintf("Var: %d = ", varindex);
            Rprintf("right_n1: %d; ", right_n1[varindex]);
            Rprintf("left_n1: %d \n", left_n1[varindex]);
            }
          
            if (x[i + 1] != x[i] && left_n >= edge) {
              left_err += ((left_n/(double)n) * negBinaryCrossEntropyMultivarCorPenalization( left_n1[varindex]/((double)left_n) ));
              right_err += ((right_n/(double)n) * negBinaryCrossEntropyMultivarCorPenalization( right_n1[varindex]/((double)right_n) ));
            }
        }
        
        
        
        //Rprintf("\n cor (left, right):\n");
        corindex = 0; // restarts, for each i sums to each position
        
        for (varindex = 0; varindex < number_of_responses; varindex++){
          for (int k = (1 + varindex) ; k < number_of_responses; k++){
            if (i == 0 && k == (1 + varindex)) left_sumprodxy[corindex] = 0; // initialization
            left_sumprodxy[corindex] += y[i][varindex]*y[i][k];
            right_sumprodxy[corindex] -= y[i][varindex]*y[i][k];
            corindex++;
          }
        }
        
        /*
          Rprintf("left_err: %lf ; ", left_err);
        Rprintf("right_err: %lf \n", right_err);
        */
          
          if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
            
            corindex = 0; // restarts, for each i sums to each position
            right_cor_error = 0;
            left_cor_error = 0;
            for (varindex = 0; varindex < number_of_responses; varindex++){
              for (int k = (1 + varindex) ; k < number_of_responses; k++){
                
                double left_cov = sqrt( (left_n*left_n1[varindex] - left_n1[varindex]*left_n1[varindex])*(left_n*left_n1[k] - left_n1[k]*left_n1[k]) );
                
                if (abs(left_cov) > 0.0000001)
                  left_split_cor_vector[corindex] = ( left_n*left_sumprodxy[corindex] - left_n1[varindex]*left_n1[k] )/left_cov;
                else left_split_cor_vector[corindex] = -2;
                
                double right_cov = sqrt( (right_n*right_n1[varindex] - right_n1[varindex]*right_n1[varindex])*(right_n*right_n1[k] - right_n1[k]*right_n1[k]) );
                
                if (abs(right_cov) > 0.0000001)
                  right_split_cor_vector[corindex] = ( right_n*right_sumprodxy[corindex] - right_n1[varindex]*right_n1[k] )/right_cov;
                else right_split_cor_vector[corindex] = -2;
                
                //Rprintf("(%lf, %lf) \n" , left_split_cor_vector[corindex],  right_split_cor_vector[corindex]);
                
                //measure the error in correlations:
                if (left_split_cor_vector[corindex] != -2){
                  left_split_cor_vector[corindex] = pow(cor_obs[corindex] - left_split_cor_vector[corindex], 2);
                } else left_split_cor_vector[corindex] = default_cor_penalization; // standard penalization
                if (right_split_cor_vector[corindex] != -2){ 
                  right_split_cor_vector[corindex] = pow(cor_obs[corindex] - right_split_cor_vector[corindex], 2);
                } else right_split_cor_vector[corindex] = default_cor_penalization;
                
                left_cor_error += left_split_cor_vector[corindex]; // todo yes can be done more efficiently
                right_cor_error += right_split_cor_vector[corindex];
                
                corindex++;
              }
            }
            
            // corrections:
            right_cor_error /= n_correlations;
            left_cor_error /= n_correlations;
            left_err /= number_of_responses;
            right_err /= number_of_responses;
            
            if (print_log){
            Rprintf("\n left_err: %lf; left_cor_error: %lf; right_err: %lf; right_cor_error: %lf \n", 
                    left_err, left_cor_error, right_err, right_cor_error);
            
            Rprintf("err0: %lf; cor_err0: %lf", err0, cor_err0);
            }
            temp = pred_penalization*err0 + cor_penalization*cor_err0 - (pred_penalization*left_err + cor_penalization*left_cor_error) -
                          (pred_penalization*right_err + cor_penalization*right_cor_error); // 
            // >0 means improvement, i.e. LESS error in split
            if (print_log){
            Rprintf("\n temp: %lf \n", temp);
              if (print_log) Rprintf(" \n DEF: %lf \n \n", default_cor_penalization);
            }
              
            
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
  binaryCrossEntropyMultivarCorPenalization_pred(double *y, double *yhat)
{
  Rprintf("Detected xpred");
  double misc = 0; //abs(y[i][0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return misc;
}













