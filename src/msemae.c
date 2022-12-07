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
  
  int msemae_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
  {
    if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
      *error = _("In msemae_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
    }
    
    *size = 2;  // 
     if (log_debug) Rprintf("Training model with MSE-MAE (%d variables) \n", *size);
    
    return 0;
  }
  
  
  
  void
  msemae_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double meany_mse, meany_mae, mse, mae;
    
    meany_mse = 0.;
    meany_mae = 0.;
    for (i = 0; i < n; i++) {
      meany_mse += y[i][0];
      meany_mae += y[i][1];
    }
    
    meany_mse = meany_mse/n;
    meany_mae = meany_mae/n;
    
    mse = 0; 
    mae = 0;
    for (i = 0; i < n; i++) {
      mse += pow((y[i][0] - meany_mse), 2);
      mae += fabs(y[i][1] - meany_mae);
    }
    
    mse /= n;
    mae /= n;
    
    value[0] = meany_mse;
    value[1] = meany_mae;
    
    *risk = n*(pow(mse, 0.5) + mae);
    if (log_debug){
      Rprintf("###\n");
      Rprintf("Node with %d elements.\n", n);
      Rprintf("(MSE, MAE) = (%lf,%lf). Error = %lf \n", 
              pow(mse, 0.5), 
              mae,
              n*(pow(mse, 0.5) + mae));
    }
  }
  
  /*
    MSE+Gamma split functon
  */
    void
  msemae_split(int n, double *y[], double *x, int nclass,
                 int edge, double *improve, double *split, int *csplit,
                 double myrisk, double *wt)
  {
    int i, j, k;
    double temp, best;
    double left_sum_mae, right_sum_mae;
    double left_sum_mse, right_sum_mse;
    int left_n, right_n;
    int direction = LEFT;
    int where = 0;
    
    
    double meany0_mse, meany0_mae;
    double mae0, mse0;
    double left_mean_mae, right_mean_mae, left_mae, right_mae;
    double left_mean_mse, right_mean_mse, left_mse, right_mse;
    
    
    right_sum_mae = 0;
    right_sum_mse = 0;
    for (i = 0; i < n; i++){
      right_sum_mse += y[i][0];
      right_sum_mae += y[i][1];
    }
    meany0_mse = right_sum_mse / n;
    meany0_mae = right_sum_mae / n;
    
    mae0 = 0;  
    mse0 = 0;
    for (i = 0; i < n; i++) {
      mse0 += pow((y[i][0] - meany0_mse), 2);
      mae0 += fabs(y[i][1] - meany0_mae);
    }
    
    mae0 /= n;
    mse0 /= n;
    right_n = n;
    
    if (nclass == 0) {/* continuous predictor */
        temp = 0;
        left_sum_mae = 0;           /* No data in left branch, to start, right_sum is sum(y) */
        left_sum_mse = 0;
        left_n = 0;
        best = 0;
          
        for (i = 0;  right_n > edge; i++) {
          left_n++;
          right_n--;
          
          left_sum_mse += y[i][0];
          right_sum_mse -= y[i][0];
          left_sum_mae += y[i][1];
          right_sum_mae -= y[i][1];
          
          if (x[i + 1] != x[i] && left_n >= edge) {
            
            
            left_mean_mse = left_sum_mse/left_n;
            right_mean_mse = right_sum_mse/right_n;
            
            left_mean_mae = left_sum_mae/left_n;
            right_mean_mae = right_sum_mae/right_n;
            
            left_mae = 0;
            left_mse = 0;
            for (k = 0; k <= i; k++){
              left_mse += pow(y[k][0] - left_mean_mse, 2);
              left_mae +=  fabs(y[k][1] - left_mean_mae);
            }
            left_mse /= left_n;
            left_mae /= left_n;
            
            right_mae = 0;
            right_mse = 0;
            for (k = n-1; k > i; k--){ 
              right_mse +=  pow(y[k][0] - right_mean_mse, 2);
              right_mae += fabs(y[k][1] - right_mean_mae); 
            }
            right_mse /= right_n;
            right_mae /= right_n;
            
            temp = n*(mae0 + pow(mse0, 0.5)) - left_n*(left_mae + pow(left_mse, 0.5)) - right_n*( right_mae + pow(right_mse, 0.5));  // >0 means improvement, i.e. LESS deviance
            
            if (log_debug){
              Rprintf("Candidate split with number of elements: (LEFT,RIGHT) = (%d,%d)\n", 
                      left_n , right_n);
              
              Rprintf("0: (RMSE, MAE) = %d x (%lf,%lf). Error = %lf\n", 
                      n,
                      pow(mse0, 0.5),
                      mae0,
                      n*(mae0 + pow(mse0, 0.5)));
              Rprintf("LEFT: (RMSE, MAE) = %d x (%lf,%lf). Error = %lf\n",
                      left_n,
                      pow(left_mse, 0.5), left_mae,
                      left_n*(left_mae + pow(left_mse, 0.5)));
              Rprintf("RIGHT: (RMSE, MAE) = %d x (%lf,%lf). Error = %lf\n \n",
                      right_n,
                      pow(right_mse, 0.5),  right_mae,
                      right_n*( right_mae + pow(right_mse, 0.5)));
              Rprintf("Improvement = %lf \n\n", temp);
            }
            
            
            if (temp > best) {
              best = temp;
              where = i;
              
              if ((left_sum_mae + left_sum_mse) < (right_sum_mae + right_sum_mse))	direction = LEFT; 
              else direction = RIGHT;
            }
          }
        }
        
          *improve = best;
          if (best > 0) {         /* found something */
              if (log_debug) Rprintf("BEST Improvement %lf ... in %d \n \n \n", best, where + 1); // I add 1 for R
            csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2;
          }
          
    }
    /*
      * Categorical predictor
    */
      else {
        error("In mae_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
      }
    
  }
  
  double
  msemae_pred(double *y, double *yhat)
  {
    Rprintf("Using mae_pred: Only first variable considered");
    double meanpred =  pow(pow( y[0] - *yhat, 2), 0.5)  ;  // y[0] is obs, *yhat is first value of prediction
    
    return meanpred;
  }
  