/*
  * The four routines for gamma distribution deviation splitting routine
*/
  
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
#include <stdbool.h>
  
static bool log_debug = true;
static int *countn;
static int *tsplit;

int MSEgammaDeviance_init(int n, double *y[], int maxcat, char **error,
                        double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaDeviationMSE_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
  }
  
  *size = 2;  // First value is mean_mse, second is mean_gamma
  return 0;
  if (log_debug) Rprintf("Training model with MSEgammaDeviance (%d variables) \n", size);
}



void
  MSEgammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i;
  double meany_mse, meany_gamma, gamma_dev, mse;
  
  meany_mse = 0.;
  meany_gamma = 0.;
  for (i = 0; i < n; i++) {
    meany_mse += y[i][0];
    meany_gamma += y[i][1];
  }
  
  meany_mse = meany_mse / n;
  meany_gamma = meany_gamma/n;
  
  gamma_dev = 0; 
  mse = 0;
  for (i = 0; i < n; i++) {
    mse += pow((y[i][0] - meany_mse), 2);
    gamma_dev += - log(y[i][1]/meany_gamma) + (y[i][1] - meany_gamma)/meany_gamma;
  }
  
  mse /= n;
  gamma_dev /= n;
  
  value[0] = meany_mse;
  value[1] = meany_gamma;
  
  *risk = n*(gamma_dev + pow(mse, 0.5));
  if (log_debug){
    Rprintf("###\n");
    Rprintf("Node with %d elements.\n", n);
    Rprintf("(MSE, GammaDev) = (%lf,%lf). Error = %lf \n", 
            pow(mse, 0.5), 
            gamma_dev,
            n*(gamma_dev + pow(mse, 0.5)));
  }
}

/*
  MSE+Gamma split functon
*/
  void
  MSEgammaDeviance_split(int n, double *y[], double *x, int nclass,
                     int edge, double *improve, double *split, int *csplit,
                     double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum_gamma, right_sum_gamma;
  double left_sum_mse, right_sum_mse;
  int left_n, right_n;
  int direction = LEFT;
  int where = 0;
  
  
  double meany0_mse, meany0_gamma;
  double gamma_dev0, mse0;
  double left_mean_gamma, right_mean_gamma, left_gamma_deviance, right_gamma_deviance;
  double left_mean_mse, right_mean_mse, left_mse, right_mse;
  
  
  right_sum_gamma = 0;
  right_sum_mse = 0;
  for (i = 0; i < n; i++){
    right_sum_mse += y[i][0];
    right_sum_gamma += y[i][1];
  }
  meany0_gamma = right_sum_gamma / n;
  meany0_mse = right_sum_mse / n;
  
  gamma_dev0 = 0;  
  mse0 = 0;
  for (i = 0; i < n; i++) {
    mse0 += pow((y[i][0] - meany0_mse), 2);
    gamma_dev0 += -log(y[i][1]/meany0_gamma) + (y[i][1] - meany0_gamma)/meany0_gamma;
  }
  
  gamma_dev0 /= n;
  mse0 /= n;
  right_n = n;
  
  if (nclass == 0) {/* continuous predictor */
      temp = 0;
      left_sum_gamma = 0;           /* No data in left branch, to start, right_sum is sum(y) */
      left_sum_mse = 0;
      left_n = 0;
      best = 0;
        
        for (i = 0;  right_n > edge; i++) {
          left_n++;
          right_n--;
          
          left_sum_mse += y[i][0];
          right_sum_mse -= y[i][0];
          left_sum_gamma += y[i][1];
          right_sum_gamma -= y[i][1];
          
          if (x[i + 1] != x[i] && left_n >= edge) {
            
            left_mean_gamma = left_sum_gamma/left_n;
            right_mean_gamma = right_sum_gamma/right_n;
            left_mean_mse = left_sum_mse/left_n;
            right_mean_mse = right_sum_mse/right_n;
            
            left_gamma_deviance = 0;
            left_mse = 0;
            for (k = 0; k <= i; k++){
              left_mse += pow(y[k][0] - left_mean_mse, 2);
              left_gamma_deviance += -log(y[k][1]/left_mean_gamma) + (y[k][1] - left_mean_gamma)/left_mean_gamma;
            }
            left_gamma_deviance /= left_n;
            left_mse /= left_n;
            
            right_gamma_deviance = 0;
            right_mse = 0;
            for (k = n-1; k > i; k--){ 
              right_mse +=  pow(y[k][0]- right_mean_mse, 2);
              right_gamma_deviance +=  -log(y[k][1]/right_mean_gamma) + (y[k][1] - right_mean_gamma)/right_mean_gamma;
            }
            
            right_gamma_deviance /= right_n;
            right_mse /= right_n;
            

            temp = n*(gamma_dev0 + pow(mse0, 0.5)) - left_n*(left_gamma_deviance + pow(left_mse, 0.5)) - right_n*(right_gamma_deviance + pow(right_mse, 0.5));  // >0 means improvement, i.e. LESS deviance
                
            if (log_debug){
              Rprintf("Candidate split with number of elements: (LEFT,RIGHT) = (%d,%d)\n", 
                      left_n , right_n);
                      
              Rprintf("0: (MSE, GammaDev) = %d x (%lf,%lf). Error = %lf\n", 
                      n,
                      pow(mse0, 0.5),
                      gamma_dev0,
                      n*(gamma_dev0 + pow(mse0, 0.5)));
              Rprintf("LEFT: (MSE, GammaDev) = %d x (%lf,%lf). Error = %lf\n",
                      left_n,
                      pow(left_mse, 0.5), left_gamma_deviance,
                      left_n*(left_gamma_deviance + pow(left_mse, 0.5)));
              Rprintf("RIGHT: (MSE, GammaDev) = %d x (%lf,%lf). Error = %lf\n \n",
                      right_n,
                      pow(right_mse, 0.5), right_gamma_deviance,
                      right_n*(right_gamma_deviance + pow(right_mse, 0.5)));
              Rprintf("Improvement = %lf \n\n", temp);
            }
                
                
            if (temp > best) {
              best = temp;
              where = i;
              
              if ((left_sum_gamma + left_sum_mse) < (right_sum_gamma + right_sum_mse))	direction = LEFT; 
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
      error("In gammaDeviation_MSE_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
    }
  
}

double
  MSEgammaDeviance_pred(double *y, double *yhat)
{
  Rprintf("Using gammaDeviation_MSE_pred... Only considers gamma deviance");
  double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}
