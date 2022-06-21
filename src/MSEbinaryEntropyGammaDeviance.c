/*
  * The four routines for MSE + BinaryCrossEntropy + gamma distribution deviation splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  #include <stdbool.h>
  
  static bool no_worsening = false;
  static bool log_debug = false;
  static int *countn;
  static int *tsplit;
  
  int MSEbinaryEntropyGammaDeviance_init(int n, double *y[], int maxcat, char **error,
                                         double *parm, int *size, int who, double *wt)
  {
    if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
      *error = _("In MSEbinaryEntropyGammaDeviance_init: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
    }
    
    *size = 3;  // First value is mean_mse; second is p; third is mean_gamma
    return 0;
    if (log_debug) Rprintf("Training model with MSEbinaryEntropyGammaDeviance (%d outputs) \n", size);
  }
  
  
  
  void
    MSEbinaryEntropyGammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double meany_mse, meany_gamma, gamma_dev, mse, binary_entropy, p;
    int n1;
    
    n1 = 0;
    meany_mse = 0.;
    meany_gamma = 0.;
    for (i = 0; i < n; i++) {
      meany_mse += y[i][0];
      meany_gamma += y[i][1];
      if (y[i][1] > 0) n1++;
    }
    
    meany_mse = meany_mse / n;
    meany_gamma = meany_gamma/n1;
    p = n1/(double)n;
    
    gamma_dev = 0; 
    mse = 0;
    for (i = 0; i < n; i++) {
      mse += pow((y[i][0] - meany_mse), 2);
      if (y[i][1] > 0) gamma_dev += - log(y[i][1]/meany_gamma) + (y[i][1] - meany_gamma)/meany_gamma;
    }
    
    if ( (n1 == n) || (n1 == 0) ) binary_entropy = 0.;
    else binary_entropy = -n1*log(p) - (n-n1)*log(1-p);
    
    mse /= n;
    gamma_dev /= n;
    binary_entropy /= n;
    
    value[0] = meany_mse;
    value[1] = p;
    value[2] = meany_gamma;
    
    *risk = n*(gamma_dev + pow(mse, 0.5) + binary_entropy);
    if (log_debug){
      Rprintf("###\n");
      Rprintf("Node with %d elements.\n", n);
      Rprintf("(MSE, BinEntropy, GammaDev) = (%lf, %lf, %lf). Error = %lf .\n", 
              pow(mse, 0.5), 
              binary_entropy,
              gamma_dev,
              n*(gamma_dev + binary_entropy + pow(mse, 0.5))
            );
    }
  }
  
  
  
  /*
    MSE+BE+Gamma split functon
  */
    void
    MSEbinaryEntropyGammaDeviance_split(int n, double *y[], double *x, int nclass,
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
    
    
    double meany0_mse, meany0_gamma, p0;
    double gamma_dev0, mse0, binary_entropy0;
    double left_mean_gamma, right_mean_gamma, left_gamma_deviance, right_gamma_deviance;
    double left_mean_mse, right_mean_mse, left_mse, right_mse;
    double left_p, right_p, left_binary_entropy, right_binary_entropy;
    int right_n1, left_n1;
    
    right_sum_gamma = 0;
    right_sum_mse = 0;
    right_n1 = 0;
    for (i = 0; i < n; i++){
      right_sum_mse += y[i][0];
      right_sum_gamma += y[i][1];
      if (y[i][1] > 0) right_n1++;
    }
    meany0_gamma = right_sum_gamma / right_n1;
    meany0_mse = right_sum_mse / n;
    p0 = right_n1/(double)n;
    
    gamma_dev0 = 0;  
    mse0 = 0;
    
    for (i = 0; i < n; i++) {
      mse0 += pow((y[i][0] - meany0_mse), 2);
      if (y[i][1] > 0){ gamma_dev0 +=  -log(y[i][1]/meany0_gamma) + (y[i][1] - meany0_gamma)/meany0_gamma; }
    }
    
    if ( (right_n1 == n) || (right_n1 == 0) ) binary_entropy0 = 0.;
    else binary_entropy0 = -right_n1*log(p0) - (n-right_n1)*log(1-p0);
    
    gamma_dev0 /= n;
    mse0 /= n;
    binary_entropy0 /= n;
    
    right_n = n;
    if (nclass == 0) {/* continuous predictor */
        temp = 0;
        left_sum_gamma = 0;           /* No data in left branch, to start, right_sum is sum(y) */
        left_sum_mse = 0;
        left_n1 = 0;
        left_n = 0;
        best = 0;
          
          for (i = 0;  right_n > edge; i++) {
            left_n++;
            right_n--;
            
            left_sum_mse += y[i][0];
            right_sum_mse -= y[i][0];
            left_sum_gamma += y[i][1];
            right_sum_gamma -= y[i][1];
            
            if (y[i][1] > 0){
              left_n1++;
              right_n1--;
            }
            
            if (x[i + 1] != x[i] && left_n >= edge) {
              left_mean_gamma = left_sum_gamma/left_n1;
              right_mean_gamma = right_sum_gamma/right_n1;
              left_mean_mse = left_sum_mse/left_n;
              right_mean_mse = right_sum_mse/right_n;
              left_p = left_n1/(double)left_n;
              right_p = right_n1/(double)right_n;
              
              left_gamma_deviance = 0;
              left_mse = 0;
              for (k = 0; k <= i; k++){
                left_mse += pow(y[k][0] - left_mean_mse, 2);
                if (y[k][1] > 0){ left_gamma_deviance += -log(y[k][1]/left_mean_gamma) + (y[k][1] - left_mean_gamma)/left_mean_gamma; }
              }
              
              if ( (left_n1 == left_n) || (left_n1 == 0) ) left_binary_entropy = 0.;
              else left_binary_entropy = -left_n1*log(left_p) - (left_n-left_n1)*log(1-left_p);
              
              
              left_gamma_deviance /= left_n;
              left_mse /= left_n;
              left_binary_entropy /= left_n;
              
              right_gamma_deviance = 0;
              right_mse = 0;
              for (k = n-1; k > i; k--){ 
                right_mse +=  pow(y[k][0]- right_mean_mse, 2);
                if (y[k][1] > 0){ right_gamma_deviance += -log(y[k][1]/right_mean_gamma) + (y[k][1] - right_mean_gamma)/right_mean_gamma; }
              }
              
              if ( (right_n1 == right_n) || (right_n1 == 0) ) right_binary_entropy = 0.;
              else right_binary_entropy = -right_n1*log(right_p) - (right_n-right_n1)*log(1-right_p);
              
              
              right_gamma_deviance /= right_n;
              right_mse /= right_n;
              right_binary_entropy /= right_n;
              
              temp = n*(gamma_dev0 + binary_entropy0 + pow(mse0, 0.5)) - left_n*(left_gamma_deviance + left_binary_entropy + pow(left_mse, 0.5)) - right_n*(right_gamma_deviance + right_binary_entropy + pow(right_mse, 0.5));  // >0 means improvement, i.e. LESS deviance
              
              if (log_debug){
                Rprintf("Candidate split with number of elements: (LEFT,RIGHT) = (%d,%d)\n", 
                        left_n , right_n);
                
                Rprintf("0: (MSE, BinEntropy, GammaDev) = %d x (%lf, %lf, %lf). Error = %lf.\n", 
                        n,
                        pow(mse0, 0.5),
                        binary_entropy0,
                        gamma_dev0,
                        n*(gamma_dev0 + binary_entropy0 + pow(mse0, 0.5))
                          );
                Rprintf("LEFT: (MSE, BinEntropy, GammaDev) = %d x (%lf,%lf,%lf). Error = %lf.\n",
                        left_n,
                        pow(left_mse, 0.5), 
                        left_binary_entropy,
                        left_gamma_deviance,
                        left_n*(left_gamma_deviance + left_binary_entropy + pow(left_mse, 0.5))
                  );
                Rprintf("RIGHT: (MSE, BinEntropy, GammaDev) = %d x (%lf,%lf,%lf). Error = %lf.\n",
                        right_n,
                        pow(right_mse, 0.5), 
                        right_binary_entropy,
                        right_gamma_deviance,
                        right_n*(right_gamma_deviance + right_binary_entropy + pow(right_mse, 0.5))
                  );
                Rprintf("Improvement = %lf \n\n", temp);
              }
              
              if (no_worsening){
                if ((n*(gamma_dev0 + binary_entropy0) > 
                       ( left_n*(left_gamma_deviance + left_binary_entropy) + 
                         right_n*(right_gamma_deviance + right_binary_entropy) 
                       )
                ) && ( n*pow(mse0, 0.5)  > (left_n*(pow(left_mse, 0.5)) + right_n*(pow(right_mse, 0.5))) ) 
                ) {
                  if (temp > best) {
                    best = temp;
                    where = i;
                    if (log_debug){ Rprintf("Current BEST Improvement %lf ... in %d \n \n \n", best, where + 1); }// I add 1 for R 
                    
                    if ((left_sum_gamma + left_sum_mse) < (right_sum_gamma + right_sum_mse))	direction = LEFT; 
                    else direction = RIGHT;
                  }
                }
              } else {
                if (temp > best) {
                  best = temp;
                  where = i;
                  if (log_debug){ Rprintf("Current BEST Improvement %lf ... in %d \n \n \n", best, where + 1); }// I add 1 for R 
                  
                  if ((left_sum_gamma + left_sum_mse) < (right_sum_gamma + right_sum_mse))	direction = LEFT; 
                  else direction = RIGHT;
                }
              }
          
          *improve = best;
          if (best > 0) {         /* found something */
            csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2;
          }
          
        }
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
      MSEbinaryEntropyGammaDeviance_pred(double *y, double *yhat)
  {
    Rprintf("Using gammaDeviation_MSE_pred... Only considers gamma deviance");
    double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
    
    return meanpred;
  }
  