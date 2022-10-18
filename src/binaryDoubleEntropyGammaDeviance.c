/*
  * The four routines for MSE + BinaryCrossEntropy + gamma distribution deviation splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  #include <stdbool.h>
  
  static bool no_worsening = false;
  static double worsening_tolerance = 0.01;
  static bool log_debug = false;
  static int *countn;
  static int *tsplit;
  
  int binaryDoubleEntropyGammaDeviance_init(int n, double *y[], int maxcat, char **error,
                                         double *parm, int *size, int who, double *wt)
  {
    if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
      *error = _("In binaryDoubleEntropyGammaDeviation: Categorical predictors detected. This is not yet implemented for gammaDeviation_init");
    }
    
    *size = 5;  // First value is p00; second is p01, third is p10, fourth is p11; fifth is mean_gamma
    return 0;
    if (log_debug) Rprintf("Training model with binaryDoubleEntropyGammaDeviation (%d outputs) \n", size);
  }
  
  
  
  void
    binaryDoubleEntropyGammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt)
  {
    // weights disabled
    int i;
    double meany_gamma, gamma_dev, binary_entropy, p00, p01, p10, p11;
    int n00, n01, n10, n11, n1;
    
    n00 = 0;
    n01 = 0;
    n10 = 0;
    n11 = 0;
    meany_gamma = 0.;
    for (i = 0; i < n; i++) {
      meany_gamma += y[i][1];
      if (y[i][0] > 0 && y[i][1] > 0) n11++;
      else if (y[i][0] == 0 && y[i][1] == 0) n00++;
      else if (y[i][0] == 0 && y[i][1] > 0) n01++;
      else {n10++;}
    }
    
    n1 = n01 + n11;
    
    if (n1 > 0) meany_gamma = meany_gamma/n1;

    
    p00 = n00/(double)n;
    p01 = n01/(double)n;
    p10 = n10/(double)n;
    p11 = n11/(double)n;
    
    gamma_dev = 0; 
    if (n1 > 0){
      for (i = 0; i < n; i++) {
        if (y[i][1] > 0) gamma_dev += -log(y[i][1]/meany_gamma) + (y[i][1] - meany_gamma)/meany_gamma;
      }
    }
    
    binary_entropy = 0.;
    if ( n00 == n  || n01 == n || n10 == n || n11 == n){
      // do nothing, perfectly split
    } else {
      if ( n00 != 0 ) binary_entropy -= n00*log(p00);
      if ( n01 != 0 ) binary_entropy -= n01*log(p01);
      if ( n10 != 0 ) binary_entropy -= n10*log(p10);
      if ( n11 != 0 ) binary_entropy -= n11*log(p11);
    }
      
    gamma_dev /= n;
    binary_entropy /= n;
    
    value[0] = p00;
    value[1] = p01;
    value[2] = p10;
    value[3] = p11;
    value[4] = meany_gamma;
    
    *risk = n*(gamma_dev + binary_entropy);
    if (log_debug){
      Rprintf("###\n");
      Rprintf("Node with %d elements.\n", n);
      Rprintf("Values: (P00 = %lf), (P01 = %lf), (P10 = %lf), (P11 = %lf), (Mean(Gamma) = %lf). \n", 
              p00, p01, p10, p11, meany_gamma);
      Rprintf("(BinEntropy, GammaDev) = (%lf, %lf). Error = %lf .\n #EndEVAL#\n", 
              binary_entropy,
              gamma_dev,
              n*(gamma_dev + binary_entropy)
      );
    }
  }
  
  
  
  /*
    ME+Gamma split functon
  */
    void
    binaryDoubleEntropyGammaDeviance_split(int n, double *y[], double *x, int nclass,
                                      int edge, double *improve, double *split, int *csplit,
                                      double myrisk, double *wt)
  {
    int i, j, k;
    double temp, best;
    double left_sum_gamma, right_sum_gamma;
    int left_n, right_n;
    int direction = LEFT;
    int where = 0;
    
    
    double meany0_gamma;
    double gamma_dev0, binary_entropy0;
    double left_mean_gamma, right_mean_gamma, left_gamma_deviance, right_gamma_deviance;
    double left_p00, left_p01, left_p10, left_p11, right_p00, right_p01, right_p10, right_p11, left_binary_entropy, right_binary_entropy;
    int left_n1, left_n00, left_n01, left_n10, left_n11, right_n1, right_n00, right_n01, right_n10, right_n11;
    
    right_n00 = 0;
    right_n01 = 0;
    right_n10 = 0;
    right_n11 = 0;
    
    left_n00 = 0;
    left_n01 = 0;
    left_n10 = 0;
    left_n11 = 0;
    
    right_sum_gamma = 0;
    for (i = 0; i < n; i++){
      right_sum_gamma += y[i][1];
      if (y[i][0] > 0 && y[i][1] > 0) right_n11++;
      else if (y[i][0] == 0 && y[i][1] == 0) right_n00++;
      else if (y[i][0] == 0 && y[i][1] > 0) right_n01++;
      else {right_n10++;}
    }
    right_n1 = right_n11 + right_n01;
    
    if (right_n1 > 0) meany0_gamma = right_sum_gamma / right_n1;
    
    right_p00 = right_n00/(double)n;
    right_p01 = right_n01/(double)n;
    right_p10 = right_n10/(double)n;
    right_p11 = right_n11/(double)n;
    
    gamma_dev0 = 0;  
    if (right_n1 > 0){
      for (i = 0; i < n; i++) {
        if (y[i][1] > 0){ gamma_dev0 +=  -log(y[i][1]/meany0_gamma) + (y[i][1] - meany0_gamma)/meany0_gamma; }
      }
    }
    
    
    binary_entropy0 = 0.;
    if ( right_n00 == n  || right_n01 == n || right_n10 == n || right_n11 == n){
      // do nothing, perfectly split
    } else {
      if (  right_n00 != 0 ) binary_entropy0 -=  right_n00*log(right_p00);
      if (  right_n01 != 0 ) binary_entropy0 -=  right_n01*log(right_p01);
      if (  right_n10 != 0 ) binary_entropy0 -=  right_n10*log(right_p10);
      if (  right_n11 != 0 ) binary_entropy0 -=  right_n11*log(right_p11);
    }
    
    gamma_dev0 /= n;
    binary_entropy0 /= n;
    
    right_n = n;
    if (nclass == 0) {/* continuous predictor */
        temp = 0;
        left_sum_gamma = 0;           /* No data in left branch, to start, right_sum is sum(y) */
          left_n1 = 0;
          left_n = 0;
          best = 0;
          
          for (i = 0;  right_n > edge; i++) {
            left_n++;
            right_n--;
            
            left_sum_gamma += y[i][1];
            right_sum_gamma -= y[i][1];
            
            if (y[i][0] > 0 && y[i][1] > 0){
              left_n1++;
              right_n1--;
              left_n11++;
              right_n11--;
            } else if (y[i][0] == 0 && y[i][1] == 0){
              left_n00++;
              right_n00--;
            } else if (y[i][0] == 0 && y[i][1] > 0){
              left_n1++;
              right_n1--;
              left_n01++;
              right_n01--;
            } else {
              left_n10++;
              right_n10--;
            }
            
            if (x[i + 1] != x[i] && left_n >= edge) {
              left_mean_gamma = left_sum_gamma/left_n1;
              right_mean_gamma = right_sum_gamma/right_n1;
              
              left_p00 = left_n00/(double)left_n;
              left_p01 = left_n01/(double)left_n;
              left_p10 = left_n10/(double)left_n;
              left_p11 = left_n11/(double)left_n;
              
              right_p00 = right_n00/(double)right_n;
              right_p01 = right_n01/(double)right_n;
              right_p10 = right_n10/(double)right_n;
              right_p11 = right_n11/(double)right_n;
              
              left_gamma_deviance = 0;
              if (left_n1 > 0){
                for (k = 0; k <= i; k++){
                  if (y[k][1] > 0){ left_gamma_deviance += -log(y[k][1]/left_mean_gamma) + (y[k][1] - left_mean_gamma)/left_mean_gamma; }
                }
              }
              
              left_binary_entropy = 0.;
              if ( left_n00 == left_n  || left_n01 == left_n || left_n10 == left_n || left_n11 == left_n){
                // do nothing, perfectly split
              } else {
                if (  left_n00 != 0 ) left_binary_entropy -=  left_n00*log(left_p00);
                if (  left_n01 != 0 ) left_binary_entropy -=  left_n01*log(left_p01);
                if (  left_n10 != 0 ) left_binary_entropy -=  left_n10*log(left_p10);
                if (  left_n11 != 0 ) left_binary_entropy -=  left_n11*log(left_p11);
              }
              
              left_gamma_deviance /= left_n;
              left_binary_entropy /= left_n;
              
              right_gamma_deviance = 0;
              if (right_n1){
                for (k = n-1; k > i; k--){ 
                  if (y[k][1] > 0){ right_gamma_deviance += -log(y[k][1]/right_mean_gamma) + (y[k][1] - right_mean_gamma)/right_mean_gamma; }
                }
              }
              
              right_binary_entropy = 0.;
              if ( right_n00 == right_n  || right_n01 == right_n || right_n10 == right_n || right_n11 == right_n){
                // do nothing, perfectly split
              } else {
                if (  right_n00 != 0 ) right_binary_entropy -=  right_n00*log(right_p00);
                if (  right_n01 != 0 ) right_binary_entropy -=  right_n01*log(right_p01);
                if (  right_n10 != 0 ) right_binary_entropy -=  right_n10*log(right_p10);
                if (  right_n11 != 0 ) right_binary_entropy -=  right_n11*log(right_p11);
              }
              
              
              right_gamma_deviance /= right_n;
              right_binary_entropy /= right_n;
              
              temp = n*(gamma_dev0 + binary_entropy0) - left_n*(left_gamma_deviance + left_binary_entropy) - right_n*(right_gamma_deviance + right_binary_entropy);  // >0 means improvement, i.e. LESS deviance
              
              if (log_debug){
                Rprintf("Candidate split with number of elements: (LEFT,RIGHT) = (%d,%d)\n", 
                        left_n , right_n);
                
                Rprintf("0: (BinEntropy, GammaDev) = %d x (%lf, %lf). Error = %lf.\n", 
                        n,
                        binary_entropy0,
                        gamma_dev0,
                        n*(gamma_dev0 + binary_entropy0)
                );
                Rprintf("LEFT: (BinEntropy, GammaDev) = %d x (%lf,%lf). Error = %lf.\n",
                        left_n,
                        left_binary_entropy,
                        left_gamma_deviance,
                        left_n*(left_gamma_deviance + left_binary_entropy)
                );
                Rprintf("RIGHT: (BinEntropy, GammaDev) = %d x (%lf,%lf). Error = %lf.\n",
                        right_n,
                        right_binary_entropy,
                        right_gamma_deviance,
                        right_n*(right_gamma_deviance + right_binary_entropy)
                );
                Rprintf("Improvement = %lf \n", temp);
              }
              
              if (no_worsening){
                if ( ( (n*gamma_dev0 + worsening_tolerance) >= (left_n*left_gamma_deviance + right_n*right_gamma_deviance) ) &&
                     ( (n*binary_entropy0 + worsening_tolerance) >= (left_n*left_binary_entropy + right_n*right_binary_entropy) ) 
                ) {
                  if (temp > best) {
                    best = temp;
                    where = i;
                    if (log_debug){ Rprintf("Current BEST Improvement %lf ... in %d \n \n \n", best, where + 1); }// I add 1 for R 
                    
                    if ((left_sum_gamma + left_binary_entropy) < (right_sum_gamma + right_binary_entropy))	direction = LEFT; 
                    else direction = RIGHT;
                  }
                } else {
                  if (log_debug){ Rprintf("Worse marginal results even though temp = %lf \n\n\n", temp); }
                }
              } else {
                if (temp > best) {
                  best = temp;
                  where = i;
                  if (log_debug){ Rprintf("Current BEST Improvement %lf ... in %d \n \n \n", best, where + 1); }// I add 1 for R 
                  
                  if ((left_sum_gamma + left_binary_entropy) < (right_sum_gamma + right_binary_entropy))	direction = LEFT; 
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
    binaryDoubleEntropyGammaDeviance_pred(double *y, double *yhat)
  {
    Rprintf("Using gammaDeviation_MENT_pred... Only considers gamma deviance");
    double meanpred = - log(y[0]/(*yhat)) +  (y[0]-(*yhat))/(*yhat) ;  // y[0] is obs, *yhat is first value of prediction
    
    return meanpred;
  }
  