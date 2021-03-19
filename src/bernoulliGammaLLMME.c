/*
  * The four routines for bernoulli-gamma Log Likelihood Moments Matching 
  * splitting routine
*/
  
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  
  static int *countn;
static int *tsplit;



double bernoulliGammaLogLikelihood(double y, double p, double shape, double rate){
  double density_bergamma;
  
  if (y == 0) density_bergamma = log(1-p);
  else density_bergamma = log(p) + shape*log(rate) + (shape-1)*log(y) - rate*y - lgamma(shape) ;
  
  return density_bergamma;
}


int bernoulliGammaLLMME_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaLLMME_init: Categorical predictors detected. This is not yet implemented for gammaLLMME");
  }
  
  *size = 4;  // First value is alpha/beta. 2nd: alpha; 3rd: beta
  return 0;
}


/*
  * The gamma Log Likelihood Moments Matching evaluation function.
* Negative Log Likelihood -> Lower is better (can be negative and positive)
* assigns:
  * *risk = Negative Log Likelihood
* value[0] = mean = shape/rate 
* value[1] = shape 
* value[2] = rate
*/
  
  void
  bernoulliGammaLLMME_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i, contn;
  double p, meany, sdy, shape, rate, bergamma_negloglikelihood;

  bergamma_negloglikelihood = 0;
  contn = 0;
  meany = 0.;
  for (i = 0; i < n; i++) {
    if (*y[i] > 0){
      meany += *y[i];
      contn += 1;
    }
  }
  
  if (contn != 0){
    meany = meany/contn;
    p = (float) contn / (float) n;
    if (contn > 1) {
      sdy = 0;
      for (i = 0; i < n; i++) {
        if (*y[i] > 0){
          sdy += (*y[i] - meany) * (*y[i] - meany);
        }
      }
      sdy /= contn-1; 
      if (sdy == 0) sdy = meany; // degenerated case: all elements equal
      shape = (meany*meany)/sdy;
      rate = shape/meany;
    } else { // degenerated case: only 1 continuous
      shape = meany;
      rate = 1;
    }
    
    for (i = 0; i < n; i++)
      bergamma_negloglikelihood -= bernoulliGammaLogLikelihood(*y[i], p, shape, rate);
    
  } else { // degenerated case: No continuous data
    p = 0;
    shape = -1;
    rate = -1;
  }
  
  /*
  Rprintf("Dev bergamma_negloglikelihood: ");
  Rprintf("%lf", bergamma_negloglikelihood);
  Rprintf("\n");  
  Rprintf("p: ");
  Rprintf("%lf", p);
  Rprintf("\n");  
  Rprintf("meany: ");
  Rprintf("%lf", meany);
  Rprintf("\n");  
  Rprintf("sdy: ");
  Rprintf("%lf", sdy);
  Rprintf("\n");  
  Rprintf("contn: ");
  Rprintf("%d", contn);
  Rprintf("\n");    
  */  
  value[0] = p*meany;
  value[1] = p;
  value[2] = shape;
  value[3] = rate;
  *risk = bergamma_negloglikelihood;
}

/*
  Bernoulli-Gamma splitting function
*/
  void
  bernoulliGammaLLMME_split(int n, double *y[], double *x, int nclass,
                 int edge, double *improve, double *split, int *csplit,
                 double myrisk, double *wt)
{
  int i, j, k;
  double temp, best;
  double left_sum, right_sum;
  int direction = LEFT;
  int where = 0;
  
  float left_n, right_n, right_contn, left_contn;
  double meany0, sdy0, shape0, rate0, bergamma_loglikelihood0, p0, alpha, beta;
  double left_mean, left_sd, right_mean, right_sd, left_bergamma_loglikelihood;
  double right_bergamma_loglikelihood, left_p, right_p;
  
  double degenerated_beta = 1.0f;
  
  right_sum = 0.0f;
  right_contn = 0.0f;
  
  for (i = 0; i < n; i++) {
    if (*y[i] > 0){
      right_sum += *y[i];
      right_contn += 1.0f;
    }
  }
  
  bergamma_loglikelihood0 = 0;
  
  if (right_contn == 0){ // everything is 0. No point in trying to split this further (shouldnt get here)
    *improve = 0;
  }
  else {
    // Compute log likelihood of root
    meany0 = right_sum/right_contn;
    sdy0 = 0;
    for (i = 0; i < n; i++){
      if (*y[i] > 0){
        sdy0 += (*y[i] - meany0) * (*y[i] - meany0);
      }
    }
    /*
    Rprintf("sdy0: ");
    Rprintf("%lf", sdy0);
    Rprintf("\n");    
    */
    if (sdy0 == 0 || right_contn < 2){ 
      bergamma_loglikelihood0 = (n-right_contn)*log(1-(right_contn/n)) + right_contn*bernoulliGammaLogLikelihood(meany0, right_contn/n, meany0, degenerated_beta);
    } else {
      p0 =  right_contn /  n;
      for (i = 0; i < n; i++)
        bergamma_loglikelihood0 += bernoulliGammaLogLikelihood(*y[i], p0, (meany0*meany0)/sdy0, meany0/sdy0);
    }
    
    // Compute improve
    
    right_n = n;
    right_p = 0;
    left_p = 0;
    
    if (nclass == 0) {/* continuous predictor */
      temp = 0;
      left_sum = 0;           /* No data in left branch, to start, right_sum is sum(y) */
      left_n = 0.0f;
      left_contn = 0.0f;
      best = 0;
        
        for (i = 0;  right_n > edge; i++) {
          left_n += 1.0f;
          right_n -= 1.0f;
          left_sum += *y[i];
          right_sum -= *y[i];
          if (*y[i] > 0){
            left_contn += 1.0f;
            right_contn -= 1.0f;
          }
          /*
          Rprintf("left_n: ");
          Rprintf("%f", left_n);
          Rprintf(", left_contn: ");
          Rprintf("%f", left_contn);
          Rprintf("; right_n: ");
          Rprintf("%f", right_n);
          Rprintf(", right_contn: ");
          Rprintf("%f", right_contn);
          Rprintf("\n");  
          Rprintf("Y: ");
          for (k = 0;  k <= i; k++) {
            Rprintf("%lf", *y[k]);
            Rprintf(",");
          }
          Rprintf("\n");
          */
          
          if (x[i + 1] != x[i] && left_n >= edge) {
            left_bergamma_loglikelihood = 0;
            if (left_contn > 0){
              left_p =  left_contn /  left_n;
              left_mean = left_sum/ left_contn;
              left_sd = 0;
              if (left_contn >= 2) for (k = 0; k <= i; k++) if (*y[k] > 0) left_sd += (*y[k] - left_mean) * (*y[k] - left_mean);
              if (left_sd != 0) {
                left_sd = left_sd/(left_contn - 1); // actually the quasivariance
                alpha = (left_mean*left_mean)/left_sd;
                beta = left_mean/left_sd;
                for (k = 0; k <= i; k++)
                  left_bergamma_loglikelihood += bernoulliGammaLogLikelihood(*y[k], left_p, alpha, beta);
                
              } else{
                if (left_p != 1) left_bergamma_loglikelihood = (left_n - left_contn)*log(1-left_p) + left_contn*bernoulliGammaLogLikelihood(left_mean, left_p, left_mean, degenerated_beta);
                else left_bergamma_loglikelihood = left_contn*bernoulliGammaLogLikelihood(left_mean, left_p, left_mean, degenerated_beta);
              }
            }
            
            right_bergamma_loglikelihood = 0;
            if (right_contn > 0){
              right_p =  right_contn /  right_n;
              right_mean = right_sum /  right_contn;
              right_sd = 0;
              if (right_contn >= 2) for (k = n-1; k > i; k--) if (*y[k] > 0) right_sd += (*y[k] - right_mean) * (*y[k] - right_mean);
              if (right_sd != 0) {
                right_sd = right_sd/ (right_contn - 1); // actually the quasivariance
                alpha = (right_mean*right_mean)/right_sd;
                beta = right_mean/right_sd;
                for (k = n-1; k > i; k--) right_bergamma_loglikelihood += bernoulliGammaLogLikelihood(*y[k], right_p, alpha, beta);
              } else{
                if (right_p != 1) right_bergamma_loglikelihood = (right_n - right_contn)*log(1-right_p) + right_contn*bernoulliGammaLogLikelihood(right_mean, right_p, right_mean, degenerated_beta);
                else right_bergamma_loglikelihood = right_contn*bernoulliGammaLogLikelihood(right_mean, right_p, right_mean, degenerated_beta);
              }
            } 
            /*
            Rprintf("right_sd: ");
            Rprintf("%lf", right_sd);
            Rprintf("\n");
            Rprintf("left_sd: ");
            Rprintf("%lf", left_sd);
            Rprintf("\n");
            Rprintf("right_p: ");
            Rprintf("%lf", right_p);
            Rprintf("\n");
            Rprintf("left_p: ");
            Rprintf("%lf", left_p);
            Rprintf("\n");
            Rprintf("right_mean: ");
            Rprintf("%lf", right_mean);
            Rprintf("\n");
            Rprintf("left_mean: ");
            Rprintf("%lf", left_mean);
            Rprintf("\n");
            Rprintf("left_bergamma_loglikelihood: ");
            Rprintf("%lf", left_bergamma_loglikelihood);
            Rprintf("\n");  
            Rprintf("right_bergamma_loglikelihood: ");
            Rprintf("%lf", right_bergamma_loglikelihood);
            Rprintf("\n");    
            Rprintf("bergamma_loglikelihood0: ");
            Rprintf("%lf", bergamma_loglikelihood0);
            Rprintf("\n");    
            */
            
            temp = left_bergamma_loglikelihood + right_bergamma_loglikelihood - bergamma_loglikelihood0;
            if (temp > best) {
              best = temp;
              where = i;
              if (left_sum < right_sum)	direction = LEFT; 
              else direction = RIGHT;
            }
          }
        }
        
        /*
        Rprintf("bergamma_loglikelihood0: ");
        Rprintf("%lf", bergamma_loglikelihood0);
        Rprintf("\n");    
        Rprintf("right_bergamma_loglikelihood: ");
        Rprintf("%lf", right_bergamma_loglikelihood);
        Rprintf("\n");    
        Rprintf("left_bergamma_loglikelihood: ");
        Rprintf("%lf", left_bergamma_loglikelihood);
        Rprintf("\n");    
        */
        
        *improve = best;
        if (best > 0) {         /* found something */
          csplit[0] = direction;
          *split = (x[where] + x[where + 1]) / 2;
        }
          
      }
  
  /*
    * Categorical predictor
  */
    else {
      error("In bernoulliGammaLLMME_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
    }
}
}

double
  bernoulliGammaLLMME_pred(double *y, double *yhat)
{
  Rprintf("Predicting out-of-sample... \n");
  double meanpred = fabs(y[0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}

         

