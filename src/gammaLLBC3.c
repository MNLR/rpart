/*
  * The four routines for gamma Log Likelihood Louzada_2018 BC3 estimator splitting routine
*/
  #include <stdbool.h>
  #include <math.h>
  #include "rpart.h"
  #include "rpartproto.h"
  
  static int *countn;
  static int *tsplit;

int gammaLLBC3_init(int n, double *y[], int maxcat, char **error,
                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("In gammaLLBC3_init: Categorical predictors detected. This is not yet implemented for gammaLLMME");
  }
  
  *size = 3;  // First value is alpha/beta. 2nd: alpha; 3rd: beta
  return 0;
}



double gammaLogLikelihood2(double y, double shape, double rate){
  double density_gamma = shape*log(rate) + (shape-1)*log(y) - rate*y - lgamma(shape) ;
  return density_gamma;
}

  void
  gammaLLBC3_eval(int n, double *y[], double *value, double *risk, double *wt)
{
  // weights disabled
  int i;
  double sumy, sumlogy, aux, shape, rate, gamma_negloglikelihood;
  bool distinct_elements = false;
  double degenerated_rate = 1;
  
  sumy = 0.;
  sumlogy = 0.;
  aux = 0.;
  for (i = 0; i < n; i++) {
    if (!distinct_elements) distinct_elements = ( *y[0] != *y[i] );

    sumy += *y[i];
    sumlogy += log(*y[i]);
    aux += (*y[i])*log(*y[i]);
  }
  
  if (distinct_elements){
    aux = n*aux - sumy*sumlogy;
    shape = n*sumy/aux;
    shape = shape - ( 3*shape - (2*shape)/( 3*(1+shape) ) - (4*shape)/( 5*(1+shape)*(1+shape) ) )/n; // bias correction
    rate = n*shape/sumy;
    
    gamma_negloglikelihood = 0;
    for (i = 0; i < n; i++)
      gamma_negloglikelihood -= gammaLogLikelihood2(*y[i], shape, rate);
    
  } else {
    shape = sumy/n;
    rate = degenerated_rate;
    
    gamma_negloglikelihood = 
      -n*( shape*log(rate) + (shape-1)*log(*y[0]) - rate*(*y[0]) - lgamma(shape) );
    
  }

  // 
    
  value[0] = shape/rate;
  value[1] = shape;
  value[2] = rate;
  *risk = gamma_negloglikelihood;
  
}

/*
  Gamma spliting functon
*/
  void
  gammaLLBC3_split(int n, double *y[], double *x, int nclass,
                 int edge, double *improve, double *split, int *csplit,
                 double myrisk, double *wt)
{
  double degenerated_rate = 1; // to be parametrized
  
  bool right_distinct_elements, left_distinct_elements;
    
  int i, j, k;
  double temp, best;
double left_sum, right_sum, right_sumlog, left_sumlog, right_sumylog, left_sumylog, aux;
  int left_n, right_n;
  int direction = LEFT;
  int where = 0;
  
  
  double shape, rate, gamma_loglikelihood0, left_aux;
  double left_gamma_logLikelihood, right_gamma_logLikelihood;
  
  right_sum = 0.;
  right_sumlog = 0.;
  right_sumylog = 0.;
  right_distinct_elements = false;
  
  for (i = 0; i < n; i++){
    if (!right_distinct_elements) right_distinct_elements = *y[0] != *y[i];
    right_sum += *y[i];
    right_sumlog += log(*y[i]);
    right_sumylog += (*y[i])*log(*y[i]);
  }

  if (!right_distinct_elements){ // Either one element or all elements are the same value. No point in spliting further
    *improve = 0;
  } else {
    aux = n*right_sumylog - right_sum*right_sumlog;
    shape = n*right_sum/aux;
    shape = shape - ( 3*shape - (2*shape)/( 3*(1+shape) ) - (4*shape)/( 5*(1+shape)*(1+shape) ) )/n; // bias correction
    
    rate = n*shape/right_sum;
 
    gamma_loglikelihood0 = 0;
    for (i = 0; i < n; i++) gamma_loglikelihood0 += gammaLogLikelihood2(*y[i], shape, rate);

    right_n = n;
    
    if (nclass == 0) {/* continuous predictor */
        temp = 0;
        left_n = 0;
        best = 0;
        
        left_sum = 0.;           /* No data in left branch, to start, right_sum is sum(y) */
        left_sumlog = 0.;
        left_sumylog = 0.;
          
        left_distinct_elements = false;
          
        for (i = 0; right_n > edge; i++) {
          left_n++;
          right_n--;
          
          left_sum += *y[i];
          left_sumlog += log(*y[i]);
          left_sumylog += (*y[i])*log(*y[i]);
          if (!left_distinct_elements) left_distinct_elements = (*y[0] != *y[i]);
          
          right_sum -= *y[i];
          right_sumlog -= log(*y[i]);
          right_sumylog -= (*y[i])*log(*y[i]);


          if (x[i + 1] != x[i] && left_n >= edge) {

              if (left_distinct_elements) {
                aux = left_n*left_sumylog - left_sum*left_sumlog;
                shape = left_n*left_sum/aux;
                shape = shape - ( 3*shape - (2*shape)/( 3*(1+shape) ) - (4*shape)/( 5*(1+shape)*(1+shape) ) )/left_n; // bias correction
                
                rate = left_n*shape/left_sum;
                
                left_gamma_logLikelihood = 0;
                for (k = 0; k <= i; k++) left_gamma_logLikelihood += gammaLogLikelihood2(*y[k], shape, rate);
              } else {
                shape = left_sum/left_n;
                left_gamma_logLikelihood = 
                  left_n*( shape*log(degenerated_rate) + (shape-1)*log(*y[0]) - degenerated_rate*(*y[0]) - lgamma(shape) );
              }
              
              right_distinct_elements = false;
              
              // Computes the estimator even if not distinct
              aux = right_n*right_sumylog - right_sum*right_sumlog;
              shape = right_n*right_sum/aux;
              shape = shape - ( 3*shape - (2*shape)/( 3*(1+shape) ) - (4*shape)/( 5*(1+shape)*(1+shape) ) )/right_n; // bias correction
              
              rate = right_n*shape/right_sum;
              

              right_gamma_logLikelihood = 0;
              for (k = n-1; k > i; k--) {
                
                if (!right_distinct_elements) right_distinct_elements = (*y[n-1] != *y[k]); // checks the right leaf distinct
                right_gamma_logLikelihood += gammaLogLikelihood2(*y[k], shape, rate);
              }
              
              if (!right_distinct_elements) { // degenerated estimator. Recompute.
                shape = right_sum/right_n;  

                right_gamma_logLikelihood = 
                  right_n*( shape*log(degenerated_rate) + (shape-1)*log(*y[n-1]) - degenerated_rate*(*y[n-1]) - lgamma(shape) );
              }  
              
              temp = left_gamma_logLikelihood + right_gamma_logLikelihood - gamma_loglikelihood0;
              
              if (temp > best) {
                    best = temp;
                    
                    where = i;
                    if (left_sum < right_sum)	direction = LEFT; 
                    else direction = RIGHT;
              }
            }
          }
          
          
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
        error("In gammaLLMME_split: Categorical predictors detected. This is not yet implemented for gammaLLMME");
      }
  }
}

double
  gammaLLBC3_pred(double *y, double *yhat)
{
  Rprintf("Predicting out-of-sample... \n");
  double meanpred = fabs(y[0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return meanpred;
}
