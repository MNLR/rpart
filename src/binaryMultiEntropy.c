
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"

static int *countn;
static int *tsplit;
static int number_of_responses;
static int number_of_probabilities;
static int frequency_vector[256];
static double informed_frequencies[8];
static int frequency_vector_LEFT[256];
static double informed_frequencies_LEFT[8];

int binaryMultiEntropy_init(int n, double *y[], int maxcat, char **error,
                                    double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  
  number_of_responses = (int)(*parm);
  if (number_of_responses > 8) 
    *error = _("More than 8 responses not allowed. 8 responses already require to inform 256 probabilities");
  
  
  *size = number_of_responses;
  
  number_of_probabilities = (int)pow((double)2, (double)number_of_responses);
  Rprintf("Training model with %d responses and %d probabilities...\n", number_of_responses, number_of_probabilities);
  
  return 0;
}

//todo discretization PRIOR!

double binaryMultiEntropyAndInformP(double **y, int n){
  int idx = 0;
  int i, varindex;
  double aux;
  double entropy = 0;
  
  for (i = 0; i < number_of_probabilities; i++) frequency_vector[i] = 0; 
  for (i = 0; i < number_of_responses; i++) informed_frequencies[i] = 0; 
      
  for (i = 0; i < n; i++) {
    idx = 0;
    
    for (varindex = 0; varindex < number_of_responses; varindex++){
      informed_frequencies[varindex] += y[i][varindex];
      idx += pow((double)2, (double)(number_of_responses - 1 - varindex))*y[i][varindex];
    }
    
    frequency_vector[idx]++;
  }
  
  // compute multientropy:
  Rprintf("\n Multiprobabilities: \n");
  for (i = 0; i < number_of_probabilities; i++){
    if (frequency_vector[i] != 0){
      aux = ((double) frequency_vector[i])/n;
      entropy -= aux*log(aux);
    }
    Rprintf("%d/%d, ", frequency_vector[i], n);
  }
  
  return entropy;
}


double addtoBinaryMultiEntropyAndReInformP(double *obs, int newn){
  // a single vector with an observation, this observation is ADDED
  int idx = 0;
  int i, varindex;
  double aux;
  double entropy = 0;
  
  for (varindex = 0; varindex < number_of_responses; varindex++){
    informed_frequencies_LEFT[varindex] += obs[varindex];
    idx += pow((double)2, (double)(number_of_responses - 1 - varindex))*obs[varindex];
  }
    
  frequency_vector_LEFT[idx]++;
  
  // compute ps

  // compute multientropy:
  for (i = 0; i < number_of_probabilities; i++){
    if (frequency_vector_LEFT[i] != 0){
      aux = ((double) frequency_vector_LEFT[i])/newn;
      entropy -= aux*log(aux);
    }
  }
  
  return entropy;
}

double substracttoBinaryMultiEntropyAndReInformP(double *obs, int newn){
  // a single vector with an observation, this observation is REMOVED 
  int idx = 0;
  int i, varindex;
  double aux;
  double entropy = 0;
  
  for (varindex = 0; varindex < number_of_responses; varindex++){
    informed_frequencies[varindex] -= obs[varindex];
    idx += pow((double)2, (double)(number_of_responses - 1 - varindex))*obs[varindex];
  }
  
  frequency_vector[idx]--;
  
   // compute multientropy:
  for (i = 0; i < number_of_probabilities; i++){
    if (frequency_vector[i] != 0){
      aux = ((double) frequency_vector[i])/newn;
      entropy -= aux*log(aux);
    }
  }
  
  return entropy;
}


void
binaryMultiEntropy_eval(int n, 
                        double *y[], double *value, double *risk, double *wt)
{
  Rprintf("\n \n EVAL, node with p: (%d elements) \n", n);
  *risk = binaryMultiEntropyAndInformP(y,n); // computes entropy and probabilities;
    
  for (int varindex = 0; varindex < number_of_responses; varindex++){
    value[varindex] = informed_frequencies[varindex]/((double)n);
    Rprintf("p = %lf, ", value[varindex]);
  }
  Rprintf("\n Risk: %lf \n", *risk);
  
  //Rprintf("\n");
  
    
}

/*
  Split function
*/
  void
  binaryMultiEntropy_split(int n, double *y[], double *x, int nclass,
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
  int right_n = n;
  
  err0 = binaryMultiEntropyAndInformP(y, n);
  
  //initializes left vectors to be able to add:
  for (i = 0; i < number_of_probabilities; i++) frequency_vector_LEFT[i] = 0; 
  for (i = 0; i < number_of_responses; i++) informed_frequencies_LEFT[i] = 0; 
  
  Rprintf("\n err0=%lf \n" , err0);
  
  if (nclass == 0) {/* continuous predictor */
    temp = 0;
    best = 0;
      
    for (i = 0;  right_n > edge; i++) {

      left_n++;
      right_n--;
      
      left_err =  addtoBinaryMultiEntropyAndReInformP(y[i], left_n);
      right_err = substracttoBinaryMultiEntropyAndReInformP(y[i], right_n);
      
      Rprintf("\n(%lf + %lf)" , left_err, right_err);

      if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
        temp = err0 - (left_err/left_n + right_err/right_n);  // >0 means improvement, i.e. LESS deviance
        
        if (temp > best) {
          //Rprintf("Improvement %lf ...\n", temp);
          best = temp;
          where = i;
          if (left_n < right_n)	direction = LEFT; 
          else direction = RIGHT;
        }
      }
    }
      

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
  binaryMultiEntropy_pred(double *y, double *yhat)
{
  Rprintf("Detected xpred");
  double misc = 0; //abs(y[i][0] - *yhat);  // y[0] is obs, *yhat is first value of prediction
  
  return misc;
}

