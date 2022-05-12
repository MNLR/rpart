
#include <math.h>
#include "rpart.h"
#include "rpartproto.h"
#include <stdbool.h>


static int *countn;
static int *tsplit;
static int number_of_responses;
//static int number_of_probabilities;
static int number_of_ys;
static double informed_frequencies[256]; // now much higher than needed (fix this)
static double cond_informed_frequencies[256]; // now much higher than needed (fix this)
static double informed_frequencies_LEFT[256];
static double cond_informed_frequencies_LEFT[256];
//static int frequency_vector[256];
//static int frequency_vector_LEFT[256];
static double default_uninformative_value = 0.5;
static int marginal_frequency_LEFT = 0;
static int marginal_frequency_RIGHT = 0;
static bool log_debug = false;
static int starting_n = 0;


int binaryMargEntropyCond_init(int n, double *y[], int maxcat, char **error,
                                double *parm, int *size, int who, double *wt)
{
  if (who == 1 && maxcat > 0) { // this code needed for categorical predictors
    *error = _("Categorical predictors detected. This is not yet implemented");
  }
  
  if (log_debug) Rprintf("init called... \n");
  
  number_of_ys = (int)(*parm);

  if (number_of_ys > 5) *error = _("More than 5 responses not allowed. 5 responses already require to inform 31 probabilities");

  number_of_responses = 0;
  for (int i = 0; i < number_of_ys; i++) number_of_responses += (int) pow((double)2, (double)i);

  *size = number_of_responses;
  //number_of_probabilities = (int)pow((double)2, (double)number_of_ys);
  
  starting_n = n;
  
  if (log_debug) Rprintf("Training model for %d variables with %d responses and %d instances...\n",
      number_of_ys, 
      number_of_responses, 
      starting_n);
  //number_of_probabilities);

  return 0;
}

double binaryEntropy(double p){
  if (p == 0 || p == 1) return 0;
  return -p*log(p) - (1-p)*log(1-p);
}



double binaryMargEntropyCondAndInformP(double **y, int n){
  int idx = 0;
  int i, varindex;
  double aux;
  double entropy = 0;
  int pos0; //, idx_entropy;
  // must be set at the start, for each split candidate (launched once on _split):
  marginal_frequency_RIGHT = 0;
  marginal_frequency_LEFT = 0;
  

  //for (i = 0; i < number_of_probabilities; i++) frequency_vector[i] = 0;
  for (i = 0; i < number_of_responses; i++) {
    informed_frequencies[i] = 0;
    cond_informed_frequencies[i] = 0;
  }

  for (i = 0; i < n; i++) {
    //idx_entropy = 0;

    for (varindex = 0; varindex < number_of_ys; varindex++){
        //idx_entropy += (int)pow((double)2, (double)varindex)*y[i][varindex];

        // for the conditional probabilities:
          if (varindex == 0){
            informed_frequencies[varindex] += y[i][varindex]; //
            cond_informed_frequencies[varindex] ++;
          } else {
            idx = 0;
            pos0 = 0;

            for (int k = 0; k < varindex; k++){
              pos0 += (int)pow((double)2, (double)k);
              idx += (int)pow((double)2, (double)k)*y[i][k];
            }

            informed_frequencies[pos0 + idx] += y[i][varindex];
            cond_informed_frequencies[pos0 + idx] ++; // (denominator)
          }
    }
    
    marginal_frequency_RIGHT += y[i][number_of_ys - 1];  // for the entropy 

    //frequency_vector[idx_entropy]++;
  }

  // compute entropy (Just for the last variable):
  // the parents are not considered
  if (log_debug) Rprintf("marginal_frequency_RIGHT: %d/%d, ", marginal_frequency_RIGHT, n);

  return binaryEntropy( ((double)marginal_frequency_RIGHT)/n );
}


double addtobinaryMargEntropyCondAndReInformP(double *obs, int newn){
  int idx = 0;
  int i, varindex;
  double aux;
  //double entropy = 0;
  int pos0;//, idx_entropy;

  for (varindex = 0; varindex < number_of_ys; varindex++){
    // for the entropy:
      //idx_entropy += (int)pow((double)2, (double)varindex)*obs[varindex];

      // for the conditional probabilities:
        if (varindex == 0){
          informed_frequencies_LEFT[varindex] += obs[varindex]; //
            cond_informed_frequencies_LEFT[varindex] ++;
        } else {
          idx = 0;
          pos0 = 0;

          for (int k = 0; k < varindex; k++){
            pos0 += (int)pow((double)2, (double)k);
            idx += (int)pow((double)2, (double)k)*obs[k];
          }

          informed_frequencies_LEFT[pos0 + idx] += obs[varindex];
          cond_informed_frequencies_LEFT[pos0 + idx] ++; // so that probabilities can be computed
        }
  }
  
  marginal_frequency_LEFT += obs[number_of_ys-1];

  //frequency_vector_LEFT[idx_entropy]++;

  // compute entropy:
  
  if (log_debug) Rprintf("marginal_frequency_LEFT: %d/%d", marginal_frequency_LEFT, newn);

  return binaryEntropy( ((double)marginal_frequency_LEFT)/newn );
}

double substracttobinaryMargEntropyCondAndReInformP(double *obs, int newn){
  int idx = 0;
  int i, varindex;
  double aux;
  //double entropy = 0;
  int pos0;//, idx_entropy;

  //idx_entropy = 0;

  for (varindex = 0; varindex < number_of_ys; varindex++){
    // for the entropy:
      //idx_entropy += (int)pow((double)2, (double)varindex)*obs[varindex];

      // for the conditional probabilities:
        if (varindex == 0){
          informed_frequencies[varindex] -= obs[varindex]; //
            cond_informed_frequencies[0] --;
        } else {
          idx = 0;
          pos0 = 0;

          for (int k = 0; k < varindex; k++){
            pos0 += (int)pow((double)2, (double)k);
            idx += (int)pow((double)2, (double)k)*obs[k];
          }

          informed_frequencies[pos0 + idx] -= obs[varindex];
          cond_informed_frequencies[pos0 + idx] --; // so that probabilities can be computed
        }
  }

  marginal_frequency_RIGHT -= obs[number_of_ys-1];
  
  //frequency_vector[idx_entropy]--;
  

  // compute entropy for last variable:

  if (log_debug) Rprintf(" - marginal_frequency_RIGHT: %d/%d, ", marginal_frequency_RIGHT, newn);
  

  return  binaryEntropy( ((double)marginal_frequency_RIGHT)/newn );
}









void
binaryMargEntropyCond_eval(int n,
                            double *y[], double *value, double *risk, double *wt)
{
  // computes entropy and probabilities;
  // n/starting_n, since otherwise the entropy is almost always higher without division
  *risk = binaryMargEntropyCondAndInformP(y,n)*((double)n/starting_n); 
  if (log_debug) Rprintf("\n \n EVAL, node with p: (%d elements) \n", n);

  if (log_debug) Rprintf("\n p = \n");
  for (int i = 0; i < number_of_responses; i++){
    if (cond_informed_frequencies[i] == 0) value[i] = default_uninformative_value;
    else{
      value[i] = informed_frequencies[i]/((double)cond_informed_frequencies[i]);
      if (log_debug) Rprintf(" %lf,", value[i]);
    }
  }
  if (log_debug) Rprintf("\n Risk: %lf \n", *risk);
}

/*
  Split function
*/
  void
binaryMargEntropyCond_split(int n, double *y[], double *x, int nclass,
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

  err0 = binaryMargEntropyCondAndInformP(y, n);

  //initializes left vectors to be able to add:
  //for (i = 0; i < number_of_probabilities; i++) frequency_vector_LEFT[i] = 0;
  for (i = 0; i < number_of_responses; i++) informed_frequencies_LEFT[i] = 0;

  if (log_debug) Rprintf("\n err0=%lf \n" , err0);

  if (nclass == 0) {/* continuous predictor */
      temp = 0;
      best = 0;

      for (i = 0; right_n > edge; i++) {
        left_n++;
        right_n--;

        left_err =  addtobinaryMargEntropyCondAndReInformP(y[i], left_n);
        right_err = substracttobinaryMargEntropyCondAndReInformP(y[i], right_n);

        if (log_debug)
          Rprintf("\n %lf - ((%d/%d)%lf + (%d/%d)%lf) = ",
                  err0 , 
                  left_n, n, left_err, 
                  right_n, n, right_err);

        if (x[i + 1] != x[i] && left_n >= edge) { // repeated since has to enter for.
          temp = err0 - ( left_err*(left_n/((double)n)) + right_err*(right_n/((double)n)) );  // >0 means improvement, i.e. LESS entropy
          //temp = 2*err0 - ( left_err + right_err );  // >0 means improvement, i.e. LESS entropy
          if (temp > best) { // if its the same improvement chooses the biggest leaf
            if (log_debug) Rprintf("%lf ...", temp);
            best = temp;
            where = i;
            if (left_n < right_n)	direction = LEFT;
            else direction = RIGHT;
          } else {if (log_debug) Rprintf("%lf (No Improvement)", temp);}
        }
        if (log_debug)  Rprintf("\n \n ");
      }


      *improve = best;
      if (best > 0) {         /* found something */
          if (log_debug) Rprintf("BEST Improvement %lf ... in %d \n", best, where);

        csplit[0] = direction;
        *split = (x[where] + x[where + 1]) / 2;
      }

  } else { // categorical
    error("Categorical predictors detected in binaryCrossEntropy");
  }

}

double
binaryMargEntropyCond_pred(double *y, double *yhat)
{
  Rprintf("Detected xpred");
  double misc = 0; //abs(y[i][0] - *yhat);  // y[0] is obs, *yhat is first value of prediction

  return misc;
}
