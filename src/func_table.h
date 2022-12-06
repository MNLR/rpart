/*
 * The table of implimented splitting functions
 *
 *  init_split   - Will be called before a tree is started.  May do very
 *                  little, but empirical Bayes like methods will need to set
 *                  some global variables.
 *  choose_split - function to find the best split
 *  eval         - function to calculate the response estimate and risk
 *  error        - Function that returns the prediction error.
 *  num_y        - Number of columns needed to represent y (usually 1)
 */

extern int anovainit(int n, double *y[], int maxcat, char **error,
		     double *parm, int *size, int who, double *wt);
extern int poissoninit(int n, double *y[], int maxcat, char **error,
		       double *parm, int *size, int who, double *wt);
extern int giniinit(int n, double *y[], int maxcat, char **error,
		    double *parm, int *size, int who, double *wt);
extern int usersplit_init(int n, double *y[], int maxcat, char **error,
			  double *parm, int *size, int who, double *wt);

extern void anovass(int n, double *y[], double *value, double *risk,
		    double *wt);
extern void poissondev(int n, double *y[], double *value, double *risk,
		       double *wt);
extern void ginidev(int n, double *y[], double *value, double *risk,
		    double *wt);
extern void usersplit_eval(int n, double *y[], double *value, double *risk,
			   double *wt);

extern void anova(int n, double *y[], double *x, int nclass,
		  int edge, double *improve, double *split, int *csplit,
		  double myrisk, double *wt);
extern void poisson(int n, double *y[], double *x, int nclass,
		    int edge, double *improve, double *split, int *csplit,
		    double myrisk, double *wt);
extern void gini(int n, double *y[], double *x, int nclass,
		 int edge, double *improve, double *split, int *csplit,
		 double myrisk, double *wt);
extern void usersplit(int n, double *y[], double *x, int nclass,
		      int edge, double *improve, double *split, int *csplit,
		      double myrisk, double *wt);

extern double anovapred(double *y, double *yhat);
extern double ginipred(double *y, double *yhat);
extern double poissonpred(double *y, double *yhat);
extern double usersplit_pred(double *y, double *yhat);

extern void gammaLLMME_split(int n, double *y[], double *x, int nclass,
                  int edge, double *improve, double *split, int *csplit,
                  double myrisk, double *wt);
extern void gammaLLMME_eval(int n, double *y[], double *value, double *risk,
                    double *wt);
extern int gammaLLMME_init(int n, double *y[], int maxcat, char **error,
                     double *parm, int *size, int who, double *wt);
extern double gammaLLMME_pred(double *y, double *yhat);

extern void gammaLLmean_split(int n, double *y[], double *x, int nclass,
                             int edge, double *improve, double *split, int *csplit,
                             double myrisk, double *wt);
extern void gammaLLmean_eval(int n, double *y[], double *value, double *risk,
                            double *wt);
extern int gammaLLmean_init(int n, double *y[], int maxcat, char **error,
                           double *parm, int *size, int who, double *wt);
extern double gammaLLmean_pred(double *y, double *yhat);

extern int bernoulliGammaLLMME_init(int n, double *y[], int maxcat, char **error,
                             double *parm, int *size, int who, double *wt);
extern void bernoulliGammaLLMME_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void bernoulliGammaLLMME_split(int n, double *y[], double *x, int nclass,
                            int edge, double *improve, double *split, int *csplit,
                            double myrisk, double *wt);
extern double  bernoulliGammaLLMME_pred(double *y, double *yhat);

extern int gammaDeviation_init(int n, double *y[], int maxcat, char **error,
                                    double *parm, int *size, int who, double *wt);
extern void gammaDeviation_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void gammaDeviation_split(int n, double *y[], double *x, int nclass,
                                      int edge, double *improve, double *split, int *csplit,
                                      double myrisk, double *wt);
extern double  gammaDeviation_pred(double *y, double *yhat);


extern int gammaLLBC3_init(int n, double *y[], int maxcat, char **error,
                           double *parm, int *size, int who, double *wt);
extern void gammaLLBC3_eval(int n, double *y[], double *value, double *risk,
                            double *wt);

extern void gammaLLBC3_split(int n, double *y[], double *x, int nclass,
                             int edge, double *improve, double *split, int *csplit,
                             double myrisk, double *wt);

extern double gammaLLBC3_pred(double *y, double *yhat);


extern int bernoulliLL_init(int n, double *y[], int maxcat, char **error,
                            double *parm, int *size, int who, double *wt);

extern void  bernoulliLL_eval(int n, double *y[], double *value, double *risk, double *wt);
  
extern void  bernoulliLL_split(int n, double *y[], double *x, int nclass,
                               int edge, double *improve, double *split, int *csplit,
                               double myrisk, double *wt);

extern double bernoulliLL_pred(double *y, double *yhat);

extern int binaryCrossEntropyGammaDeviation_init(int n, double *y[], int maxcat, char **error,
                               double *parm, int *size, int who, double *wt);
extern void binaryCrossEntropyGammaDeviation_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryCrossEntropyGammaDeviation_split(int n, double *y[], double *x, int nclass,
                                 int edge, double *improve, double *split, int *csplit,
                                 double myrisk, double *wt);
extern double  binaryCrossEntropyGammaDeviation_pred(double *y, double *yhat);

extern int binaryCrossEntropy_init(int n, double *y[], int maxcat, char **error,
                                   double *parm, int *size, int who, double *wt);
extern void binaryCrossEntropy_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryCrossEntropy_split(int n, double *y[], double *x, int nclass,
                                                   int edge, double *improve, double *split, int *csplit,
                                                   double myrisk, double *wt);
extern double  binaryCrossEntropy_pred(double *y, double *yhat);

extern int binaryCrossEntropyMultivar_init(int n, double *y[], int maxcat, char **error,
                                   double *parm, int *size, int who, double *wt);
extern void binaryCrossEntropyMultivar_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryCrossEntropyMultivar_split(int n, double *y[], double *x, int nclass,
                                     int edge, double *improve, double *split, int *csplit,
                                     double myrisk, double *wt);
extern double  binaryCrossEntropyMultivar_pred(double *y, double *yhat);


extern int binaryCrossEntropyMultivarCorPenalization_init(int n, double *y[], int maxcat, char **error,
                                           double *parm, int *size, int who, double *wt);
extern void binaryCrossEntropyMultivarCorPenalization_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryCrossEntropyMultivarCorPenalization_split(int n, double *y[], double *x, int nclass,
                                             int edge, double *improve, double *split, int *csplit,
                                             double myrisk, double *wt);
extern double  binaryCrossEntropyMultivarCorPenalization_pred(double *y, double *yhat);


extern int binaryMultiEntropy_init(int n, double *y[], int maxcat, char **error,
                                                          double *parm, int *size, int who, double *wt);
extern void binaryMultiEntropy_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryMultiEntropy_split(int n, double *y[], double *x, int nclass,
                                                            int edge, double *improve, double *split, int *csplit,
                                                            double myrisk, double *wt);
extern double binaryMultiEntropy_pred(double *y, double *yhat);


extern int binaryMultiEntropyCond_init(int n, double *y[], int maxcat, char **error,
                                   double *parm, int *size, int who, double *wt);
extern void binaryMultiEntropyCond_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryMultiEntropyCond_split(int n, double *y[], double *x, int nclass,
                                     int edge, double *improve, double *split, int *csplit,
                                     double myrisk, double *wt);
extern double binaryMultiEntropyCond_pred(double *y, double *yhat);


extern int binaryMargEntropyCond_init(int n, double *y[], int maxcat, char **error,
                                       double *parm, int *size, int who, double *wt);
extern void binaryMargEntropyCond_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void binaryMargEntropyCond_split(int n, double *y[], double *x, int nclass,
                                         int edge, double *improve, double *split, int *csplit,
                                         double myrisk, double *wt);
extern double binaryMargEntropyCond_pred(double *y, double *yhat);


extern int multiBinaryGammaEntropy_init(int n, double *y[], int maxcat, char **error,
                                       double *parm, int *size, int who, double *wt);
extern void multiBinaryGammaEntropy_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void multiBinaryGammaEntropy_split(int n, double *y[], double *x, int nclass,
                                         int edge, double *improve, double *split, int *csplit,
                                         double myrisk, double *wt);
extern double multiBinaryGammaEntropy_pred(double *y, double *yhat);


extern int MSEgammaDeviance_init(int n, double *y[], int maxcat, char **error,
                               double *parm, int *size, int who, double *wt);
extern void MSEgammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void MSEgammaDeviance_split(int n, double *y[], double *x, int nclass,
                                 int edge, double *improve, double *split, int *csplit,
                                 double myrisk, double *wt);
extern double  MSEgammaDeviance_pred(double *y, double *yhat);


extern int MSEbinaryEntropyGammaDeviance_init(int n, double *y[], int maxcat, char **error,
                                 double *parm, int *size, int who, double *wt);
extern void MSEbinaryEntropyGammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void MSEbinaryEntropyGammaDeviance_split(int n, double *y[], double *x, int nclass,
                                   int edge, double *improve, double *split, int *csplit,
                                   double myrisk, double *wt);
extern double  MSEbinaryEntropyGammaDeviance_pred(double *y, double *yhat);

extern int  binaryDoubleEntropyGammaDeviance_init(int n, double *y[], int maxcat, char **error,
                                              double *parm, int *size, int who, double *wt);
extern void  binaryDoubleEntropyGammaDeviance_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void  binaryDoubleEntropyGammaDeviance_split(int n, double *y[], double *x, int nclass,
                                                int edge, double *improve, double *split, int *csplit,
                                                double myrisk, double *wt);
extern double  binaryDoubleEntropyGammaDeviance_pred(double *y, double *yhat);

extern int  multiMSE_init(int n, double *y[], int maxcat, char **error,
                                                  double *parm, int *size, int who, double *wt);
extern void  multiMSE_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void  multiMSE_split(int n, double *y[], double *x, int nclass,
                                                    int edge, double *improve, double *split, int *csplit,
                                                    double myrisk, double *wt);
extern double  multiMSE_pred(double *y, double *yhat);

extern int  mae_init(int n, double *y[], int maxcat, char **error,
                          double *parm, int *size, int who, double *wt);
extern void  mae_eval(int n, double *y[], double *value, double *risk, double *wt);
extern void  mae_split(int n, double *y[], double *x, int nclass,
                            int edge, double *improve, double *split, int *csplit,
                            double myrisk, double *wt);
extern double  mae_pred(double *y, double *yhat);


static struct {
    int (*init_split) ();
    void (*choose_split) ();
    void (*eval) ();
    double (*error) ();
} func_table[] = {
    {anovainit, anova, anovass, anovapred},
    {poissoninit, poisson, poissondev, poissonpred},
    {giniinit, gini, ginidev, ginipred},
    {usersplit_init, usersplit, usersplit_eval, usersplit_pred},
    {gammaLLMME_init, gammaLLMME_split, gammaLLMME_eval, gammaLLMME_pred},
    {gammaLLmean_init, gammaLLmean_split, gammaLLmean_eval, gammaLLmean_pred},
    {bernoulliGammaLLMME_init, bernoulliGammaLLMME_split, bernoulliGammaLLMME_eval, bernoulliGammaLLMME_pred},
    {gammaDeviation_init, gammaDeviation_split, gammaDeviation_eval, gammaDeviation_pred},
    {gammaLLBC3_init, gammaLLBC3_split, gammaLLBC3_eval, gammaLLBC3_pred},
    {bernoulliLL_init, bernoulliLL_split, bernoulliLL_eval, bernoulliLL_pred},
    {binaryCrossEntropyGammaDeviation_init, binaryCrossEntropyGammaDeviation_split, binaryCrossEntropyGammaDeviation_eval, binaryCrossEntropyGammaDeviation_pred},
    {binaryCrossEntropy_init, binaryCrossEntropy_split, binaryCrossEntropy_eval, binaryCrossEntropy_pred},
    {binaryCrossEntropyMultivar_init, binaryCrossEntropyMultivar_split, binaryCrossEntropyMultivar_eval, binaryCrossEntropyMultivar_pred},
    {binaryCrossEntropyMultivarCorPenalization_init, binaryCrossEntropyMultivarCorPenalization_split, binaryCrossEntropyMultivarCorPenalization_eval, binaryCrossEntropyMultivarCorPenalization_pred},
    {binaryMultiEntropy_init, binaryMultiEntropy_split, binaryMultiEntropy_eval, binaryMultiEntropy_pred},
    {binaryMultiEntropyCond_init, binaryMultiEntropyCond_split, binaryMultiEntropyCond_eval, binaryMultiEntropyCond_pred},
    {binaryMargEntropyCond_init, binaryMargEntropyCond_split, binaryMargEntropyCond_eval, binaryMargEntropyCond_pred},
    {multiBinaryGammaEntropy_init, multiBinaryGammaEntropy_split, multiBinaryGammaEntropy_eval, multiBinaryGammaEntropy_pred},
    {MSEgammaDeviance_init, MSEgammaDeviance_split, MSEgammaDeviance_eval, MSEgammaDeviance_pred},
    {MSEbinaryEntropyGammaDeviance_init, MSEbinaryEntropyGammaDeviance_split, MSEbinaryEntropyGammaDeviance_eval, MSEbinaryEntropyGammaDeviance_pred},
    {binaryDoubleEntropyGammaDeviance_init, binaryDoubleEntropyGammaDeviance_split, binaryDoubleEntropyGammaDeviance_eval, binaryDoubleEntropyGammaDeviance_pred},
    {multiMSE_init, multiMSE_split, multiMSE_eval, multiMSE_pred},
    {mae_init, mae_split, mae_eval, mae_pred}
};

#define NUM_METHODS 23       /* size of the above structure */
