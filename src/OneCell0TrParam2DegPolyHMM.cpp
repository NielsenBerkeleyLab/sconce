#include "OneCell0TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
OneCell0TrParam2DegPolyHMM::OneCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy) : OneCell0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 0, 3, 0, 1) { // 0 transition params to est, 3 fixedTrParams, 0 fixedLibs, 1 branch
}
OneCell0TrParam2DegPolyHMM::OneCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : OneCell3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches) {
  this->maxNumBFGSStarts = 0;
}

OneCell0TrParam2DegPolyHMM::~OneCell0TrParam2DegPolyHMM() {
  // TODO
}

/*
 ********
 * accessors and mutators
 ********
 */
/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
/*
 * transitionParams = [t]
 */
double OneCell0TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  double beta  = gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1);
  double lambda = gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2);

  double t = gsl_vector_get(transitionParams, 0);

  // save into paramsToEst
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX, gsl_vector_get(transitionParams, 0));

  return OneCell3TrParam2DegPolyHMM::setTransition(dest, this->getAlpha(), beta, lambda, t); // call parent class's setTransition helper method
}

void OneCell0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);

  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t)); // set T
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(s));
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCell0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);

  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T)); // set t
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
}

/*
 ********
 * functions
 ********
 */

/*
 * calls bfgs, returns a bestGuessHMM
 */
OneCell0TrParam2DegPolyHMM* OneCell0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  // create new HMM with the best guess parameters and return it
  OneCell0TrParam2DegPolyHMM* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, verbose);
  return bestGuessHMM;
}

