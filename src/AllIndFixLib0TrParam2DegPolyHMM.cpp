#include "AllIndFixLib0TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllIndFixLib0TrParam2DegPolyHMM::AllIndFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy) : AllIndFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 0, 3, depths->size(), depths->size(), depths->size() * 1) { // 0 shared transition params to est, 3 fixed transition param (alpha/beta/lambda), n fixed libs, n HMMs, n * 1 branches
}
AllIndFixLib0TrParam2DegPolyHMM* AllIndFixLib0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy) {
  AllIndFixLib0TrParam2DegPolyHMM* hmm = new AllIndFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy);
  hmm->makeHMMs();
  return hmm;
}

AllIndFixLib0TrParam2DegPolyHMM::AllIndFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst) : AllInd0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numSharedTrParamsToEst, numFixedTrParams, numFixedLibs, numHMMs, numBranchesToEst) {
}

void AllIndFixLib0TrParam2DegPolyHMM::makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams) {
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  gsl_vector* currFixedParams = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->depthsVec->size(); hmmIdx++) {
    currDepths = new std::vector<DepthPair*>();
    currDepths->push_back((*this->depthsVec)[hmmIdx]);

    currFixedParams = gsl_vector_alloc(1 + 3); // 1 lib, alpha/beta/lambda
    // set libs
    gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx));

    // set shared transition params
    gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0)); // alpha
    gsl_vector_set(currFixedParams, 2, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1)); // beta
    gsl_vector_set(currFixedParams, 3, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2)); // lambda
 
    OneCellFixLib0TrParam2DegPolyHMM* hmm = new OneCellFixLib0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy());
    (*this->hmmVec)[hmmIdx] = hmm;

    // do rest of HMM set up (usually happens in main.cpp)
    currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
    gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
    hmm->setMeanVarianceFn(currMeanVarCoefVec);
    hmm->setTransition(transitionParams); // this vector isn't saved anywhere
    hmm->setAlpha(this->getAlpha());

    this->setLibScalingFactor(hmmIdx, hmm->getLibScalingFactor(0));
  }
}

AllIndFixLib0TrParam2DegPolyHMM::~AllIndFixLib0TrParam2DegPolyHMM() {
  // TODO
}

void AllIndFixLib0TrParam2DegPolyHMM::setLibScalingFactor(int cellIdx, double lib) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, lib);
  (*this->hmmVec)[cellIdx]->setLibScalingFactor(0, lib);
}

/*
 * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
 */
void AllIndFixLib0TrParam2DegPolyHMM::setAllLibScalingFactors(double libScalingFactor) {
  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
    gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, libScalingFactor);
    (*this->hmmVec)[i]->setLibScalingFactor(0, libScalingFactor);
  }
}
double AllIndFixLib0TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

// Optimizable methods
/*
 * method to extract HMM specific variables and set those variables for each HMM,
 * using each HMM's own setParamsToEst method.
 * params is in probability space
 */
double AllIndFixLib0TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);

  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);
  double currTBFGS = 0;

  double status = 0;
  // for each HMM, call that HMM's setParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {

    // get apropriate branch lengths (there is 1 branch lengths stored per HMM)
    currTBFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + hmmIdx);

    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
  }
  return status;
}

/*
 * convert probabilitiy space values in src to BFGS space values in dest.
 * See OneCellFixLib0TrParam2DegPolyHMM versions and comments for constraints/calculations;
 * these are the same, just scaled up
 */
void AllIndFixLib0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {

  // branch lengths
  double t = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(t)); // set T
  }
}

void AllIndFixLib0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  // branch lengths
  double T = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T)); // set t1
  }
}

/*
 * same as setParamsToEst, but using simParamsToEst
 */
void AllIndFixLib0TrParam2DegPolyHMM::setSimParamsToEst(gsl_vector* params) {
  this->simParamsToEst = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simParamsToEst, params);

  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);
  double currTBFGS = 0;

  // for each HMM, call that HMM's setSimParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // get apropriate branch lengths (there is 1 branch lengths stored per HMM)
    currTBFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + hmmIdx);
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setSimParamsToEst on subclass
    (*this->hmmVec)[hmmIdx]->setSimParamsToEst(currHMMParams);
  }
  gsl_vector_free(currHMMParams);
}

/*
 * dimensions all match, just calls setSimFixedParams for each HMM (which allocs its own vector)
 */
void AllIndFixLib0TrParam2DegPolyHMM::setSimFixedParams(gsl_vector* params) {
  this->simFixedParams = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simFixedParams, params);
  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = (*this->hmmVec)[0]->FIXED_TRANSITION_PROB_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumFixedParams();

  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);
  double alpha = this->getAlpha();
  double beta  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double currLib = 0;

  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // lib
    currLib = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx);
    gsl_vector_set(currHMMParams, hmmLibIdx, currLib);

    gsl_vector_set(currHMMParams, hmmTrIdx + 0, alpha);
    gsl_vector_set(currHMMParams, hmmTrIdx + 1, beta);
    gsl_vector_set(currHMMParams, hmmTrIdx + 2, lambda);
    (*this->hmmVec)[hmmIdx]->setSimFixedParams(currHMMParams);
  }
}

void AllIndFixLib0TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  gsl_vector_set_zero(initGuess);
  if(iter == 0) {
    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.225); // set t
    }
  }
  else if(iter == 1) {
    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.13); // set t
    }
  }
  else if(iter == 2) {
    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.02); // set t
    }
  }
}

