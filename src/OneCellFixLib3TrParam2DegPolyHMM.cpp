#include "OneCellFixLib3TrParam2DegPolyHMM.hpp"


OneCellFixLib3TrParam2DegPolyHMM::OneCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy) : OneCellFixLib3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 2, 1, 1, 1) { // 2 numTrParamsToEst (beta, lambda), 1 numFixedTrParams (alpha), 1 numFixedLibs, 1 numBranches
}
OneCellFixLib3TrParam2DegPolyHMM::OneCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : OneCell3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches) {
  this->maxNumBFGSStarts = 0;

  // set up totalLogEmissionLookup
  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  this->totalLogEmissionLookup = new std::vector<gsl_matrix*>();
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
    this->totalLogEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size())); // transposed
    gsl_matrix_set_zero((*this->totalLogEmissionLookup)[chrIdx]);
  }
}

OneCellFixLib3TrParam2DegPolyHMM::~OneCellFixLib3TrParam2DegPolyHMM() {
  for(std::vector<gsl_matrix*>::iterator it = this->totalLogEmissionLookup->begin(); it != this->totalLogEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->totalLogEmissionLookup;
}

void OneCellFixLib3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
  this->updateTotalLogEmissionLookup();
}
double OneCellFixLib3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 * currChrDepthsVec is ignored. Assumes updateTotalLogEmissionLookup is called
 * every time lib sizes are updated
 */
double OneCellFixLib3TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  gsl_matrix* currChrEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
  return gsl_matrix_get(currChrEmissionMat, depthIdx, stateIdx); // transposed
}

/*
 * updates totalLogEmissionLookup, where each entry corresponds to one chr's matrix.
 * Each matrix is rows: depthIdx, cols: states (trying out since  matrices are stored row major, and we iter through all states in a window before moving to the next window). Each entry is calculated by
 * OneCell3Tr2DegPolyHMM's getTotalLogEmissionProb
 */
void OneCellFixLib3TrParam2DegPolyHMM::updateTotalLogEmissionLookup() {
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  int cellIdx = 0;
  std::string currChr;
  gsl_matrix* currChrEmissionMat = nullptr;

  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChrEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
    currChr = (*chrVec)[chrIdx];
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    for(unsigned int stateIdx = 0; stateIdx < currChrEmissionMat->size2; stateIdx++) { // transposed
      for(unsigned int depthIdx = 0; depthIdx < currChrEmissionMat->size1; depthIdx++) {
        gsl_matrix_set(currChrEmissionMat, depthIdx, stateIdx, OneCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx)); // transposed
      }
    }
  }
}

void OneCellFixLib3TrParam2DegPolyHMM::setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  HMM::setMeanVarianceFn(meanVarianceCoefVec);
  this->updateTotalLogEmissionLookup();
}

/*
 * same as OneCell3TrParam2DegPolyHMM, but lib sizes are skipped
 */
void OneCellFixLib3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  double beta = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t)); // set T
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCellFixLib3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y)); // set beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z)); // set lambda
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T)); // set t
}

/*
 * calls bfgs, returns calling object
 */
OneCellFixLib3TrParam2DegPolyHMM* OneCellFixLib3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  OneCellFixLib3TrParam2DegPolyHMM* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, verbose);
  return bestGuessHMM;
}

