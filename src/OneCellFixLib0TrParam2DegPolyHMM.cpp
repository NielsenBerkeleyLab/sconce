#include "OneCellFixLib0TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
OneCellFixLib0TrParam2DegPolyHMM::OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy) : OneCellFixLib0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 0, 3, 1, 1) { // fixedParams, 0 transition params to est, 3 fixedTrParams, 1 fixedLibs, 1 branch
}
OneCellFixLib0TrParam2DegPolyHMM::OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : OneCell0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches) {
  this->maxNumBFGSStarts = 0;

  // set up totalLogEmissionLookup
  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  this->totalLogEmissionLookup = new std::vector<gsl_matrix*>();
  this->totalEmissionLookup = new std::vector<gsl_matrix*>();
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
    this->totalLogEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size()));
    this->totalEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size()));
    gsl_matrix_set_zero((*this->totalLogEmissionLookup)[chrIdx]);
    gsl_matrix_set_zero((*this->totalEmissionLookup)[chrIdx]);
  }
}

OneCellFixLib0TrParam2DegPolyHMM::~OneCellFixLib0TrParam2DegPolyHMM() {
  for(std::vector<gsl_matrix*>::iterator it = this->totalLogEmissionLookup->begin(); it != this->totalLogEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  for(std::vector<gsl_matrix*>::iterator it = this->totalEmissionLookup->begin(); it != this->totalEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->totalLogEmissionLookup;
  delete this->totalEmissionLookup;
}

/*
 ********
 * accessors and mutators
 ********
 */
void OneCellFixLib0TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
  this->updateAllTotalEmissionLookup();
}
double OneCellFixLib0TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 * same as OneCellFixLib3TrParam2DegPolyHMM
 */
double OneCellFixLib0TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  gsl_matrix* currChrEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
  return gsl_matrix_get(currChrEmissionMat, depthIdx, stateIdx);
}
double OneCellFixLib0TrParam2DegPolyHMM::getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  gsl_matrix* currChrEmissionMat = (*this->totalEmissionLookup)[chrIdx];
  return gsl_matrix_get(currChrEmissionMat, depthIdx, stateIdx);
}

/*
 * same as OneCellFixLib3TrParam2DegPolyHMM
 */
void OneCellFixLib0TrParam2DegPolyHMM::updateAllTotalEmissionLookup() {
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  int cellIdx = 0;
  std::string currChr;
  gsl_matrix* currChrLogEmissionMat = nullptr;
  gsl_matrix* currChrEmissionMat = nullptr;
  double currLogEmi = 0;

  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChrLogEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
    currChrEmissionMat = (*this->totalEmissionLookup)[chrIdx];
    currChr = (*chrVec)[chrIdx];
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    for(unsigned int stateIdx = 0; stateIdx < currChrLogEmissionMat->size2; stateIdx++) {
      for(unsigned int depthIdx = 0; depthIdx < currChrLogEmissionMat->size1; depthIdx++) {
        currLogEmi = OneCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);
        gsl_matrix_set(currChrLogEmissionMat, depthIdx, stateIdx, currLogEmi);
        gsl_matrix_set(currChrEmissionMat, depthIdx, stateIdx, exp(currLogEmi));
      }
    }
  }
}
void OneCellFixLib0TrParam2DegPolyHMM::setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  HMM::setMeanVarianceFn(meanVarianceCoefVec);
  this->updateAllTotalEmissionLookup();
}

/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
void OneCellFixLib0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t));
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCellFixLib0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T)); // set t
}

/*
 * calls bfgs, returns a bestGuessHMM
 */
OneCellFixLib0TrParam2DegPolyHMM* OneCellFixLib0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  // create new HMM with the best guess parameters and return it
  OneCellFixLib0TrParam2DegPolyHMM* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  bestGuessHMM->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, verbose);
  return bestGuessHMM;
}

void OneCellFixLib0TrParam2DegPolyHMM::miscFunctions() {
  // do nothing. This is usually used to kick off viterbi decoding for lib size est
  // but in this class, lib sizes are fixed
}

