#include "OneCellFixLib0TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
OneCellFixLib0TrParam2DegPolyHMM::OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy) : OneCellFixLib0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 0, 3, 1, 1) { // fixedParams, 0 transition params to est, 3 fixedTrParams, 1 fixedLibs, 1 branch
  //this->updateTotalLogEmissionLookup(); // call here, not in other ctor, since it relies on polymorphism to call correct setTransition() // Mon 27 Jul 2020 07:36:54 PM PDT causing problems because setMeanVarianceFn hasn't been called yet
}
OneCellFixLib0TrParam2DegPolyHMM::OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : OneCell0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches) {
  this->maxNumBFGSStarts = 0;

  // set up totalLogEmissionLookup
  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  this->totalLogEmissionLookup = new std::vector<gsl_matrix*>();
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
    this->totalLogEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size()));
    gsl_matrix_set_zero((*this->totalLogEmissionLookup)[chrIdx]);
  }
  //this->updateTotalLogEmissionLookup();

  /*// set fixedParams with library sizes
  double totalAvgDiploidDepth = (*this->depths)[0]->getTotalDiploidDepth();
  this->setLibScalingFactor(0, (*this->depths)[0]->getTotalTumorDepth() / totalAvgDiploidDepth);
  this->setLibScalingFactor(1, (*this->depths)[1]->getTotalTumorDepth() / totalAvgDiploidDepth);*/ // Mon 27 Jan 2020 07:10:45 PM PST passed fixedParams already contains lib sizes
}

OneCellFixLib0TrParam2DegPolyHMM::~OneCellFixLib0TrParam2DegPolyHMM() {
  for(std::vector<gsl_matrix*>::iterator it = this->totalLogEmissionLookup->begin(); it != this->totalLogEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->totalLogEmissionLookup;
}

/*
 ********
 * accessors and mutators
 ********
 */
void OneCellFixLib0TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
  this->updateTotalLogEmissionLookup();
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

/*
 * same as OneCellFixLib3TrParam2DegPolyHMM
 */
void OneCellFixLib0TrParam2DegPolyHMM::updateTotalLogEmissionLookup() {
  //std::cout << "OneCellFixLib0TrParam2DegPolyHMM::updateTotalLogEmissionLookup" << std::endl;
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  int cellIdx = 0;
  std::string currChr;
  gsl_matrix* currChrEmissionMat = nullptr;

  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChrEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
    currChr = (*chrVec)[chrIdx];
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    for(unsigned int stateIdx = 0; stateIdx < currChrEmissionMat->size2; stateIdx++) {
      for(unsigned int depthIdx = 0; depthIdx < currChrEmissionMat->size1; depthIdx++) {
        gsl_matrix_set(currChrEmissionMat, depthIdx, stateIdx, OneCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx));
      }
    }
  }
}
void OneCellFixLib0TrParam2DegPolyHMM::setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  HMM::setMeanVarianceFn(meanVarianceCoefVec);
  this->updateTotalLogEmissionLookup();
}

/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
void OneCellFixLib0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //double d = (double) (*this->depths)[0]->maxWindowSize;
  //double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX) / d;
  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(-(d * t) / (d * t - 1))); // set T
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t));
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCellFixLib0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  //double c = 1.0 / (1 + exp(T));
  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T) * c); // set t
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
  //std::cerr << "in OneCellFixLib0TrParam2DegPolyHMM::miscFunctions(), doing nothing" << std::endl;
}

