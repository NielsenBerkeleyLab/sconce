#include "AllInd0TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllInd0TrParam2DegPolyHMM::AllInd0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy) : AllInd0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 0, 3, 0, depths->size(), depths->size() * 1) { // 0 shared transition params to est, 3 fixed transition param (alpha/beta/lambda), 0 fixed libs, n HMMs, n * 1 branches
}
AllInd0TrParam2DegPolyHMM* AllInd0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy) {
  AllInd0TrParam2DegPolyHMM* hmm = new AllInd0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy);
  hmm->makeHMMs();
  return hmm;
}

AllInd0TrParam2DegPolyHMM::AllInd0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst) : AllInd3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numSharedTrParamsToEst, numFixedTrParams, numFixedLibs, numHMMs, numBranchesToEst) {
}

/*
 * helper method to prep for HMM set up
 */
void AllInd0TrParam2DegPolyHMM::makeHMMs() {
  // prep for HMM set up
  gsl_vector* meanVarianceCoefVec = HMM::createMeanVarianceCoefVec();

  // prep transition mat
  gsl_vector* transitionParams = gsl_vector_alloc(1);
  gsl_vector_set(transitionParams, 0, 0.05);  // t

  // save into paramsToEst
  for(int i = 0; i < this->NUM_CELLS; i++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, 0)); // t
  }

  // process all pairs (call subclass's makeHMMs)
  makeHMMs(meanVarianceCoefVec, transitionParams);
  gsl_vector_free(meanVarianceCoefVec);
  gsl_vector_free(transitionParams);

  this->miscFunctions();
}

void AllInd0TrParam2DegPolyHMM::makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams) {
  std::cout << "in AllInd0TrParam2DegPolyHMM::makeHMMs" << std::endl;
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  gsl_vector* currFixedParams = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->depthsVec->size(); hmmIdx++) {
    currDepths = new std::vector<DepthPair*>();
    currDepths->push_back((*this->depthsVec)[hmmIdx]);
    currFixedParams = gsl_vector_alloc(this->fixedParams->size);
    gsl_vector_memcpy(currFixedParams, this->fixedParams);
    OneCell0TrParam2DegPolyHMM* hmm = new OneCell0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy());
    (*this->hmmVec)[hmmIdx] = hmm;

    // do rest of HMM set up (usually happens in main.cpp)
    currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
    gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
    hmm->setMeanVarianceFn(currMeanVarCoefVec);
    hmm->setTransition(transitionParams); // this vector isn't saved anywhere
    hmm->setLibScalingFactorsToTotalRatio();
    hmm->setAlpha(this->getAlpha());

    this->setLibScalingFactor(hmmIdx, hmm->getLibScalingFactor(0));
  }
}

AllInd0TrParam2DegPolyHMM::~AllInd0TrParam2DegPolyHMM() {
  // TODO
}


/*
 * same as AllInd3TrParam2DegPolyHMM::setAllTransition but skips over beta/lambda terms, and only takes branches
 * for this->paramsToEst and for initializing all HMMs to same branches
 * branches should be [t] (ie shared for all HMMs)
 */
double AllInd0TrParam2DegPolyHMM::setAllBranches(gsl_vector* branches) {
  double t = gsl_vector_get(branches, 0);

  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx, t);
    status += (*this->hmmVec)[hmmIdx]->setTransition(branches);
  }

  return status;
}


// Optimizable methods
/*
 * method to extract HMM specific variables and set those variables for each HMM,
 * using each HMM's own setParamsToEst method.
 * params is in probability space
 */
double AllInd0TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);

  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  double currLibBFGS = 0;
  double currTBFGS = 0;

  double status = 0;
  // for each HMM, call that HMM's setParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // get appropriate lib
    currLibBFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx);

    // get apropriate branch lengths (there is 1 branch lengths stored per HMM)
    currTBFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + hmmIdx);

    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmLibIdx, currLibBFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
  }
  return status;
}

/*
 * convert probabilitiy space values in src to BFGS space values in dest.
 * See OneCell0TrParam2DegPolyHMM versions and comments for constraints/calculations;
 * these are the same, just scaled up
 */
void AllInd0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  // lib scaling factors
  double r = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, log(r));
  }

  // branch lengths
  double t = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(t)); // set T
  }
}

void AllInd0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  // lib scaling factors
  double w = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, exp(w));
  }

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
void AllInd0TrParam2DegPolyHMM::setSimParamsToEst(gsl_vector* params) {
  this->simParamsToEst = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simParamsToEst, params);

  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  double currLibBFGS = 0;
  double currTBFGS = 0;

  // for each HMM, call that HMM's setSimParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // get appropriate libs
    if(this->NUM_LIBS_TO_EST > 0) {
      currLibBFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx);
      gsl_vector_set(currHMMParams, hmmLibIdx, currLibBFGS);
    }

    // get apropriate branch lengths (there is 1 branch lengths stored per HMM)
    currTBFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + hmmIdx);
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    (*this->hmmVec)[hmmIdx]->setSimParamsToEst(currHMMParams);
  }
  gsl_vector_free(currHMMParams);
}

/*
 * dimensions all match, just calls setSimFixedParams for each HMM (which allocs its own vector)
 */
void AllInd0TrParam2DegPolyHMM::setSimFixedParams(gsl_vector* params) {
  this->simFixedParams = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simFixedParams, params);

  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    (*this->hmmVec)[hmmIdx]->setSimFixedParams(this->simFixedParams);
  }
}

/*
 * because no params are shared between any HMMs (libs and branches only appear once), each call to bfgs is independent and can be easily parallelized.
 * This is based on AllPairsFixLib0TrParam2DegPolyHMM::bfgs
 */
AllInd0TrParam2DegPolyHMM* AllInd0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  // create new HMM with the best guess parameters and return it
  AllInd0TrParam2DegPolyHMM* bestGuessOptim = this;
  gsl_vector* savedBestEstParams = bestGuessOptim->getParamsToEst();

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double initTotalLik = this->getLogLikelihood();
  gsl_vector* currInitGuess = gsl_vector_alloc((*this->hmmVec)[0]->getNumParamsToEst());
  gsl_vector* currBestEst = nullptr;
  int currInitGuessIdx = 0;
  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // extract relevant lib and branch length to est
    currInitGuessIdx = 0;
    if(this->NUM_LIBS_TO_EST > 0) {
      gsl_vector_set(currInitGuess, currInitGuessIdx, gsl_vector_get(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx));
      currInitGuessIdx++;
    }
    gsl_vector_set(currInitGuess, currInitGuessIdx, gsl_vector_get(initGuess, this->BRANCH_LENGTH_START_IDX + hmmIdx));

    // call each HMM's bfgs
    if(verbose) {
      std::cout << "\nCalling BFGS on HMM " << hmmIdx << "(" << (*this->sampleList)[hmmIdx] << "):" << std::endl;
      std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before each bfgs
    }
    HMM* currBestHMM = (HMM*)(*this->hmmVec)[hmmIdx]->bfgs(currInitGuess, verbose);

    // save params into bestGuessOptim
    currBestEst = currBestHMM->getParamsToEst();
    if(this->NUM_LIBS_TO_EST > 0) {
      gsl_vector_set(savedBestEstParams, hmmIdx, gsl_vector_get(currBestEst, hmmLibIdx));
    }
    gsl_vector_set(savedBestEstParams, hmmIdx + this->NUM_LIBS_TO_EST, gsl_vector_get(currBestEst, hmmBranchIdx)); // offset indexing by number of libs estimated
    (*bestGuessOptim->hmmVec)[hmmIdx] = currBestHMM;
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
  double totalLik = this->getLogLikelihood();
  printf("STAGE 2 TOTAL LOGLIKELIHOOD FOUND: %.40f\n", totalLik);
  printf("CHANGE IN STAGE 2 TOTAL LIKELIHOOD: %.40f\n", totalLik - initTotalLik);
  printf("STAGE 2 TOTAL TIME (sec): %.10f\n", elapsedSec);
  printf("DONE WITH STAGE 2 TOTAL (AllPairsFixLib0TrParam2DegPolyHMM):\n\n");
  this->print(stdout);
  std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file

  return bestGuessOptim;
}

void AllInd0TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  gsl_vector_set_zero(initGuess);
  if(iter == 0) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
    }

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.225); // set t
    }
  }
  else if(iter == 1) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, .75);
    }

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.13); // set t
    }
  }
  else if(iter == 2) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.25);
    }

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.02); // set t
    }
  }
}

