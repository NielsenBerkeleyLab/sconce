#include "AllInd3TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllInd3TrParam2DegPolyHMM::AllInd3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy) : AllInd3TrParam2DegPolyHMM(depths, sampleList, nullptr, maxPloidy, 2, 1, 0, depths->size(), depths->size() * 1) { // 2 shared transition params to est (beta/lambda), 1 fixed transition param (alpha), 0 fixed libs, n numHMMs, n * 1 branches
}
AllInd3TrParam2DegPolyHMM* AllInd3TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy) {
  AllInd3TrParam2DegPolyHMM* hmm = new AllInd3TrParam2DegPolyHMM(depths, sampleList, maxPloidy);
  hmm->makeHMMs();
  return hmm;
}

AllInd3TrParam2DegPolyHMM::AllInd3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst) : AllCells3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numSharedTrParamsToEst, numFixedTrParams, numFixedLibs, numHMMs, numBranchesToEst) {

  // set hmmNames (no need to copy, they're all the same)
  this->hmmNames = this->sampleList;
}

/*
 * helper method to prep for HMM set up
 */
void AllInd3TrParam2DegPolyHMM::makeHMMs() {
  // prep for HMM set up
  gsl_vector* meanVarianceCoefVec = HMM::createMeanVarianceCoefVec();

  // prep transition mat
  gsl_vector* transitionParams = gsl_vector_alloc(3);
  //gsl_vector_set(transitionParams, 0, 0.2); // beta
  //gsl_vector_set(transitionParams, 1, 10); // lambda
  //gsl_vector_set(transitionParams, 2, 0.05);  // t
  gsl_vector_set(transitionParams, 0, 0.01); // beta
  gsl_vector_set(transitionParams, 1, 500); // lambda
  gsl_vector_set(transitionParams, 2, 0.01);  // t

  // save into paramsToEst
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(transitionParams, 0)); // beta
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(transitionParams, 1)); // lambda
  for(int i = 0; i < this->NUM_CELLS; i++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, 2)); // t
  }

  // process all pairs (call subclass's makeHMMs)
  makeHMMs(meanVarianceCoefVec, transitionParams);
  gsl_vector_free(meanVarianceCoefVec);
  gsl_vector_free(transitionParams);

  this->miscFunctions();
}

void AllInd3TrParam2DegPolyHMM::makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams) {
  std::cout << "in AllInd3TrParam2DegPolyHMM::makeHMMs" << std::endl;
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  for(unsigned int hmmIdx = 0; hmmIdx < this->depthsVec->size(); hmmIdx++) {
    currDepths = new std::vector<DepthPair*>();
    currDepths->push_back((*this->depthsVec)[hmmIdx]);
    OneCell3TrParam2DegPolyHMM* hmm = new OneCell3TrParam2DegPolyHMM(currDepths, this->getKploidy());
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

AllInd3TrParam2DegPolyHMM::~AllInd3TrParam2DegPolyHMM() {
  // TODO
}

void AllInd3TrParam2DegPolyHMM::setLibScalingFactor(int cellIdx, double lib) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, lib);
  (*this->hmmVec)[cellIdx]->setLibScalingFactor(0, lib);
}

/*
 * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
 */
void AllInd3TrParam2DegPolyHMM::setAllLibScalingFactors(double libScalingFactor) {
  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
    gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, libScalingFactor);
    (*this->hmmVec)[i]->setLibScalingFactor(0, libScalingFactor);
  }
}
double AllInd3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 * sets transition matrices for all HMMs to be these shared transitionParams.
 * Also saves into this->paramsToEst. Useful for initializing all HMMs at the same transition params
 * (setParamsToEst requires individual branch lengths, libs, etc).
 * //transitionParams should be [beta, gamma, t] (ie shared for all HMMs)
 * transitionParams should be [beta, lambda, t] (ie shared for all HMMs)
 */
double AllInd3TrParam2DegPolyHMM::setAllTransition(gsl_vector* transitionParams) {
  // save beta and lambda
  double beta = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0); // idx should be 0
  //double gamma = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double lambda = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, beta);
  //gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gamma);
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, lambda);

  double t = gsl_vector_get(transitionParams, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 0); // idx should be 2
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx, t);
    status += (*this->hmmVec)[hmmIdx]->setTransition(transitionParams);
  }
  return status;
}


// Optimizable methods
/*
 * method to extract HMM specific variables and set those variables for each HMM,
 * using each HMM's own setParamsToEst method.
 * params is in probability space
 */
double AllInd3TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);
  //std::cout << "AllInd3TrParam2DegPolyHMM setParamsToEst params: " << std::endl;
  //printColVector(params);

  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = (*this->hmmVec)[0]->TRANSITION_PROB_START_IDX;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);
  //std::cout << hmmLibIdx << ", " << hmmTrIdx << ", " << hmmBranchIdx << ", " << currHMMParams->size << std::endl;
  //std::cout << this->LIB_SIZE_SCALING_FACTOR_START_IDX << ", " << this->SHARED_TRANSITION_PROB_START_IDX << ", " << this->BRANCH_LENGTH_START_IDX << std::endl;

  double currLibBFGS = 0;
  double betaBFGS  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double gammaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double lambdaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double currTBFGS = 0;

  double status = 0;
  // for each HMM, call that HMM's setParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    //// alloc new currHMMParams for each HMM (old paramsToEst for each HMM is freed in HMM::setParamsToEst) Thu 30 Jan 2020 10:50:38 AM PST changed to be memcpy
    //currHMMParams = gsl_vector_alloc(numHMMParams);


    // get appropriate libs
    currLibBFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx);

    // get apropriate branch lengths (there is 1 branch lengths stored per HMM)
    currTBFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + hmmIdx);

    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << lambdaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmLibIdx, currLibBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 0, betaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 1, gammaBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 1, lambdaBFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    //std::cout << "IN subclass::SETPARAMSTOEST, currHMMParams:" << std::endl;
    //printColVector(currHMMParams);
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
    //(*this->hmmVec)[hmmIdx]->print(stdout);
  }
  //gsl_vector_free(currHMMParams); // each HMM needs its own copy; don't free this here
  //std::cout << "end of subclass:setParamsToEst" << std::endl;
  //this->print(stdout);
  //std::cout << "#########" << std::endl;
  /*if(status != GSL_SUCCESS) {
    std::cerr << "AllInd3TrParam2DegPolyHMM::setParamsToEst caught nan status, returning nan" << std::endl;
  }*/
  return status;
}

/*
 * convert probabilitiy space values in src to BFGS space values in dest.
 * See OneCell3TrParam2DegPolyHMM versions and comments for constraints/calculations;
 * these are the same, just scaled up
 */
void AllInd3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllInd3TrParam2DegPolyHMM::convertProbToParam" << std::endl;
  // lib scaling factors
  double r = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, log(r / 10.0));// lib10
  }

  //// shared transition parameters
  //double a = this->getAlpha();
  //double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = (b * this->getKploidy() + g - 1);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z

  //// branch lengths
  //double d = (double) (*this->depthsVec)[0]->maxWindowSize;
  //double t = 0;
  //for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
  //  t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0) / d;
  //  c = (d * t - 1);
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(-(d * t) / c)); // set T
  //}

  // shared transition parameters
  double beta = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(lambda / 1.0e6)); // set z

  // branch lengths
  double t = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(t)); // set T
  }
}

void AllInd3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllInd3TrParam2DegPolyHMM::convertParamToProb" << std::endl;
  // lib scaling factors
  double w = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, exp(w) * 10.0); // lib10
  }

  //// shared transition parameters
  //double a = this->getAlpha();
  //double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = 1 - 2.0*a + exp(y) + exp(z);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma

  //// pairwise branch lengths
  //double T = 0;
  //for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
  //  T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
  //  c = 1.0 / (1 + exp(T));
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T) * c); // set t1
  //}

  // shared transition parameters
  double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z)); // lambda
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z) * 1e6); // lambda

  // pairwise branch lengths
  double T = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T)); // set t1
  }
}

/*
 * same as setParamsToEst, but using simParamsToEst
 */
void AllInd3TrParam2DegPolyHMM::setSimParamsToEst(gsl_vector* params) {
  this->simParamsToEst = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simParamsToEst, params);

  int hmmLibIdx = (*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = (*this->hmmVec)[0]->TRANSITION_PROB_START_IDX;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  double currLibBFGS = 0;
  double betaBFGS  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double gammaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double lambdaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
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

    gsl_vector_set(currHMMParams, hmmTrIdx + 0, betaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 1, gammaBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 1, lambdaBFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    (*this->hmmVec)[hmmIdx]->setSimParamsToEst(currHMMParams);
  }
  gsl_vector_free(currHMMParams);
}

/*
 * assumes only libs and alpha are ever fixed
 */
void AllInd3TrParam2DegPolyHMM::setSimFixedParams(gsl_vector* params) {
  this->simFixedParams = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simFixedParams, params);

  int hmmLibIdx = 0;
  int hmmTrIdx = (*this->hmmVec)[0]->FIXED_TRANSITION_PROB_START_IDX;
  int numHMMParams = (*this->hmmVec)[0]->getNumFixedParams();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  double currLibBFGS = 0;
  double alpha = 0;

  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // get appropriate libs
    if(this->NUM_FIXED_LIBS > 0) {
      currLibBFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx);
      gsl_vector_set(currHMMParams, hmmLibIdx, currLibBFGS);
    }
    alpha = gsl_vector_get(params, this->FIXED_TRANSITION_PROB_START_IDX);
    gsl_vector_set(currHMMParams, hmmTrIdx, alpha);
    (*this->hmmVec)[hmmIdx]->setSimFixedParams(currHMMParams);
  }
}

/*
 * this function currently does nothing Tue 28 Jul 2020 06:57:17 PM PDT
 */
void AllInd3TrParam2DegPolyHMM::miscFunctions() {
}

void AllInd3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  //std::cout << "AllInd3TrParam2DegPolyHMM::setInitGuessNthTime iter: " << iter << ", numTotalRuns: " << numTotalRuns << ", this->baumWelchParamResults->size(): " << this->baumWelchParamResults->size() << std::endl;
  gsl_vector_set_zero(initGuess);
  if(this->baumWelchParamResults != nullptr && iter < (int) this->baumWelchParamResults->size()) {
  //if(this->baumWelchParamResults != nullptr) {
    gsl_vector_memcpy(initGuess, (*this->baumWelchParamResults)[iter]); // save iter'th set of params into initGuess
    bool bwParamAdjusted = false; // bwAdj
    for(unsigned int paramIdx = 0; paramIdx < initGuess->size; paramIdx++) {
      double bwParamValue = gsl_vector_get(initGuess, paramIdx);
      if(bwParamValue < 1e-4) {
        bwParamAdjusted = true;
        bwParamValue = 1e-3;
        gsl_vector_set(initGuess, paramIdx, bwParamValue);
      }
    }
    if(bwParamAdjusted) {
      std::cout << "adjusted initGuess values from baumWelchParamResults < 1e-4 to be 1e-3" << std::endl;
      printColVector(initGuess);
    }

    // if passed -1, this is an emergency reset (default params) because the prior transition params failed for this run. reset lib size to libRatio and try different transition params
    //if(numTotalRuns == -1) {
    if(numTotalRuns < 0) {
      /*// libs
      if(this->NUM_LIBS_TO_EST > 0) {
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(0);
          gsl_vector_set(initGuess, initGuessIdx, lib);
          initGuessIdx += 1;
        }
      }*/
      gsl_vector_memcpy(initGuess, (*this->baumWelchParamResults)[iter]); // save iter'th set of params into initGuess to get libs

      // bw libs
      if(numTotalRuns == -1) {
        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      else if(numTotalRuns == -2) {
        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -3) {
        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      // scaled libRatio
      else if(-7 <= numTotalRuns && numTotalRuns <= -4) {
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        double currScaler = 0;
        switch(numTotalRuns % 4) {
          case 0:
            currScaler = 0.5;
            break;
          case -1:
            currScaler = 1;
            break;
          case -2:
            currScaler = 1.5;
            break;
          case -3:
            currScaler = 2;
            break;
        }
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * currScaler;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      // scaled libRatio
      else if(-11 <= numTotalRuns && numTotalRuns <= -8) {
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        double currScaler = 0;
        switch(numTotalRuns % 4) {
          case 0:
            currScaler = 0.5;
            break;
          case -1:
            currScaler = 1;
            break;
          case -2:
            currScaler = 1.5;
            break;
          case -3:
            currScaler = 2;
            break;
        }
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * currScaler;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      // scaled libRatio
      else if(-15 <= numTotalRuns && numTotalRuns <= -12) {
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        double currScaler = 0;
        switch(numTotalRuns % 4) {
          case 0:
            currScaler = 0.5;
            break;
          case -1:
            currScaler = 1;
            break;
          case -2:
            currScaler = 1.5;
            break;
          case -3:
            currScaler = 2;
            break;
        }
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * currScaler;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }

      /*// check if default params match this->initGuessCopyBeforeBFGS (ie baum welch and ls failed, so those set initGuess to orig params. don't to reset to same thing). only checks beta and lambda, assumes t are the same
      if(numTotalRuns == -2 || (compareDoubles(gsl_vector_get(this->initGuessCopyBeforeBFGS, this->SHARED_TRANSITION_PROB_START_IDX + 0), 0.01) && // beta
        compareDoubles(gsl_vector_get(this->initGuessCopyBeforeBFGS, this->SHARED_TRANSITION_PROB_START_IDX + 0), 500))) { // lambda
        std::cout << "Detected beta and lambda values from default params and this->initGuessCopyBeforeBFGS match, or secondary default params requested. Using secondary default params" << std::endl;
        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -3) { // assume that by now, there's something wrong with the lib, boost it upwards
        // lib scaling factors
        for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
          gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.5);
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -4) { // assume that by now, there's something wrong with the lib, boost it upwards
        // lib scaling factors
        for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
          gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 2);
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      else if(numTotalRuns == -5) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      else if(numTotalRuns == -6) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -7) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      else if(numTotalRuns == -8) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      else if(numTotalRuns == -9) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -10) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 1.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      else if(numTotalRuns == -11) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 2;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      else if(numTotalRuns == -12) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 2;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -13) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 2;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      else if(numTotalRuns == -14) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 0.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }
      else if(numTotalRuns == -15) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 0.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
        }
      }
      else if(numTotalRuns == -16) { // assume that by now, there's something wrong with the lib, try multiples of libRatio
        // lib scaling factors
        int initGuessIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
        double lib = 0;
        for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
          for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
            lib = (*this->hmmVec)[hmmIdx]->calcLibScalingFactorsToTotalRatio(cellIdx) * 0.5;
            gsl_vector_set(initGuess, initGuessIdx, lib);
            initGuessIdx += 1;
          }
        }

        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
        }
      }
      else {
        // shared transition parameters
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
        gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

        // branch lengths
        for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
        }
      }*/
    }

    return;
  }
  if(iter == 0) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
    }

    // shared transition parameters
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 500); // lambda

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.01); // set t
    }
  }
  else if(iter == 1) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.5);
    }

    // shared transition parameters
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.1); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 100); // lambda

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.001); // set t
    }
  }
  else if(iter == 2) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 2);
    }

    // shared transition parameters
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 1000); // lambda

    // branch lengths
    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + cellIdx, 0.05); // set t
    }
  }
}

/*
 * internal method to actually calculate the sum of squared residuals for least squares. taken out of setBaumWelchInitGuess.
 * needs to be overridden by subclasses if there are different num and ordering of params
 */
double AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) {
  // convert v from BFGS space to prob space (like regular convertParamToProb)
  gsl_vector* probs = gsl_vector_alloc(v->size);
  this->baumWelchLeastSquares_convertParamToProb(probs, v);

  // first check if any probs are too small (bfgs optim step gets stuck when values are too close to 0)
  double probMin = gsl_vector_min(probs);
  //if(probMin < 1e-3 || gsl_isnan(probMin)) {
  //if(probMin < 1e-4 || gsl_isnan(probMin)) {
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
  //if(probMin < 1e-6 || gsl_isnan(probMin)) {
  //if(probMin < 1e-8 || gsl_isnan(probMin)) { // bwe8
    return GSL_NAN;
  }

  // shortcut for any rates becoming too large
  double probMax = gsl_vector_max(probs);
  if(probMax > 10000.0 || gsl_isnan(probMax)) { // pme4
    return GSL_NAN;
  }

  // pull out relevant params (like setParamsToEst)
  int hmmTrIdx = (*this->hmmVec)[0]->TRANSITION_PROB_START_IDX - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst() - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
  gsl_vector* currHMMProbs = gsl_vector_alloc(numHMMParams);
  //gsl_vector* oldHMMProbs = gsl_vector_alloc((*this->hmmVec)[0]->getNumParamsToEst());

  double beta  = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  double lambda = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double currT = 0;

  double sumSqResid = 0;
  double status = 0;
  double unifTransitionMatrix = 0;
  // for each HMM
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // get apropriate branch length
    currT = gsl_vector_get(probs, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + hmmIdx);

    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << gammaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
    // set everything into currHMMProbs
    gsl_vector_set(currHMMProbs, hmmTrIdx + 0, beta);
    gsl_vector_set(currHMMProbs, hmmTrIdx + 1, lambda);
    gsl_vector_set(currHMMProbs, hmmBranchIdx + 0, currT);
    //std::cout << beta << ", " << lambda << ", " << currT << std::endl;
    //printColVector(currHMMProbs);

    // call TwoCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs)
    //gsl_vector_memcpy(oldHMMProbs, (*this->hmmVec)[hmmIdx]->getParamsToEst());
    status = (*this->hmmVec)[hmmIdx]->baumWelchLeastSquares_f(currHMMProbs);
    if(gsl_isnan(status)) {
      fprintf(stderr, "ERROR: could not set transition matrix in baumWelchLeastSquares_f\t|\t");
      printRowVector(stderr, currHMMProbs);
      sumSqResid = status;
      //(*this->hmmVec)[hmmIdx]->setParamsToEst(oldHMMProbs); // reset to old params if they failed Wed 26 May 2021 05:28:45 PM PDT don't do this
      break;
    }

    // check if these params led to uniform matrix
    unifTransitionMatrix = (*this->hmmVec)[hmmIdx]->checkStateValidity((*this->hmmVec)[hmmIdx]->getBaumWelchTransitionMat());
    if(gsl_isnan(unifTransitionMatrix)) {
      fprintf(stderr, "ERROR: params resulted in uniform transition matrix in baumWelchLeastSquares_f\t|\t");
      printRowVector(stderr, currHMMProbs);
      sumSqResid = GSL_NAN;
      break;
    }

    sumSqResid += status;
  }
  gsl_vector_free(probs);
  gsl_vector_free(currHMMProbs);
  //gsl_vector_free(oldHMMProbs);
  return sumSqResid;
}

/*
 * this is to be called during setBaumWelchInitGuess, in order to save params
 * estimated after least squares back into this->paramsToEst and initGuess.
 *
 * Each HMM already has params set correctly, just need to save them in this object. varsToEst_probSpace doesn't contain libs (varsToESt_probSpace is from least squares), so those need to be saved first (if any)
 *
 * varsToEst_probSpace = [beta, gamma, t_cell0, t_cell1, ..., t_cellN]
 * initGuess = paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t_cell0, t_cell1, ..., t_cellN]
 */
void AllInd3TrParam2DegPolyHMM::saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) {
  gsl_vector* estsWithLibs = gsl_vector_alloc(this->getNumParamsToEst());
  gsl_vector_set_zero(estsWithLibs);
  int estsWithLibsIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
  double lib = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    lib = (*this->hmmVec)[hmmIdx]->getLibScalingFactor(0);
    this->setLibScalingFactor(hmmIdx, lib);
    if(this->NUM_LIBS_TO_EST > 0) {
      gsl_vector_set(estsWithLibs, estsWithLibsIdx, lib);
      estsWithLibsIdx += 1;
    }
  }
  // copy the rest of varsToEst_probSpace
  for(unsigned int varsToEstIdx = 0; varsToEstIdx < varsToEst_probSpace->size; varsToEstIdx++, estsWithLibsIdx++) {
    gsl_vector_set(estsWithLibs, estsWithLibsIdx, gsl_vector_get(varsToEst_probSpace, varsToEstIdx));
  }
  this->setParamsToEst(estsWithLibs);
  gsl_vector_memcpy(initGuess, estsWithLibs);
  //gsl_vector_free(varsToEst_probSpace);
  gsl_vector_free(estsWithLibs);
}

// same as convertProbToParam, just without library sizes
void AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  ////std::cout << "AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam" << std::endl;
  //// shared transition parameters
  //double a = this->getAlpha();
  //double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = (b * this->getKploidy() + g - 1);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z

  //// branch lengths
  //double d = (double) (*this->depthsVec)[0]->maxWindowSize;
  //double t = 0;
  //for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
  //  t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0) / d;
  //  c = (d * t - 1);
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(-(d * t) / c)); // set T
  //}

  //std::cout << "AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam" << std::endl;
  // shared transition parameters
  double beta = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  //double lambda = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double lambda = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, log(beta)); // set y
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, log(lambda / 1.0e6)); // set z

  // branch lengths
  double t = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + cellIdx);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + cellIdx, log(t)); // set T
  }
}

void AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  ////std::cout << "AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb" << std::endl;
  //// shared transition parameters
  //double a = this->getAlpha();
  //double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = 1 - 2.0*a + exp(y) + exp(z);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma

  //// pairwise branch lengths
  //double T = 0;
  //for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
  //  T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
  //  c = 1.0 / (1 + exp(T));
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T) * c); // set t1
  //}

  //std::cout << "AllInd3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb" << std::endl;
  // shared transition parameters
  double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, exp(y)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, exp(z)); // lambda
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, exp(z) * 1e6); // lambda

  // pairwise branch lengths
  double T = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + cellIdx);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + cellIdx, exp(T)); // set t1
  }
}

