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
  //gsl_vector_set(transitionParams, 0, 0.002); // beta
  //gsl_vector_set(transitionParams, 1, 0.04); // gamma
  gsl_vector_set(transitionParams, 0, 0.05);  // t

  // save into paramsToEst
  //gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(transitionParams, 0)); // beta
  //gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(transitionParams, 1)); // gamma
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

//void AllInd0TrParam2DegPolyHMM::setLibScalingFactor(int cellIdx, double lib) {
//  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, lib);
//  (*this->hmmVec)[cellIdx]->setLibScalingFactor(0, lib);
//}
//
///*
// * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
// */
//void AllInd0TrParam2DegPolyHMM::setAllLibScalingFactors(double libScalingFactor) {
//  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
//    gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, libScalingFactor);
//    (*this->hmmVec)[i]->setLibScalingFactor(0, libScalingFactor);
//  }
//}
//double AllInd0TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
//  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
//}

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

    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << gammaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmLibIdx, currLibBFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx, currTBFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    //std::cout << "IN subclass::SETPARAMSTOEST, currHMMParams:" << std::endl;
    //printColVector(currHMMParams);
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
    //(*this->hmmVec)[hmmIdx]->print(stdout);
  }
  return status;
}

/*
 * convert probabilitiy space values in src to BFGS space values in dest.
 * See OneCell0TrParam2DegPolyHMM versions and comments for constraints/calculations;
 * these are the same, just scaled up
 */
void AllInd0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllInd0TrParam2DegPolyHMM::convertProbToParam" << std::endl;
  // lib scaling factors
  double r = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, log(r));
  }

  // branch lengths
  //double d = (double) (*this->depthsVec)[0]->maxWindowSize;
  double t = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    //t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0) / d;
    //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(-(d * t) / (d * t - 1))); // set T
    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(t)); // set T
  }
}

void AllInd0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllInd0TrParam2DegPolyHMM::convertParamToProb" << std::endl;
  // lib scaling factors
  double w = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, exp(w));
  }

  // branch lengths
  double T = 0;
  //double c = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    //T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
    //c = 1.0 / (1 + exp(T));
    //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T) * c); // set t1
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

///*
// * this function currently does nothing Tue 28 Jul 2020 06:57:17 PM PDT
// */
//void AllInd0TrParam2DegPolyHMM::miscFunctions() {
//}

/*
 * because no params are shared between any HMMs (libs and branches only appear once), each call to bfgs is independent and can be easily parallelized.
 * This is based on AllPairsFixLib0TrParam2DegPolyHMM::bfgs
 */
AllInd0TrParam2DegPolyHMM* AllInd0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  //std::cout << "AllInd0TrParam2DegPolyHMM::bfgs" << std::endl;
  /*AllInd0TrParam2DegPolyHMM* bestGuessOptim = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, verbose);
  return bestGuessOptim;*/

  // create new HMM with the best guess parameters and return it
  AllInd0TrParam2DegPolyHMM* bestGuessOptim = this;
  gsl_vector* savedBestEstParams = bestGuessOptim->getParamsToEst();
  /*gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, verbose);*/

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

///*
// * internal method to actually calculate the sum of squared residuals for least squares. taken out of setBaumWelchInitGuess.
// * needs to be overridden by subclasses if there are different num and ordering of params
// */
//double AllInd0TrParam2DegPolyHMM::baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) {
//  // convert v from BFGS space to prob space (like regular convertParamToProb)
//  gsl_vector* probs = gsl_vector_alloc(v->size);
//  this->baumWelchLeastSquares_convertParamToProb(probs, v);
//
//  // pull out relevant params (like setParamsToEst)
//  int hmmTrIdx = (*this->hmmVec)[0]->TRANSITION_PROB_START_IDX - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
//  int hmmBranchIdx = (*this->hmmVec)[0]->BRANCH_LENGTH_START_IDX - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
//  int numHMMParams = (*this->hmmVec)[0]->getNumParamsToEst() - (*this->hmmVec)[0]->NUM_LIBS_TO_EST;
//  gsl_vector* currHMMProbs = gsl_vector_alloc(numHMMParams);
//
//  double betaBFGS  = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
//  double gammaBFGS = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
//  double currTBFGS = 0;
//
//  double sumSqResid = 0;
//  // for each HMM
//  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
//    // get apropriate branch length
//    currTBFGS = gsl_vector_get(probs, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + hmmIdx);
//
//    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << gammaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
//    // set everything into currHMMProbs
//    gsl_vector_set(currHMMProbs, hmmTrIdx + 0, betaBFGS);
//    gsl_vector_set(currHMMProbs, hmmTrIdx + 1, gammaBFGS);
//    gsl_vector_set(currHMMProbs, hmmBranchIdx + 0, currTBFGS);
//
//    // call TwoCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs)
//    sumSqResid += (*this->hmmVec)[hmmIdx]->baumWelchLeastSquares_f(currHMMProbs);
//  }
//  gsl_vector_free(probs);
//  gsl_vector_free(currHMMProbs);
//  return sumSqResid;
//}
//
///*
// * this is to be called during setBaumWelchInitGuess, in order to save params
// * estimated after least squares back into this->paramsToEst and initGuess.
// *
// * Each HMM already has params set correctly, just need to save them in this object. varsToEst_probSpace doesn't contain libs (varsToESt_probSpace is from least squares), so those need to be saved first (if any)
// *
// * varsToEst_probSpace = [beta, gamma, t_cell0, t_cell1, ..., t_cellN]
// * initGuess = paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t_cell0, t_cell1, ..., t_cellN]
// */
//void AllInd0TrParam2DegPolyHMM::saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) {
//  gsl_vector* estsWithLibs = gsl_vector_alloc(this->getNumParamsToEst());
//  gsl_vector_set_zero(estsWithLibs);
//  int estsWithLibsIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
//  double lib = 0;
//  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
//    lib = (*this->hmmVec)[hmmIdx]->getLibScalingFactor(0);
//    this->setLibScalingFactor(hmmIdx, lib);
//    if(this->NUM_LIBS_TO_EST > 0) {
//      gsl_vector_set(estsWithLibs, estsWithLibsIdx, lib);
//      estsWithLibsIdx += 1;
//    }
//  }
//  // copy the rest of varsToEst_probSpace
//  for(unsigned int varsToEstIdx = 0; varsToEstIdx < varsToEst_probSpace->size; varsToEstIdx++, estsWithLibsIdx++) {
//    gsl_vector_set(estsWithLibs, estsWithLibsIdx, gsl_vector_get(varsToEst_probSpace, varsToEstIdx));
//  }
//  this->setParamsToEst(estsWithLibs);
//  gsl_vector_memcpy(initGuess, estsWithLibs);
//  gsl_vector_free(varsToEst_probSpace);
//  gsl_vector_free(estsWithLibs);
//}
//
//// same as convertProbToParam, just without library sizes
//void AllInd0TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
//  //std::cout << "AllInd0TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam" << std::endl;
//  // shared transition parameters
//  double a = this->getAlpha();
//  double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
//  double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
//  double c = (b * this->getKploidy() + g - 1);
//  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
//  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z
//
//  // branch lengths
//  double d = (double) (*this->depthsVec)[0]->maxWindowSize;
//  double t = 0;
//  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
//    t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx + 0) / d;
//    c = (d * t - 1);
//    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx, log(-(d * t) / c)); // set T
//  }
//}
//
//void AllInd0TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
//  //std::cout << "AllInd0TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb" << std::endl;
//  // shared transition parameters
//  double a = this->getAlpha();
//  double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
//  double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
//  double c = 1 - 2.0*a + exp(y) + exp(z);
//  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
//  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma
//
//  // pairwise branch lengths
//  double T = 0;
//  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
//    T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + cellIdx);
//    c = 1.0 / (1 + exp(T));
//    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + cellIdx + 0, exp(T) * c); // set t1
//  }
//}


