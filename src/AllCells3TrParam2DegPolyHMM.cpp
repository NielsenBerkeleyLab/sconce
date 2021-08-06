#include "AllCells3TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllCells3TrParam2DegPolyHMM::AllCells3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst) :
                                        MAX_PLOIDY(maxPloidy),
                                        NUM_CELLS(depths->size()),
                                        NUM_HMMS(numHMMs),
                                        NUM_LIBS_TO_EST(this->NUM_CELLS - numFixedLibs),
                                        NUM_SHARED_TRANSITION_PARAMS_TO_EST(numSharedTrParamsToEst),
                                        NUM_BRANCH_LENGTHS_TO_EST(numBranchesToEst),
                                        NUM_FIXED_LIBS(numFixedLibs),
                                        NUM_FIXED_SHARED_TRANSITION_PARAMS(numFixedTrParams),
                                        LIB_SIZE_SCALING_FACTOR_START_IDX(0),
                                        SHARED_TRANSITION_PROB_START_IDX(this->NUM_CELLS - numFixedLibs), // libs go first, estimating one for each "unfixed lib" cell
                                        BRANCH_LENGTH_START_IDX(this->SHARED_TRANSITION_PROB_START_IDX + numSharedTrParamsToEst),
                                        FIXED_TRANSITION_PROB_START_IDX(numFixedLibs) // if no libs are fixed, then this starts at 0
                                        {
  this->depthsVec = depths;
  this->sampleList = sampleList;
  this->hmmVec = new std::vector<HMM*>(this->NUM_HMMS);
  this->hmmNames = nullptr;
  this->baumWelchParamResults = nullptr;
  this->paramsToEst = gsl_vector_alloc(this->NUM_LIBS_TO_EST + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST); // lib scaling for each cell + b/L + t1/t2/t3 for each pair
  gsl_vector_set_all(this->paramsToEst, 1.0);
  if(fixedParams == nullptr) {
    this->fixedParams = gsl_vector_alloc(numFixedLibs + numFixedTrParams);
    gsl_vector_set_zero(this->fixedParams);
  }
  else {
    this->fixedParams = fixedParams;
  }
  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst());
  this->maxNumBFGSStarts = 0;

  this->setAlpha(0.1); // arbitrary starting value
}

AllCells3TrParam2DegPolyHMM::~AllCells3TrParam2DegPolyHMM() {
  // TODO
}

int AllCells3TrParam2DegPolyHMM::getKploidy() const {
  return this->MAX_PLOIDY;
}

std::vector<HMM*>* AllCells3TrParam2DegPolyHMM::getHMMs() {
  return this->hmmVec;
}

void AllCells3TrParam2DegPolyHMM::print(FILE* stream) {
  // print each HMM
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    fprintf(stream, "HMM %i (%s):\n", hmmIdx, (*this->hmmNames)[hmmIdx].c_str());
    (*this->hmmVec)[hmmIdx]->print(stream);
    fprintf(stream, "\n");
  }

  fprintf(stream, "\n\n");

  // print paramsToEst
  fprintf(stream, "AllCells_TrParam2DegPolyHMM paramsToEst:\n");
  printColVector(stream, this->paramsToEst);

  // print fixedParams
  fprintf(stream, "AllCells_TrParam2DegPolyHMM fixedParams:\n");
  printColVector(stream, this->fixedParams);

  // print total loglikelihood
  fprintf(stream, "AllCells_TrParam2DegPolyHMM total loglikelihood: %.10f\n\n", this->runForwardAlg());
}

std::vector<std::string>* AllCells3TrParam2DegPolyHMM::getSampleList() {
  return this->sampleList;
}

void AllCells3TrParam2DegPolyHMM::setSampleList(std::vector<std::string>* sampleList) {
  if(this->sampleList != nullptr && this->sampleList != sampleList) {
    delete this->sampleList;
  }
  this->sampleList = sampleList;
}

std::vector<gsl_vector*>* AllCells3TrParam2DegPolyHMM::getBaumWelchParamResults() const {
  return this->baumWelchParamResults;
}
void AllCells3TrParam2DegPolyHMM::setBaumWelchParamResults(std::vector<gsl_vector*>* paramResults) {
  this->baumWelchParamResults = paramResults;
}

double AllCells3TrParam2DegPolyHMM::getAlpha() const {
  return gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX);
}

void AllCells3TrParam2DegPolyHMM::setAlpha(double alpha) {
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX, alpha);
}

void AllCells3TrParam2DegPolyHMM::setAllMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  gsl_vector* currMeanVarCoefVec = nullptr;
  HMM* hmm = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    hmm = (*this->hmmVec)[hmmIdx];
    currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
    gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
    hmm->setMeanVarianceFn(currMeanVarCoefVec);
  }
}

// Optimizable methods
/*
 * same as HMM::checkOptimProbValidity
 */
double AllCells3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  //return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  //std::cout << "AllCells3TrParam2DegPolyHMM::checkOptimProbValidity" << std::endl;
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  double storedLib = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    //storedLib = gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    storedLib = gsl_vector_get(this->initGuessCopyBeforeBFGS, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    //if(currLibSizeScalingFactor < 5e-3 || currLibSizeScalingFactor > 1e3) {
    //if(currLibSizeScalingFactor < 5e-3 || currLibSizeScalingFactor > 1e10) {
    //std::cout << "change in lib is: " << storedLib << ", " << currLibSizeScalingFactor << ", " <<std::abs(storedLib - currLibSizeScalingFactor) / storedLib << std::endl;
    //if(gsl_isinf(currLibSizeScalingFactor) || gsl_isnan(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
    //if(gsl_isinf(currLibSizeScalingFactor) || gsl_isnan(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2 || std::abs(storedLib - currLibSizeScalingFactor) / storedLib > 0.25) { // libguard can't go more than .25% in either direction
    if(gsl_isinf(currLibSizeScalingFactor) || gsl_isnan(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2 || (storedLib - currLibSizeScalingFactor) / storedLib > 0.25) { // libguard can get arbitrarily large but can't get more than .25% smaller
      //std::cerr << "shortcut for bad lib sizes" << std::endl;
      return GSL_NAN;
    }
  }

  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  //if(compareDoubles(probMin, 0, 1e-5) || probMin < 0 || gsl_isnan(probMin)) {
  //if(gsl_isnan(probMin)) {
  //if(probMin < 1e-3 || gsl_isnan(probMin)) {
  //if(probMin < 1e-4 || gsl_isnan(probMin)) {
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
  //if(probMin < 1e-10 || gsl_isnan(probMin)) { // optVal
    //std::cerr << "shortcut for any probabilities becoming too small: " << probMin << std::endl;
    //printRowVector(stderr, (gsl_vector*)v);
    //printRowVector(stderr, probs);
    return GSL_NAN;
  }

  // shortcut for any rates becoming too large
  double probMax = gsl_vector_max(probs);
  if(probMax > 10000.0 || gsl_isnan(probMax)) { // pme4
    return GSL_NAN;
  }
  return 0;
}

/*
 * function to check if the state of the HMMs is currently ok.
 * That is, returns nan if at least one HMM has a completely
 * uniform transition matrix, 0 o/w
 */
double AllCells3TrParam2DegPolyHMM::checkStateValidity(double epsilon) const {
  HMM* currHMM = nullptr;
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    currHMM = (*this->hmmVec)[hmmIdx];
    status = status + currHMM->checkStateValidity(epsilon);
  }
  return status;
}


/*
 * Function to check for transient states
 * (ex state 1 is only seen in 0->1->2 or 2->1->0 transitions)
 * Returns nan if transient states are found, 0 o/w
 */
double AllCells3TrParam2DegPolyHMM::checkForTransientStates() {
  HMM* currHMM = nullptr;
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    currHMM = (*this->hmmVec)[hmmIdx];
    status = status + currHMM->checkForTransientStates();
  }
  std::cout << "AllCells3TrParam2DegPolyHMM::checkForTransientStates: " << status << std::endl;
  return status;
}

void AllCells3TrParam2DegPolyHMM::setBaumWelchInitGuess(gsl_vector* initGuess, int numBWIters, int numLibStarts, double libStartValue) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = 0;
  double totalTime = 0;
  double initTotalLik = 0;
  double bwTotalLik = 0;

  gsl_vector* libMultipliers = gsl_vector_alloc(numLibStarts);
  if(numLibStarts == 2) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
  }
  else if(numLibStarts == 3) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
    //gsl_vector_set(libMultipliers, 2, 0.5);
    //gsl_vector_set(libMultipliers, 2, 0.25); // half of orig
    gsl_vector_set(libMultipliers, 2, 4);
  }
  else if(numLibStarts == 4) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
    gsl_vector_set(libMultipliers, 2, 4);
    //gsl_vector_set(libMultipliers, 3, 0.5);
    //gsl_vector_set(libMultipliers, 2, 2); // quadruple orig
    //gsl_vector_set(libMultipliers, 3, 0.125); // half of orig
    gsl_vector_set(libMultipliers, 3, 0.5);
    //gsl_vector_set(libMultipliers, 0, 0.5);
    //gsl_vector_set(libMultipliers, 1, 1);
    //gsl_vector_set(libMultipliers, 2, 2);
    //gsl_vector_set(libMultipliers, 3, 4);
  }
  else {
    gsl_vector_set_all(libMultipliers, libStartValue);
  }
  this->baumWelchParamResults = new std::vector<gsl_vector*>();

  gsl_matrix* bwLibResults = gsl_matrix_alloc(numLibStarts, this->NUM_LIBS_TO_EST); // record results from each bw run. rows: runs, cols: cells
  gsl_matrix_set_zero(bwLibResults);
  //std::vector<gsl_matrix*>* bestBwTransitionMats = new std::vector<gsl_matrix*>();
  //std::vector<gsl_vector*>* bestBwInitProbs = new std::vector<gsl_vector*>();
  std::vector<std::vector<gsl_matrix*>*>* bwTransitionMats = new std::vector<std::vector<gsl_matrix*>*>(); // one tr mat per hmm per libStart
  std::vector<std::vector<gsl_vector*>*>* bwInitProbs = new std::vector<std::vector<gsl_vector*>*>();
  std::vector<bool>* bwHasTransientStates = new std::vector<bool>();
  std::vector<double>* bwTotalLiks = new std::vector<double>();

  gsl_vector* origParamsToEst = gsl_vector_alloc(this->getNumParamsToEst());
  gsl_vector_memcpy(origParamsToEst, this->paramsToEst);
  gsl_vector* origTransitionParams = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  for(unsigned origTrIdx = 0; origTrIdx < origTransitionParams->size; origTrIdx++) {
    gsl_vector_set(origTransitionParams, origTrIdx, gsl_vector_get(origParamsToEst, this->SHARED_TRANSITION_PROB_START_IDX + origTrIdx));
  }
  //std::cout << "ORIG TRANSITION PARAMS" << std::endl;
  //printColVector(origTransitionParams);
  for(unsigned int libIdx = 0; libIdx < libMultipliers->size; libIdx++) {
    // reset to original transition params. libs set below to multiple of prev iter
    this->setAllTransition(origTransitionParams);

    // for each HMM, set the lib start for each cell to the current libStart value
    for(unsigned int hmmIdx = 0, i = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS; cellIdx++) {
        // if this is the first one, use the libRatio as the starting point. else, scale the prev iter's final lib
        if(libIdx == 0) {
          double libRatio = currHMM->calcLibScalingFactorsToTotalRatio(cellIdx);
          //currHMM->setLibScalingFactor(cellIdx, gsl_vector_get(libMultipliers, libIdx));
          currHMM->setLibScalingFactor(cellIdx, gsl_vector_get(libMultipliers, libIdx) * libRatio); // scale libRatio by multiplier
          gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, gsl_vector_get(libMultipliers, libIdx) * libRatio);
        }
        else {
          //double prevLib = currHMM->getLibScalingFactor(cellIdx); // scaling relative to consecutive runs
          double prevLib = gsl_matrix_get(bwLibResults, 0, hmmIdx * currHMM->NUM_LIBS_TO_EST + cellIdx); // scaling relative to first run
          currHMM->setLibScalingFactor(cellIdx, gsl_vector_get(libMultipliers, libIdx) * prevLib); // scale prevLib by multiplier
          gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, gsl_vector_get(libMultipliers, libIdx) * prevLib);
        }
        i++;
      }
    }

    std::cout << "Starting AllCells3TrParam2DegPolyHMM::setBaumWelchInitGuess, with libMultiplier " << gsl_vector_get(libMultipliers, libIdx) << " (RUN " << libIdx + 1 << " OF " << libMultipliers->size << ")" << std::endl;
    this->print(stdout);

    // first call baum welch on each HMM
    initTotalLik = this->getLogLikelihood();

    begin = std::chrono::steady_clock::now();
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      std::cout << "STARTING BAUM WELCH for HMM " << hmmIdx << ", WITH LIBMULTIPLIER " << gsl_vector_get(libMultipliers, libIdx) << " (RUN " << libIdx + 1 << " OF " << libMultipliers->size << ")" << std::endl;
      std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before each baum welch
      currHMM->runBaumWelch(numBWIters);
      currHMM->setUpBaumWelchLeastSquares();

      // set lib size from bw into this->paramsToEst TODO check this works for multiple cells, only tested for allIndCells as of Tue 04 May 2021 05:08:49 PM PDT
      for(int cellIdx = 0; cellIdx < currHMM->NUM_LIBS_TO_EST; cellIdx++) {
        double currLib = currHMM->getLibScalingFactor(cellIdx);
        gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + (hmmIdx * currHMM->NUM_LIBS_TO_EST) + cellIdx, currLib);
        gsl_matrix_set(bwLibResults, libIdx, (hmmIdx * currHMM->NUM_LIBS_TO_EST) + cellIdx, currLib);
      }
      //(*this->hmmVec)[hmmIdx]->saveViterbiDecodedCNA("test_bw_vitDec_" + (*this->hmmNames)[hmmIdx] + "_k" + std::to_string(this->MAX_PLOIDY) + "_run" + std::to_string(libIdx) + ".txt");

      //(*this->hmmVec)[hmmIdx]->print(stderr);
    }
    //std::cout << "Finished baum welch, checking for transient states, for libMultiplier " << gsl_vector_get(libMultipliers, libIdx) << " (RUN " << libIdx + 1 << " OF " << libMultipliers->size << ")" << std::endl;
    //this->print(stdout);
    bwTotalLik = this->getLogLikelihood();
    // rule out any runs with transient states
    double hasTransientStates = this->checkForTransientStates();
    //bool saveThisRun = false;
    /*if(gsl_isnan(hasTransientStates) && libIdx != libMultipliers->size - 1) {
      std::cout << "WARNING: detected transient states after baum welch in at least one HMM for run " << libIdx + 1 << ". NOT saving transition mat; doubling lib and rerunning baum welch." << std::endl;
      saveThisRun = false;
      //continue;
    }
    else {
      if(libIdx == libMultipliers->size - 1 && gsl_isnan(hasTransientStates)) {
        if(bestLibIdx == -1) {
          std::cout << "WARNING: detected transient states after baum welch in at least one HMM for run " << libIdx + 1 << ", but no more libMultipliers and haven't found a best lib yet, so saving transition mat" << std::endl;
          saveThisRun = true;
        }
        else {
          std::cout << "WARNING: detected transient states after baum welch in at least one HMM for run " << libIdx + 1 << ". NOT saving, and no more libMultipliers" << std::endl;
          saveThisRun = false;
        }
      }
      // no transient states detected, run with this param set if it has the best ll
      else if(bwTotalLik > bestBwLik) {
        std::cout << "No transient states detected after baum welch, and this run has a better loglikelihood (" << bwTotalLik << ") than seen so far (" << bestBwLik << "), saving transition mat from run " << libIdx + 1 << std::endl;
        saveThisRun = true;
      }
      else {
        std::cout << "No transient states detected after baum welch, and this run has a lower loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). NOT saving transition mat from run " << libIdx + 1 << std::endl;
        saveThisRun = false;
      }
    }*/ // up to v11

    /*int diffInLlThreshold = 100; // threshold for "close" ll
    // if this is the first run, just save
    if(libIdx == 0) {
      std::cout << "First run, saving run " << libIdx + 1 << " with loglikelihood (" << bwTotalLik << ")";
      if(gsl_isnan(hasTransientStates)) {
        std::cout << ", even though detected transient states." << std::endl;
      }
      else {
        std::cout << ", and no transient states detected." << std::endl;
      }
      saveThisRun = true;
    }
    // if this one is worse than both best and second best, straightforward reject
    else if(bwTotalLik < bestBwLik && bwTotalLik < secondBestBwLik) {
      std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best (" << bestBwLik << ") and second best (" << secondBestBwLik << ") seen so far. Not saving." <<  std::endl;
      saveThisRun = false;
    }
    else {
      // if this run has better ll than best two
      if(bwTotalLik > bestBwLik && bwTotalLik > secondBestBwLik) {
        // if no transient states, not sus. save
        if(!gsl_isnan(hasTransientStates)) {
          std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "), and has no transient states. Saving." <<  std::endl;
          saveThisRun = true;
          secondBestBwLik = bestBwLik; // previously best lik is now second best
        }
        // else, this one is sus due to transient states. is the current best nearby or not sus?
        else {
          // current best is nearby, keep current best
          if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold) {
            std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). BUT this has transient states and current best is nearby (" << bwTotalLik - bestBwLik << "). Not saving." <<  std::endl;
            saveThisRun = false;
            secondBestBwLik = bwTotalLik; // record this one as second best bc we're not using it, even though ll is higher
          }
          // current best is not sus, keep current best
          else if(!bestHasTransientStates) {
            std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). BUT this has transient states and current best does not have transient states. Not saving." <<  std::endl;
            saveThisRun = false;
            secondBestBwLik = bwTotalLik; // record this one as second best bc we're not using it, even though ll is higher
          }
          else {
            std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). BUT this has transient states BUT current best is far (" << bwTotalLik - bestBwLik << ") and current best has transient states. Saving." <<  std::endl;
            saveThisRun = true;
            secondBestBwLik = bestBwLik; // previously best lik is now second best
          }
        }
      }
      // else, this run has a worse ll than top, but is still better than second best. should current top run be disqualified?
      else {
        // if top run is not sus, don't replace
        if(!bestHasTransientStates) {
          std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). Current best has no transient states. Not saving." <<  std::endl;
          saveThisRun = false;
          // if this run is better than the current second best, then record it
          if(bwTotalLik > secondBestBwLik) {
            secondBestBwLik = bwTotalLik;
          }
        }
        // else, current best is sus. should we replace with this one?
        else {
          // this one is not sus, replace
          if(!gsl_isnan(hasTransientStates)) {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). BUT current best has transient states and this one does not. Saving." <<  std::endl;
            saveThisRun = true;
            secondBestBwLik = bestBwLik; // despite being worse in ll, the prev best is now second
          }
          // this one is nearby, replace
          else if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold) {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). Even though this has transient states, current best is nearby (" << bwTotalLik - bestBwLik << "). Saving." <<  std::endl;
            saveThisRun = true;
            secondBestBwLik = bestBwLik; // despite being worse in ll, the prev best is now second
          }
          // else, don't save
          else {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). This has transient states and current best is far (" << bwTotalLik - bestBwLik << "). Not saving." <<  std::endl;
            saveThisRun = false;
            // if this run is better than the current second best, then record it
            if(bwTotalLik > secondBestBwLik) {
              secondBestBwLik = bwTotalLik;
            }
          }
        }
      }
    }*/
    /*else {
      // if this run has a better ll
      if(bwTotalLik > bestBwLik) {
        // if this run has a better lik and has no transient states, then is straightforwardly better. save
        if(!gsl_isnan(hasTransientStates)) {
          std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "), and has no transient states. Saving." <<  std::endl;
          saveThisRun = true;
        }
        // else, this run has a better lik but has transient states
        else {
          // if diff is small (<100 ll units; ie this one is a little bit better), then this lib is prob wrong (ie lib is prob est to be 0.5 instead of 1). do not save
          if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold) {
            std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). But, has transient states AND difference in loglikelihood is small (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
            saveThisRun = false;
          }
          // else, has transient states but is much much better
          else {
            // if prev best has transient states, and this one has transient states, but is much much better, this lib is prob better. save
            if(bestHasTransientStates) {
              std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). But, has transient states AND difference in loglikelihood is large (" << bwTotalLik - bestBwLik << ") AND previous best had transient states. Saving." <<  std::endl;
              saveThisRun = true;
            }
            // if prev best has no transient states, and this one has transient states, but is much much better, this lib is suspicious. don't save
            else {
              std::cout << "Run " << libIdx + 1 << " has a better loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). But, has transient states AND difference in loglikelihood is large (" << bwTotalLik - bestBwLik << ") AND previous best did NOT have transient states. NOT saving." <<  std::endl;
              saveThisRun = false;
            }
          }
        }
      }
      // this is a worse ll. But, if difference is small, might be correct
      else {
        // if prev best had transient states, prev best might be wrong and need to be replaced
        if(bestHasTransientStates) {
          // if diff in ll is small && if has no transient states, this run is probably right. save
          if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold){
            if(!gsl_isnan(hasTransientStates)) {
              std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << "). But, this run has no transient states AND previous best has transient states AND difference in loglikelihood is small (" << bwTotalLik - bestBwLik << "). Saving." <<  std::endl;
              saveThisRun = true;
            }
            // if diff in ll is small && if has transient states, there's no compelling reason to say this run is better. don't save.
            else {
          //else if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold && gsl_isnan(hasTransientStates)) {
              std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND this run has transient states AND previous best has transient states AND difference in loglikelihood is small (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
              saveThisRun = false;
            }
          }
          // else, large diff
          else {
            // if diff in ll is large BUT prev run is sus and this one isn't, this run is probably right. save.
            if(!gsl_isnan(hasTransientStates)) {
              std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND this run has NO transient states AND previous best has transient states AND difference in loglikelihood is large (" << bwTotalLik - bestBwLik << "). Saving." <<  std::endl;
              saveThisRun = true;
            }
            // else, diff in ll is large AND prev run is sus AND this one is too, this run is probably wrong. don't save.
            else {
              std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND this run has transient states AND previous best has transient states AND difference in loglikelihood is large (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
              saveThisRun = false;
            }
          }
        }
        // else, prev didn't have transient states, prob isn't wrong
        else {
          // if diff in ll is small (ie this one is a little bit worse) && if has no transient states, then this run is prob wrong. save
          if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold && !gsl_isnan(hasTransientStates)) {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND this run has no transient states AND previous best has no transient states AND difference in loglikelihood is small (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
            saveThisRun = false;
          }
          // if diff in ll is small (ie this one is a little bit worse) && if has transient states, then this run is prob wrong. don't save
          else if(std::abs(bwTotalLik - bestBwLik) < diffInLlThreshold && gsl_isnan(hasTransientStates)) {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND this run has transient states AND previous best has no transient states AND difference in loglikelihood is small (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
            saveThisRun = false;
          }
          // if diff in ll is large (ie this one is much worse), then this run is prob wrong. don't save
          else {
            std::cout << "Run " << libIdx + 1 << " has a worse loglikelihood (" << bwTotalLik << ") than best seen so far (" << bestBwLik << ") AND previous best has no transient states AND difference in loglikelihood is large (" << bwTotalLik - bestBwLik << "). NOT saving." <<  std::endl;
            saveThisRun = false;
          }
        }
      }
    }*/
    /*if(saveThisRun) {
      bestBwLik = bwTotalLik;
      bestLibIdx = libIdx;
      bestHasTransientStates = gsl_isnan(hasTransientStates);
      bestBwTransitionMats->clear();
      bestBwInitProbs->clear();
      for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
        HMM* currHMM = (*this->hmmVec)[hmmIdx];
        gsl_matrix* currTransMat = currHMM->getTransition();
        gsl_matrix* transitionMatToSave = gsl_matrix_alloc(currTransMat->size1, currTransMat->size2);
        gsl_matrix_memcpy(transitionMatToSave, currTransMat);
        bestBwTransitionMats->push_back(transitionMatToSave);

        gsl_vector* currInitProb = currHMM->getInitProb();
        gsl_vector* initProbToSave = gsl_vector_alloc(currInitProb->size);
        gsl_vector_memcpy(initProbToSave, currInitProb);
        bestBwInitProbs->push_back(initProbToSave);
      }
    }*/
    //(*bwTransitionMats)[libIdx] = new std::vector<gsl_matrix*>();
    //(*bwInitProbs)[libIdx] = new std::vector<gsl_vector*>();
    bwTransitionMats->push_back(new std::vector<gsl_matrix*>());
    bwInitProbs->push_back(new std::vector<gsl_vector*>());
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      gsl_matrix* currTransMat = currHMM->getTransition();
      gsl_matrix* transitionMatToSave = gsl_matrix_alloc(currTransMat->size1, currTransMat->size2);
      gsl_matrix_memcpy(transitionMatToSave, currTransMat);
      (*bwTransitionMats)[libIdx]->push_back(transitionMatToSave);

      gsl_vector* currInitProb = currHMM->getInitProb();
      gsl_vector* initProbToSave = gsl_vector_alloc(currInitProb->size);
      gsl_vector_memcpy(initProbToSave, currInitProb);
      (*bwInitProbs)[libIdx]->push_back(initProbToSave);
    }
    bwHasTransientStates->push_back(gsl_isnan(hasTransientStates));
    bwTotalLiks->push_back(bwTotalLik);


    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    totalTime = elapsedSec;
    printf("BAUM WELCH TOTAL LOGLIKELIHOOD FOUND: %.40f\n", bwTotalLik);
    printf("CHANGE IN BAUM WELCH TOTAL LIKELIHOOD: %.40f\n", bwTotalLik - initTotalLik);
    printf("BAUM WELCH TOTAL TIME (sec): %.10f\n", elapsedSec);
    printf("DONE WITH BAUM WELCH TOTAL (AllCells3TrParam2DegPolyHMM):\n\n");
    this->print(stdout);
    std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before starting least squares bfgs
  } // end iterating over libMultipliers


    // figure out which libMultiplier result to send
    int bestLibIdx = -1;
    int secondBestLibIdx = -1;
    double bestBwLik = -std::numeric_limits<double>::max();
    double secondBestBwLik = -std::numeric_limits<double>::max(); // not strictly second best by ll, could be higher in ll but have transient states/be suspicious
    for(unsigned int bwIdx = 0; bwIdx < bwTotalLiks->size(); bwIdx++) {
      bwTotalLik = (*bwTotalLiks)[bwIdx];
      //std::cout << "checking bwIdx " << bwIdx << ", bwTotalLik = " << bwTotalLik << ", bestLibIdx so far: " << bestLibIdx << " (" << bestBwLik << "), secondBestLibIdx: " << secondBestLibIdx << " (" << secondBestBwLik << ")" << std::endl;
      // bwTotalLik > bestBwLik > secondBestBwLik
      // bwTotalLik is better than both bests seen so far, this is now the best
      if(bwTotalLik > bestBwLik && bwTotalLik > secondBestBwLik) {
        secondBestBwLik = bestBwLik;
        secondBestLibIdx = bestLibIdx;
        bestBwLik = bwTotalLik;
        bestLibIdx = bwIdx;
      }
      // bestBwLik > bwTotalLik > secondBestBwLik
      // bwTotalLik is in the middle, new second best
      else if(bwTotalLik < bestBwLik && bwTotalLik > secondBestBwLik) {
        secondBestBwLik = bwTotalLik;
        secondBestLibIdx = bwIdx;
      }
      // bestBwLik > secondBestBwLik > bwTotalLik
      // bwTotalLik is worse, skip over
      else {
        continue;
      }
    }
    bool bestHasTransientStates = (*bwHasTransientStates)[bestLibIdx];
    bool secondBestHasTransientStates = (*bwHasTransientStates)[secondBestLibIdx];

    std::cout << "Run " << bestLibIdx + 1 << " has the best loglikelihood (" << bestBwLik << "), followed by run " << secondBestLibIdx + 1 << " (" << secondBestBwLik << "), diff (" << std::abs(secondBestBwLik - bestBwLik) << ")." << std::endl;
    // if bestHasTransientStates, might get replaced
    int diffInLlThreshold = 100; // threshold for "close" ll
    if(bestHasTransientStates) {
      std::cout << "Best loglik has transient states, ";
      // if second best is nearby or is not suspicious, replace best
      if(std::abs(secondBestBwLik - bestBwLik) < diffInLlThreshold || !secondBestHasTransientStates) {
        if(std::abs(secondBestBwLik - bestBwLik) < diffInLlThreshold) {
          std::cout << "diff is small, ";
        }
        if(!secondBestHasTransientStates) {
          std::cout << "second best does NOT have transient states, ";
        }
        std::cout << "so replacing with second best." << std::endl;
        bestLibIdx = secondBestLibIdx;
        bestBwLik = secondBestBwLik;
        bestHasTransientStates = secondBestHasTransientStates;
      }
      // else, no compelling reason to replace
      else {
        std::cout << "but second best ";
        if(std::abs(secondBestBwLik - bestBwLik) >= diffInLlThreshold) {
          std::cout << "is far, ";
        }
        if(secondBestHasTransientStates) {
          std::cout << "has transient states, ";
        }
        std::cout << "so NOT replacing with second best." << std::endl;
      }
    }
    // else, just keep best
    else {
      std::cout << "Best loglik does NOT have transient states, using best." << std::endl;
    }
    std::vector<gsl_matrix*>* bestBwTransitionMats = (*bwTransitionMats)[bestLibIdx];
    std::vector<gsl_vector*>* bestBwInitProbs = (*bwInitProbs)[bestLibIdx];




    // send the best ll transition mats and lib sizes to least squares
    for(unsigned int hmmIdx = 0, i = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      gsl_matrix* savedTransMat = (*bestBwTransitionMats)[hmmIdx];
      currHMM->setTransition(savedTransMat);
      gsl_vector* savedInitProb = (*bestBwInitProbs)[hmmIdx];
      currHMM->setInitProb(savedInitProb);
      for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS; cellIdx++) {
        double storedLib = gsl_matrix_get(bwLibResults, bestLibIdx, hmmIdx * currHMM->NUM_LIBS_TO_EST + cellIdx);
        currHMM->setLibScalingFactor(cellIdx, storedLib);
        gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, storedLib);
        i++;
      }
    }
    std::cout << "RESETTING TO BEST TRANSITION MATRICES AND LIB SIZE SCALING FACTORS, FROM RUN " << bestLibIdx + 1 << ":" << std::endl;
    this->print(stdout);


    /*
     * check if lib size needs to be doubled.
     * lib size needs to be doubled if the viterbi decoding has transient states
     * (ex state 1 is only seen in 0->1->2 or 2->1->0 transitions)
     */
    ////double hasTransientStates = this->checkForTransientStates();
    //hasTransientStates = this->checkForTransientStates();
    //if(gsl_isnan(hasTransientStates) && libIdx != libMultipliers->size - 1) {
    //  std::cout << "WARNING: detected transient states after baum welch in at least one HMM for run " << libIdx + 1 << ". Doubling lib and rerunning baum welch." << std::endl;
    //  continue;
    //}
    //else {
    //  // no transient states detected, run with this param set
    //  if(libIdx == libMultipliers->size - 1 && gsl_isnan(hasTransientStates)) {
    //    std::cout << "WARNING: detected transient states after baum welch in at least one HMM for run " << libIdx + 1 << ", but no more libMultipliers, so continuing on to least squares" << std::endl;
    //  }
    //  else {
    //    std::cout << "No transient states detected after baum welch, continuing to least squares for run " << libIdx + 1 << std::endl;
    //  }

    // then solve overdetermined system of equations for each HMM
    this->doLeastSquares(initGuess);

    double lsTotalLik = this->getLogLikelihood();
    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    totalTime += elapsedSec;
    printf("CHANGE IN LEAST SQUARES LIKELIHOOD FROM STARTING POINT: %.40f\n", lsTotalLik - initTotalLik);
    printf("CHANGE IN LEAST SQUARES LIKELIHOOD FROM BAUM WELCH: %.40f\n", lsTotalLik - bestBwLik);
    printf("BAUM WELCH WITH LEAST SQUARES TOTAL TIME (sec): %.10f\n", totalTime);
    printf("DONE WITH AllCells3TrParam2DegPolyHMM BAUM WELCH LEAST SQUARES\n");

    // if optim failed (as determined by a uniform transition matrix), try again. else, save params
    double allValidTransitionMatrices = this->checkStateValidity();
    if(gsl_isnan(allValidTransitionMatrices)) {
      std::cout << "WARNING: at least one HMM had an invalid transition matrix from this baum welch + least squares parameter set:" << std::endl;
      //printColVector(paramsToSave);
      this->print(stdout);
      //std::cout << "Attempting to run baum welch again, run " << libIdx + 1 << " failed after least squares" << std::endl;
      //continue;
    }
    gsl_vector* paramsToSave = gsl_vector_alloc(initGuess->size);
    gsl_vector_memcpy(paramsToSave, initGuess);
    this->baumWelchParamResults->push_back(paramsToSave);

    this->print(stdout);
    std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file after ending least squares bfgs

    //  break;
    //}
  //} // end iterating over libMultipliers


  if(this->baumWelchParamResults->size() == 0) {
    std::cout << "WARNING: All baum welch + least squares runs failed. Using original params" << std::endl;
    std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file after ending least squares bfgs
    gsl_vector* origParamsToEstCopy = gsl_vector_alloc(origParamsToEst->size);
    gsl_vector_memcpy(origParamsToEstCopy, origParamsToEst);
    this->baumWelchParamResults->push_back(origParamsToEstCopy);
    this->setParamsToEst(origParamsToEstCopy);
  }
  std::cout << "The baum welch parameter sets are:" << std::endl;
  for(unsigned int i = 0; i < this->baumWelchParamResults->size(); i++) {
    if((*this->baumWelchParamResults)[i] != nullptr) {
      printColVector((*this->baumWelchParamResults)[i]);
    }
  }
  gsl_vector_memcpy(initGuess, (*this->baumWelchParamResults)[0]); // save first set of params into initGuess

  // clean up
  //for(unsigned int i = 0; i < unsortedBaumWelchParamResults->size(); i++) {
  //  gsl_vector_free((*unsortedBaumWelchParamResults)[i]);
  //}
  //delete unsortedBaumWelchParamResults; // delete vector, but keep contents (which are now stored in this->baumWelchParamResults)
  //for(unsigned int i = 0; i < unsortedAveragePloidies->size(); i++) {
  //  gsl_matrix_free((*unsortedAveragePloidies)[i]);
  //}
  //delete unsortedAveragePloidies;
  bestBwTransitionMats->clear();
  bestBwInitProbs->clear();
  gsl_vector_free(libMultipliers);
  gsl_vector_free(origParamsToEst);
  gsl_matrix_free(bwLibResults);
}
double AllCells3TrParam2DegPolyHMM::doLeastSquares(gsl_vector* initGuess, int maxNumAttempts) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double overallStatus = GSL_NAN;
  double status = GSL_CONTINUE;
  double initSumSqResid = 0;
  double currSumSqResid = 0;
  double prevSumSqResid = 0;
  double deltaSumSqResid = 0;
  double totalSumSqResid = 0;
  double finalSumSqResid = 0;
  int iter = 0;
  int countTooClose = 0;
  std::chrono::steady_clock::time_point lsBegin;
  std::chrono::steady_clock::time_point lsEnd;
  double lsElapsedTime = 0;

  /*
   * structure of varsToEst
   * [b]
   * [L]
   * [pair0, t1]
   * [pair0, t2]
   * [pair0, t3]
   * ...
   * [pairN, t3]
   */
  gsl_vector* varsToEst = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  gsl_vector* varsToEst_probSpace = gsl_vector_alloc(varsToEst->size); // for printing in prob space

  for(int lsAttempt = 1; lsAttempt <= maxNumAttempts; lsAttempt++) {
    std::cout << "STARTING LEAST SQUARES, (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ")" << std::endl;

    if(lsAttempt == 1) {
      gsl_vector_set(varsToEst_probSpace, 0, 0.01); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 500); // lambda
      //gsl_vector_set(varsToEst_probSpace, 1, 15); // lambda
      //gsl_vector_set(varsToEst_probSpace, 0, 0.001); // beta
      //gsl_vector_set(varsToEst_probSpace, 1, 1500); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        //gsl_vector_set(varsToEst_probSpace, 2, 0.5); // t (for single independent cells)
        //gsl_vector_set(varsToEst_probSpace, 2, 0.001); // t (for single independent cells)
        gsl_vector_set(varsToEst_probSpace, 2, 0.01); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.01); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.01); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.01); // t3
        }
      }
      /*gsl_vector_set(varsToEst_probSpace, 0, 0.025); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 10); // lambda
      for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
        gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.025); // t1
        gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.025); // t2
        gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.025); // t3
      }*/
      //gsl_vector_set(varsToEst_probSpace, 1, 100);
    }
    else if(lsAttempt == 2) {
      gsl_vector_set(varsToEst_probSpace, 0, 0.1); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 100); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        gsl_vector_set(varsToEst_probSpace, 2, 0.001); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.001); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.001); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.001); // t3
        }
      }
    }
    else {
      gsl_vector_set(varsToEst_probSpace, 0, 0.05); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 1000); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        gsl_vector_set(varsToEst_probSpace, 2, 0.05); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.05); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.05); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.05); // t3
        }
      }
    }

    this->baumWelchLeastSquares_convertProbToParam(varsToEst, varsToEst_probSpace);

    gsl_multimin_function_fdf my_func; // which function to minimize
    my_func.df = &baumWelchLeastSquares_df; // gradient of function
    my_func.fdf = &baumWelchLeastSquares_fdf; // how to set both f and df

    my_func.f = &baumWelchLeastSquares_f; // function itself
    my_func.params = this;
    my_func.n = varsToEst->size;

    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2; // which algorithm to use
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, my_func.n);
    gsl_multimin_fdfminimizer_set(s, &my_func, varsToEst, 1e-1, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance

    status = GSL_CONTINUE;
    currSumSqResid = 0;
    deltaSumSqResid = 0;
    totalSumSqResid = 0;
    iter = 0;
    countTooClose = 0;
    lsElapsedTime = 0;

    initSumSqResid = baumWelchLeastSquares_f(varsToEst, this);
    prevSumSqResid = initSumSqResid;
    printf("BW LEAST SQUARES BFGS INITIAL currSumSqResid: %.20f\n", initSumSqResid);

    while(status == GSL_CONTINUE) {
      lsBegin = std::chrono::steady_clock::now();

      status = gsl_multimin_fdfminimizer_iterate(s); // do one iter
      currSumSqResid = gsl_multimin_fdfminimizer_minimum(s);

      lsEnd = std::chrono::steady_clock::now();
      lsElapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(lsEnd - lsBegin).count() / 1000000.0;
      deltaSumSqResid = std::abs(prevSumSqResid - currSumSqResid);
      totalSumSqResid += deltaSumSqResid;
      prevSumSqResid = currSumSqResid;

      printf("ON BW LEAST SQUARES BFGS ITER %d, currSumSqResid %.20f, time elapsed (sec) %.5f, total change in currSumSqResid %.5f\n", iter, currSumSqResid, lsElapsedTime, totalSumSqResid);

      if(status) {
        std::cout << "STATUS IS: " << gsl_strerror(status) << std::endl;
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s));
        printRowVector(varsToEst_probSpace);
        break;
      }
      //status = gsl_multimin_test_gradient(s->gradient, 1e-4);
      //status = gsl_multimin_test_gradient(s->gradient, 1e-6);
      status = gsl_multimin_test_gradient(s->gradient, 1e-8);

      if(status == GSL_SUCCESS) {
        printf("SUCCESS: BW+LS minimum found at: ");
        //printRowVector(varsToEst);
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s));
        printRowVector(varsToEst_probSpace);
      }
      if(iter >= 500) {
      //if(iter >= 2) {
      //if(iter >= 1) {
        std::cout << "STATUS IS: minimum not found in 500 iterations" << std::endl;
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s)); // save anyway
        printRowVector(varsToEst_probSpace);
        break;
      }

      if(deltaSumSqResid < 1e-5) {
      //if(deltaSumSqResid < 1e-8) {
        countTooClose++;
        if(countTooClose >= 3) {
          std::cout << "STATUS IS: converged by consecutive small change in sumSqResid" << std::endl;
          this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s)); // save anyway
          printRowVector(varsToEst_probSpace);
          break;
        }
      } else {
        countTooClose = 0;
      }

      iter++;
    }
    gsl_multimin_fdfminimizer_free(s);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    printf("LEAST SQUARES TIME (sec): %.10f\n", elapsedSec);

    // if ls succeeded (>0 change in currSumSqResid and no uniform matrix), then break. else, try again
    finalSumSqResid = currSumSqResid;
    deltaSumSqResid = std::abs(initSumSqResid - finalSumSqResid);
    double allValidTransitionMatrices = this->checkStateValidity();
    if(deltaSumSqResid > 0 && !gsl_isnan(allValidTransitionMatrices)) {
      std::cout << "LEAST SQUARES (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ") SUCCEEDED, BREAKING\n" << std::endl;
      overallStatus = 0;

      // save results. delegate to subclass since depends on number of cells to save libs
      this->saveBaumWelchEstsIntoParamsToEst(varsToEst_probSpace, initGuess);
      double lsTotalLik = this->getLogLikelihood();
      printf("LEAST SQUARES LOGLIKELIHOOD FOUND: %.40f\n", lsTotalLik);

      break;
    }
    else {
      std::cout << "LEAST SQUARES (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ") FAILED due to ";
      if(gsl_isnan(allValidTransitionMatrices)) {
        std::cout << "at least one invalid transition matrix";
      }
      else {
        std::cout << "no change in currSumSqResid";
      }
      std::cout << ". Trying again\n" << std::endl;
    }
  }

  return overallStatus;
}
double AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(const gsl_vector* v, void* params) {
  AllCells3TrParam2DegPolyHMM* hmm = (AllCells3TrParam2DegPolyHMM*) params;
  double sumSqResid = hmm->baumWelchLeastSquares_calcSumSqResid(v); // delegate to subclass, since depends on param ordering
  return sumSqResid;
}
void AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_df(const gsl_vector* v, void* params, gsl_vector* df) {
  // approx using finite difference method:
  // f'(x) = [f(x + h) - f(x)] / h
  //double h = 1e-8; // step size from testOptim
  //double h = 1e-6; // step size from Optimizable
  //double h = 1e-5;
  double h = 1e-4; // bwh4
  //double h = 1e-3;
  //double h = 1e-2;
  /*double f_x = AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(v, params);
  if(gsl_isnan(f_x)) {
    std::cerr << "AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f is nan" << std::endl;
    gsl_vector_set_all(df, GSL_NAN);
    return;
  }*/

  double f_x_ph = 0; // f(x+h)
  double f_x_mh = 0; // f(x-h)
  double x = 0;
  double derivApprox = 0;
  gsl_vector_set_zero(df);
  gsl_vector* currVec = gsl_vector_alloc(v->size);
  gsl_vector_set_zero(currVec);
  for(unsigned int i = 0; i < v->size; i++) {
    gsl_vector_memcpy(currVec, v);
    //x_h = gsl_vector_get(v, i) + h;
    x = gsl_vector_get(v, i);
    gsl_vector_set(currVec, i, x + h);
    f_x_ph = AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(currVec, params);

    if(gsl_isnan(f_x_ph)) {
      derivApprox = GSL_NAN;
    } else {
      //derivApprox = (f_x_ph - f_x) / h; // forward difference

      gsl_vector_set(currVec, i, x - h);
      f_x_mh = AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(currVec, params);
      if(gsl_isnan(f_x_mh)) {
        derivApprox = GSL_NAN;
      } else {
        derivApprox = (f_x_ph - f_x_mh) / (2*h); // central difference*/
      }
    }
    gsl_vector_set(df, i, derivApprox);
  }

  // BFGS printing
  gsl_vector* probs = gsl_vector_alloc(v->size);
  ((AllCells3TrParam2DegPolyHMM*) params)->baumWelchLeastSquares_convertParamToProb(probs, v);
  fprintf(stderr, "%.20f \t", AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(v, params));
  fprintf(stderr, "|\t");
  for(unsigned int idx = 0; idx < probs->size; idx++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(probs, idx));
  }
  fprintf(stderr, "|\t");
  for(unsigned int idx = 0; idx < v->size; idx++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(v, idx));
  }
  fprintf(stderr, "|\t");
  for(unsigned int idx = 0; idx < df->size; idx++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(df, idx));
  }
  fprintf(stderr, "\n");
  gsl_vector_free(probs);
}
void AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
  *f = baumWelchLeastSquares_f(v, params);
  baumWelchLeastSquares_df(v, params, df);
}

AllCells3TrParam2DegPolyHMM* AllCells3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  AllCells3TrParam2DegPolyHMM* bestGuessOptim = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, verbose);
  return bestGuessOptim;
}

double AllCells3TrParam2DegPolyHMM::getLogLikelihood() {
  return this->runForwardAlg();
}

double AllCells3TrParam2DegPolyHMM::getViterbiLogLikelihood() {
  double vitLL = 0;
  double currLL = 0;
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    currLL = (*itr)->getViterbiLogLikelihood();
    vitLL += currLL;
  }
  return vitLL;
}

/*
 * run forward alg on each pair of cells, summing loglikelihood over all pairs and returning that
 */
double AllCells3TrParam2DegPolyHMM::runForwardAlg() {
  double totalLogLikelihood = 0;
  double currLL = 0;
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    currLL = (*itr)->runForwardAlg();
    totalLogLikelihood += currLL;
  }
  return totalLogLikelihood;
}

/*
 * call each HMM's viterbiDecode
 */
void AllCells3TrParam2DegPolyHMM::viterbiDecodeAll() {
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    (*itr)->viterbiDecode();
    // TODO HERE write counted transitions matrix to file here
  }
}

/*
 * call each HMM's saveViterbiDecodedCNA method. Each HMM's output is numbered from 0
 */
void AllCells3TrParam2DegPolyHMM::saveAllViterbiDecodedCNA(std::string filename) {
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    std::string hmmName = (*this->hmmNames)[hmmIdx];
    boost::replace_all(hmmName, ",", "__");
    (*this->hmmVec)[hmmIdx]->saveViterbiDecodedCNA(filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".viterbiDecoded");
  }

  std::vector<gsl_vector*>* BFGSParamResults = this->getBFGSParamResults();
  if(BFGSParamResults != nullptr) {
    gsl_vector* bestParams = gsl_vector_alloc(this->getNumParamsToEst());
    gsl_vector_memcpy(bestParams, this->getParamsToEst());

    for(unsigned int i = 1; i <= BFGSParamResults->size(); i++) {
      this->setParamsToEst((*BFGSParamResults)[i - 1]);
      std::cout << "SAVING RESULTS FOR i: " << i << std::endl;
      this->print(stdout);
      for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
        std::string hmmName = (*this->hmmNames)[hmmIdx];
        boost::replace_all(hmmName, ",", "__");
        // redo viterbi decoding under these params before saving
        (*this->hmmVec)[hmmIdx]->viterbiDecode();
        (*this->hmmVec)[hmmIdx]->saveViterbiDecodedCNA(filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + "__run" + std::to_string(i) + ".viterbiDecoded");
      }
    }
    // reset to previously best params
    this->setParamsToEst(bestParams);
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      (*this->hmmVec)[hmmIdx]->viterbiDecode();
    }
    gsl_vector_free(bestParams);
  }
}
