#include "HMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
//HMM::HMM(std::vector<DepthPair*>* depths, int numTransitionParamsToEst, int numCells, int numBranches, int maxPloidy) :
HMM::HMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int numTransitionParamsToEst, int numCells, int numBranches, int maxPloidy, int numFixedTrParams, int numFixedLibs) :
                                        MAX_PLOIDY(maxPloidy),
                                        NUM_CELLS(numCells),
                                        NUM_LIBS_TO_EST(numCells - numFixedLibs),
                                        NUM_TRANSITION_PARAMS_TO_EST(numTransitionParamsToEst),
                                        NUM_BRANCH_LENGTHS_TO_EST(numBranches),
                                        NUM_FIXED_LIBS(numFixedLibs),
                                        NUM_FIXED_TRANSITION_PARAMS(numFixedTrParams),
                                        LIB_SIZE_SCALING_FACTOR_START_IDX(0),
                                        //TRANSITION_PROB_START_IDX(numCells),
                                        TRANSITION_PROB_START_IDX(numCells - numFixedLibs), // libs go first, estimating one for each "unfixed lib" cell
                                        BRANCH_LENGTH_START_IDX(this->TRANSITION_PROB_START_IDX + numTransitionParamsToEst),
                                        FIXED_TRANSITION_PROB_START_IDX(numFixedLibs) // if no libs are fixed, then this starts at 0
                                        {
  // member variables
  // alphabet consists of 0..maxObserved. assumes only ints
  this->depths = depths;
  int maxObserved = -1;
  int currTumorMax = -1;
  int currDiploidMax = -1;
  for(std::vector<DepthPair*>::iterator it = this->depths->begin(); it != this->depths->end(); ++it) {
    currTumorMax = (*it)->maxTumorDepth;
    if(currTumorMax > maxObserved) {
      maxObserved = currTumorMax;
    }
    currDiploidMax = (*it)->maxDiploidDepth;
    if(currDiploidMax > maxObserved) {
      maxObserved = currDiploidMax;
    }
  }
  std::vector<int> fullAlphabet(maxObserved + 1); // add one to make inclusive
  std::iota(fullAlphabet.begin(), fullAlphabet.end(), 0);
  this->alphabet = new std::set<int>(fullAlphabet.begin(), fullAlphabet.end());

  this->states = createStateSpace(this->NUM_CELLS, this->MAX_PLOIDY);
  this->adjBinLabels = createStateSpace(2, this->MAX_PLOIDY);
  //this->transitionStrings = nullptr;
  this->rateMatrixStrings = nullptr;

  this->transition = gsl_matrix_alloc(this->states->size(), this->states->size());
  gsl_matrix_set_zero(this->transition);

  //this->rateMatrixQ = gsl_matrix_alloc(this->states->size(), this->states->size());
  //this->timeDepMatrixP = gsl_matrix_alloc(this->states->size(), this->states->size());
  // contiuous time matrices are looking at adj bins, so the dimensions are (k+1,k+1) x (k+1,k+1)
  int contTimeMatSize = (this->MAX_PLOIDY + 1) * (this->MAX_PLOIDY + 1);
  this->rateMatrixQ    = gsl_matrix_alloc(contTimeMatSize, contTimeMatSize);
  this->timeDepMatrixP = gsl_matrix_alloc(contTimeMatSize, contTimeMatSize);
  gsl_matrix_set_zero(this->rateMatrixQ);
  gsl_matrix_set_zero(this->timeDepMatrixP);

  this->initProb = gsl_vector_alloc(this->states->size());
  gsl_vector_set_zero(this->initProb);

  this->generator = nullptr;

  this->optimSuccess = false;

  // intermediates
  this->allocIntermediates();
  //this->setMeanVarianceFn(gsl_vector_alloc(3)); // intercept, slope, first order coef TODO move this into subclass
  this->meanVarianceCoefVec = nullptr; // Thu 23 Jul 2020 05:44:08 PM PDT setMeanVarianceFn changes ptrs

  // params to estimate with BFGS: a library scaling factor for each cell and the transition params
  //std::cout << "IN HMM CTOR: " << this->NUM_CELLS << ", " << numFixedLibs << ", " << this->NUM_TRANSITION_PARAMS_TO_EST << ", " << numFixedTrParams << ", " << this->NUM_BRANCH_LENGTHS_TO_EST << std::endl;
  //this->paramsToEst = gsl_vector_alloc(this->NUM_CELLS - numFixedLibs + this->NUM_TRANSITION_PARAMS_TO_EST - numFixedTrParams + this->NUM_BRANCH_LENGTHS_TO_EST);
  this->paramsToEst = gsl_vector_alloc(this->NUM_CELLS - numFixedLibs + this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  //gsl_vector_set_zero(this->paramsToEst);
  gsl_vector_set_all(this->paramsToEst, 1.0);
  if(fixedParams == nullptr) {
    this->fixedParams = gsl_vector_alloc(numFixedLibs + numFixedTrParams);
    gsl_vector_set_zero(this->fixedParams);
  }
  else {
    this->fixedParams = fixedParams;
  }
  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst());
}
void HMM::allocIntermediates() {
  //std::cerr << "ALLOCATING INTERMEDIATES" << std::endl;
  this->prevForwardCol = gsl_vector_alloc(this->states->size());
  this->currForwardCol = gsl_vector_alloc(this->states->size());
  this->nextBackwardCol = gsl_vector_alloc(this->states->size());
  this->currBackwardCol = gsl_vector_alloc(this->states->size());
  this->transitionTranspose = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  int numSteadyStateEVals = this->transitionTranspose->size1;
  this->ssEval = gsl_vector_complex_alloc(numSteadyStateEVals);
  this->ssEvec = gsl_matrix_complex_alloc(numSteadyStateEVals, numSteadyStateEVals);
  this->ssEigenWorkspace = gsl_eigen_nonsymmv_alloc(numSteadyStateEVals);

  int numRateEvals = this->rateMatrixQ->size1;
  this->rateEval = gsl_vector_complex_alloc(numRateEvals);
  this->rateEvec = gsl_matrix_complex_alloc(numRateEvals, numRateEvals);
  this->rateEigenWorkspace = gsl_eigen_nonsymmv_alloc(numRateEvals);
  this->rateRealEvecMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
  this->rateDiagMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
  this->rateLUdecompMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
  this->ratePerm = gsl_permutation_alloc(numRateEvals);
  this->rateInverseMat = gsl_matrix_alloc(numRateEvals, numRateEvals);

  this->forwardMatVec = new std::vector<gsl_matrix*>();
  this->backwardMatVec = new std::vector<gsl_matrix*>();
  this->forBackMatVec = new std::vector<gsl_matrix*>();
  this->backTraceVec = new std::vector<gsl_matrix*>();
  this->scalingVecVec = new std::vector<gsl_vector*>();
  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();

  for(unsigned int i = 0; i < chrVec->size(); i++) {
    currChr = (*chrVec)[i];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
    //this->forwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
    //this->backwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
    //this->forBackMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
    //this->backTraceVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
    //this->scalingVecVec->push_back(gsl_vector_alloc(currNumWindows+1));
    this->forwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
    this->backwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
    this->forBackMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
    this->backTraceVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
    this->scalingVecVec->push_back(gsl_vector_alloc(currNumWindows));
    gsl_matrix_set_zero((*this->forwardMatVec)[i]);
    gsl_matrix_set_zero((*this->backwardMatVec)[i]);
    gsl_matrix_set_zero((*this->forBackMatVec)[i]);
    gsl_matrix_set_zero((*this->backTraceVec)[i]);
    gsl_vector_set_zero((*this->scalingVecVec)[i]);
  }

  //this->fillDiploidDepthPloidyPreCalc();
  this->numTransitionsMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_matrix_set_zero(this->numTransitionsMat);

}

/*void HMM::fillDiploidDepthPloidyPreCalc() {
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  this->diploidDepthPloidyPreCalc = new std::vector<double>(firstDepthPair->numWindows * (this->MAX_PLOIDY + 1));
  int preCalcIdx = 0;
  int currPloidy = -1;
  double currDiploidDepth = -1;
  double diploidPloidyProd = 0;
  std::vector<double>* currDiploidDepthsVec = nullptr;
  // for each chr
  for(std::vector<std::string>::iterator chrItr = firstDepthPair->chrVec->begin(); chrItr != firstDepthPair->chrVec->end(); ++chrItr) {
    currDiploidDepthsVec = (*firstDepthPair->chrToDiploidDepthMap)[*chrItr];
    // for each region in chr
    for(unsigned int dipDepthVecIdx = 0; dipDepthVecIdx < currDiploidDepthsVec->size(); dipDepthVecIdx++) {
      currDiploidDepth = (*currDiploidDepthsVec)[dipDepthVecIdx];
      // for each possible ploidy, calc ploidy * diploidDepth / 2
      for(currPloidy = 0; currPloidy <= this->MAX_PLOIDY; currPloidy++) {
        diploidPloidyProd = currPloidy * currDiploidDepth / 2.0;
        (*this->diploidDepthPloidyPreCalc)[preCalcIdx] = diploidPloidyProd;
        preCalcIdx++;
      }
    }
  }
  //for(int i = 0; i < this->diploidDepthPloidyPreCalc->size(); i++) {
  //  std::cout << (*this->diploidDepthPloidyPreCalc)[i] << ", ";
  //}
  //std::cout << std::endl;
}*/

//// deep copy ctor
//HMM::HMM(const HMM& otherHMM) :
//        Optimizable(otherHMM),
//        MAX_PLOIDY(otherHMM.MAX_PLOIDY),
//        NUM_CELLS(otherHMM.NUM_CELLS),
//        NUM_LIBS_TO_EST(otherHMM.NUM_LIBS_TO_EST),
//        NUM_TRANSITION_PARAMS_TO_EST(otherHMM.NUM_TRANSITION_PARAMS_TO_EST),
//        NUM_BRANCH_LENGTHS_TO_EST(otherHMM.NUM_BRANCH_LENGTHS_TO_EST),
//        NUM_FIXED_LIBS(otherHMM.NUM_FIXED_LIBS),
//        NUM_FIXED_TRANSITION_PARAMS(otherHMM.NUM_FIXED_TRANSITION_PARAMS),
//        LIB_SIZE_SCALING_FACTOR_START_IDX(otherHMM.LIB_SIZE_SCALING_FACTOR_START_IDX),
//        TRANSITION_PROB_START_IDX(otherHMM.TRANSITION_PROB_START_IDX),
//        BRANCH_LENGTH_START_IDX(otherHMM.BRANCH_LENGTH_START_IDX),
//        FIXED_TRANSITION_PROB_START_IDX(otherHMM.FIXED_TRANSITION_PROB_START_IDX) {
//  //std::cout << "IN HMM DEEP COPY CTOR. BRANCH_LENGTH_START_IDX: " << this->BRANCH_LENGTH_START_IDX << std::endl;
//  // member variables
//  this->depths = new std::vector<DepthPair*>(*(otherHMM.depths)); // TODO not sure correct copy ctor is called here
//  this->states = new std::vector<std::string>(*(otherHMM.states));
//  this->alphabet = new std::set<int>(*(otherHMM.alphabet));
//
//  gsl_matrix* mat = otherHMM.getTransition();
//  this->transition = gsl_matrix_alloc(mat->size1, mat->size2);
//  gsl_matrix_memcpy(this->transition, mat);
//
//  this->generator = nullptr; // TODO what's a good way to copy seed info? is it worth it?
//
//  this->optimSuccess = false;
//
//  // make a copy of transitionStrings
//  std::string* currTrStr = nullptr;
//  this->transitionStrings = new std::vector<std::vector<std::string*>*>(this->states->size());
//  for(unsigned int i = 0; i < this->states->size(); i++) {
//    (*this->transitionStrings)[i] = new std::vector<std::string*>(this->states->size());
//    for(unsigned int j = 0; j < this->states->size(); j++) {
//      currTrStr = new std::string(*(*(*otherHMM.transitionStrings)[i])[j]);
//      (*(*this->transitionStrings)[i])[j] = currTrStr;
//    }
//  }
//
//  gsl_vector* vec = otherHMM.getInitProb();
//  this->initProb = gsl_vector_alloc(vec->size);
//  gsl_vector_memcpy(this->initProb, vec);
//
//  // intermediates
//  this->forwardMatVec = new std::vector<gsl_matrix*>();
//  this->backwardMatVec = new std::vector<gsl_matrix*>();
//  this->forBackMatVec = new std::vector<gsl_matrix*>();
//  this->backTraceVec = new std::vector<gsl_matrix*>();
//  this->scalingVecVec = new std::vector<gsl_vector*>();
//  std::string currChr;
//  int currNumWindows = -1;
//  std::vector<std::string>* chrVec = this->getChrVec();
//
//  for(unsigned int i = 0; i < chrVec->size(); i++) {
//    currChr = (*chrVec)[i];
//    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
//    //this->forwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
//    //this->backwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
//    //this->forBackMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
//    //this->backTraceVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows+1));
//    this->forwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
//    this->backwardMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
//    this->forBackMatVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
//    this->backTraceVec->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
//    this->scalingVecVec->push_back(gsl_vector_alloc(currNumWindows+1));
//    gsl_matrix_set_zero((*this->forwardMatVec)[i]);
//    gsl_matrix_set_zero((*this->backwardMatVec)[i]);
//    gsl_matrix_set_zero((*this->forBackMatVec)[i]);
//    gsl_matrix_set_zero((*this->backTraceVec)[i]);
//    gsl_vector_set_zero((*this->scalingVecVec)[i]);
//  }
//
//  this->prevForwardCol = gsl_vector_alloc(this->states->size());
//  this->currForwardCol = gsl_vector_alloc(this->states->size());
//  this->nextBackwardCol = gsl_vector_alloc(this->states->size());
//  this->currBackwardCol = gsl_vector_alloc(this->states->size());
//  this->transitionTranspose = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
//  int numSteadStateEVals = this->transitionTranspose->size1;
//  this->ssEval = gsl_vector_complex_alloc(numSteadStateEVals);
//  this->ssEvec = gsl_matrix_complex_alloc(numSteadStateEVals, numSteadStateEVals);
//  gsl_vector_complex_memcpy(this->ssEval, otherHMM.eval);
//  gsl_matrix_complex_memcpy(this->ssEvec, otherHMM.evec);
//  this->ssEigenWorkspace = gsl_eigen_nonsymmv_alloc(numSteadStateEVals);
//
//  int numRateEvals = this->rateMatrixQ->size1;
//  this->rateEval = gsl_vector_complex_alloc(numRateEvals);
//  this->rateEvec = gsl_matrix_complex_alloc(numRateEvals, numRateEvals);
//  this->rateEigenWorkspace = gsl_eigen_nonsymmv_alloc(numRateEvals);
//  this->rateRealEvecMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
//  this->rateDiagMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
//  this->rateLUdecompMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
//  this->ratePerm = gsl_permutation_alloc(numRateEvals);
//  this->rateInverseMat = gsl_matrix_alloc(numRateEvals, numRateEvals);
//
//  this->meanVarianceCoefVec = gsl_vector_alloc(otherHMM.meanVarianceCoefVec->size);
//  gsl_vector_memcpy(this->meanVarianceCoefVec, otherHMM.meanVarianceCoefVec);
//
//  //this->diploidDepthPloidyPreCalc = new std::vector<double>(*(otherHMM.diploidDepthPloidyPreCalc));
//
//  /*// params to estimate with BFGS
//  this->paramsToEst = gsl_vector_alloc(otherHMM.getNumParamsToEst());
//  gsl_vector_memcpy(this->paramsToEst, otherHMM.paramsToEst);
//  this->fixedParams = gsl_vector_alloc(otherHMM.fixedParams->size);
//  gsl_vector_memcpy(this->fixedParams, otherHMM.fixedParams);
//  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst());*/ // Wed 12 Feb 2020 05:05:17 PM PST set in Optimizable deep copy ctor
//
//  this->numTransitionsMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
//  gsl_matrix_memcpy(this->numTransitionsMat, otherHMM.numTransitionsMat);
//}

// destructor
HMM::~HMM() {
  // TODO fix the DepthPair destructor. something to do with sharing data structures, probably
  for(std::vector<DepthPair*>::iterator it = this->depths->begin(); it != this->depths->end(); ++it) {
    delete *it;
  }
  delete this->alphabet;
  gsl_matrix_free(this->transition);
  gsl_matrix_free(this->rateMatrixQ);
  gsl_matrix_free(this->timeDepMatrixP);
  /*for(unsigned int i = 0; i < this->states->size(); i++) {
    for(unsigned int j = 0; j < this->states->size(); j++) {
      delete (*(*this->transitionStrings)[i])[j];
    }
    delete (*this->transitionStrings)[i];
  }
  delete this->transitionStrings;*/
  //for(unsigned int i = 0; i < this->states->size(); i++) {
  //  for(unsigned int j = 0; j < this->states->size(); j++) {
  for(unsigned int i = 0; i < this->rateMatrixQ->size1; i++) {
    for(unsigned int j = 0; j < this->rateMatrixQ->size2; j++) {
      delete (*(*this->rateMatrixStrings)[i])[j];
    }
    delete (*this->rateMatrixStrings)[i];
  }
  delete this->rateMatrixStrings;
  delete this->states;
  delete this->adjBinLabels;
  delete this->generator;

  for(std::vector<gsl_matrix*>::iterator it = this->forwardMatVec->begin(); it != this->forwardMatVec->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->forwardMatVec;

  for(std::vector<gsl_matrix*>::iterator it = this->backwardMatVec->begin(); it != this->backwardMatVec->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->backwardMatVec;

  for(std::vector<gsl_matrix*>::iterator it = this->forBackMatVec->begin(); it != this->forBackMatVec->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->forBackMatVec;

  for(std::vector<gsl_matrix*>::iterator it = this->backTraceVec->begin(); it != this->backTraceVec->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->backTraceVec;

  for(std::vector<gsl_vector*>::iterator it = this->scalingVecVec->begin(); it != this->scalingVecVec->end(); ++it) {
    gsl_vector_free(*it);
  }
  delete this->scalingVecVec;

  gsl_vector_free(this->initProb);
  gsl_vector_free(this->paramsToEst);
  gsl_vector_free(this->meanVarianceCoefVec);
  gsl_vector_free(this->bestOptimLLInitGuess);
  gsl_vector_free(this->prevForwardCol);
  gsl_vector_free(this->currForwardCol);
  gsl_vector_free(this->nextBackwardCol);
  gsl_vector_free(this->currBackwardCol);

  gsl_eigen_nonsymmv_free(this->ssEigenWorkspace);
  gsl_vector_complex_free(this->ssEval);
  gsl_matrix_complex_free(this->ssEvec);
  gsl_matrix_free(this->transitionTranspose);

  gsl_eigen_nonsymmv_free(this->rateEigenWorkspace);
  gsl_vector_complex_free(this->rateEval);
  gsl_matrix_complex_free(this->rateEvec);
  gsl_matrix_free(this->rateRealEvecMat);
  gsl_matrix_free(this->rateDiagMat);
  gsl_matrix_free(this->rateLUdecompMat);
  gsl_matrix_free(this->rateInverseMat);
  gsl_permutation_free(this->ratePerm);

  delete this->fixedParams;
  gsl_matrix_free(this->numTransitionsMat);
}

/*
 ********
 * accessors and mutators
 ********
 */
std::vector<DepthPair*>* HMM::getDepths() {
  return this->depths;
}
void HMM::setStates(std::vector<std::string>* states) {
  this->states = states;
}
std::vector<std::string>* HMM::getStates() const {
  return this->states;
}
void HMM::setAlphabet(std::set<int>* alphabet) {
  this->alphabet = alphabet;
}
std::set<int>* HMM::getAlphabet() const {
  return this->alphabet;
}
void HMM::setTransition(gsl_matrix* transition) {
  this->transition = transition;
}
double HMM::setTransition() {
  // subset out just transition probs from paramsToEst
  gsl_vector* tmpTrProbs = gsl_vector_alloc(this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  int tmpTrProbsIdx = 0;
  //printRowVector(this->paramsToEst);
  for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
    //std::cout << "from: " << TRANSITION_PROB_START_IDX + i  << " to: " << tmpTrProbsIdx << std::endl;
    gsl_vector_set(tmpTrProbs, tmpTrProbsIdx, gsl_vector_get(this->paramsToEst, TRANSITION_PROB_START_IDX + i));
    tmpTrProbsIdx++;
  }
  for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
    //std::cout << "from: " << BRANCH_LENGTH_START_IDX + i  << " to: " << tmpTrProbsIdx << std::endl;
    gsl_vector_set(tmpTrProbs, tmpTrProbsIdx, gsl_vector_get(this->paramsToEst, BRANCH_LENGTH_START_IDX + i));
    tmpTrProbsIdx++;
  }
  //printRowVector(tmpTrProbs);
  double status = this->setTransition(tmpTrProbs); // call subclass overridden version
  /*if(status != GSL_SUCCESS) {
    std::cerr << "HMM::setTransition caught nan status, returning nan" << std::endl;
  }*/
  gsl_vector_free(tmpTrProbs);
  return status;
}
// always saves into this->transition; delegates to subclass
double HMM::setTransition(gsl_vector* transitionParams) {
  double status = this->setTransition(this->transition, transitionParams);
  //std::cout << "HMM::setTransition this->setTransition(this->transition, transitionParams) STATUS: " << status << std::endl;
  if(status != GSL_SUCCESS) {
    return status;
  }
  return this->setInitProbSteadyState();
}
/*
 * Method to set rateMatrixQ in continuous time 2 adj loci model where
 * q_(i,j),(i',j') = mutation process forward in time of going from ploidy (i,j) to (i',j'), with rate lambda = P(event affects both bins)
 *              = { lambda * beta         if any same CNA > 1 (ie (i',j') \in {(i+n,j+n), or (i-n,j-n), n>1})
 *              = { lambda * (alpha+beta) if any same CNA = +-1 (ie (i',j') \in {(i+n,j+n), or (i-n,j-n), n=1})
 *              = { beta                  if any CNA > 1 affects only one cell (ie (i',j') \in {(i+n,j), (i-n,j), (i,j+n), (i,j-n), n>1})
 *              = { alpha+beta            if any CNA = +-1 affects only one cell (ie (i',j') \in {(i+n,j), (i-n,j), (i,j+n), (i,j-n), n=1})
 *              = { 0                     o/w (ie if n takes i or j < 0 or > k)
 *
 * rateParams = [alpha, beta, lambda]
 */
//void HMM::setRateMatrixQ(gsl_vector* rateParams) {
void HMM::setRateMatrixQ(double alpha, double beta, double lambda) {
  //double alpha = gsl_vector_get(rateParams, 0);
  //double beta = gsl_vector_get(rateParams, 1);
  //double lambda = gsl_vector_get(rateParams, 2);

  // if haven't saved the string represetation of the rate matrix yet, save it now
  bool saveStrs = false;
  std::string* currTrStr = nullptr;
  if(this->rateMatrixStrings == nullptr) {
    saveStrs = true;
    //this->rateMatrixStrings = new std::vector<std::vector<std::string*>*>(this->states->size());
    //for(unsigned int i = 0; i < this->states->size(); i++) {
    this->rateMatrixStrings = new std::vector<std::vector<std::string*>*>(this->rateMatrixQ->size1);
    for(unsigned int i = 0; i < this->rateMatrixQ->size1; i++) {
      //(*this->rateMatrixStrings)[i] = new std::vector<std::string*>(this->states->size());
      //for(unsigned int j = 0; j < this->states->size(); j++) {
      (*this->rateMatrixStrings)[i] = new std::vector<std::string*>(this->rateMatrixQ->size2);
      for(unsigned int j = 0; j < this->rateMatrixQ->size2; j++) {
        currTrStr = new std::string();
        *currTrStr += "(";
        (*(*this->rateMatrixStrings)[i])[j] = currTrStr;
      }
    }
  }

  // zero out everything
  gsl_matrix_set_zero(this->rateMatrixQ);

  int frI = 0;
  int frJ = 0;
  int toI = 0;
  int toJ = 0;
  double currRate = 0;
  //for(unsigned int row = 0; row < this->states->size(); row++) {
  for(unsigned int row = 0; row < this->rateMatrixQ->size1; row++) {
    frI = getIndvPloidyFromStateIdx(0, row);
    frJ = getIndvPloidyFromStateIdx(1, row);
    for(unsigned int col = 0; col < this->rateMatrixQ->size2; col++) {
      currTrStr = (*(*this->rateMatrixStrings)[row])[col];
      currRate = 0;
      toI = getIndvPloidyFromStateIdx(0, col);
      toJ = getIndvPloidyFromStateIdx(1, col);

      // no movement (main diagonal); will be set to 0-rowsum later
      if((frI == toI) && (frJ == toJ)){
        currRate = 0;
        if(saveStrs) {
          *currTrStr += "r)";
        }
        continue;
      }

      /*// Thu 21 Jan 2021 11:58:47 AM PST
      // if in a 0 state, can't get out of it (ie if you lose a segment, can't duplicate nothing into something)
      // ex 1,0 -> 1,1
      //    0,2 -> 1,2
      if((frJ == 0 && toJ != 0) || (frI == 0 && toI != 0)) { // TODO if we do this, then need to have some way of accounting for 0,0 -> 0,0 being an absorbing state, since that makes the matrix noninvertible (necessary to calculate the transition matrix)
        currRate = 0;
        if(saveStrs) {
          *currTrStr += "0)";
        }
      }*/
      // if move in tandem (include lambda term)
      // ex 0,1 -> 1,2  (+1)
      //    1,0 -> 2,1  (+1)
      //    2,3 -> 0,1  (-2)
      //    3,2 -> 1,0  (-2)
      //std::cout << row << "," << col << ": (" << frI << "," << frJ << ")->(" << toI << "," << toJ << ")\t" << frI-toI << ", " << frJ - toJ << ", " <<toI - frI << ", " <<toJ - frJ << std::endl;
      if((frI - toI) == (frJ - toJ)){// || (toI - frI) == (toJ - frJ)) {
        //std::cout << frI-toI << ", " << frJ - toJ << ", " <<toI - frI << ", " <<toJ - frJ << std::endl;
        // move only one in tandem (include alpha term)
        if(std::abs(frI - toI) == 1) {
          currRate += alpha;
          if(saveStrs) {
            *currTrStr += "a+";
          }
        }
        // add beta and lambda terms
        currRate += beta;
        currRate *= lambda;
        if(saveStrs) {
          *currTrStr += "b)*L";
        }
        //std::cout << "lambda term: (" << frI << "," << frJ << ")->(" << toI << "," << toJ << ")" << std::endl;
      }
      // else only one bin moves
      else if(frI == toI || frJ == toJ) {
        // move only one (include alpha term)
        if(std::abs(frI - toI) == 1 || std::abs(frJ - toJ) == 1) {
          currRate += alpha;
          if(saveStrs) {
            *currTrStr += "a+";
          }
      //std::cout << row << "," << col << ": (" << frI << "," << frJ << ")->(" << toI << "," << toJ << ")\t" << frI-toI << ", " << frJ - toJ << ", " <<toI - frI << ", " <<toJ - frJ << std::endl;
        }
        // add beta term
        currRate += beta;
        if(saveStrs) {
          *currTrStr += "b)";
        }
      }
      // else, some impossible transition
      else {
        currRate = 0;
        if(saveStrs) {
          *currTrStr += "0)";
        }
      }
      gsl_matrix_set(this->rateMatrixQ, row, col, currRate);
    }
  } 

  // set diagonals to be -rowsum so entire row sums to 0 (see https://en.wikipedia.org/wiki/Continuous-time_Markov_chain#Definition)
  for(unsigned int row = 0; row < this->rateMatrixQ->size1; row++) {
    double rowSum = 0;
    for(unsigned int col = 0; col < this->rateMatrixQ->size1; col++) {
      rowSum += gsl_matrix_get(this->rateMatrixQ, row, col);
    }
    gsl_matrix_set(this->rateMatrixQ, row, row, 0 - rowSum);
  }
}

/*
 * Function to set time dependent transition matrix P = exp(Q*t)
 * where t=time. Assumes rateMatrixQ and all intermediates are set up already.
 * returns status from matrixExponential
 */
double HMM::setTimeDepMatrixP(double time) {
  return this->setTimeDepMatrixP(this->timeDepMatrixP, time);
}
double HMM::setTimeDepMatrixP(gsl_matrix* destMat, double time) {
  double status = matrixExponential(destMat, this->rateMatrixQ, time, this->rateEigenWorkspace, this->rateEval, this->rateEvec, this->rateRealEvecMat, this->rateDiagMat, this->rateLUdecompMat, this->ratePerm, this->rateInverseMat);

  if(status != GSL_SUCCESS) {
    //std::cerr << "HMM::setTimeDepMatrixP matrixExponential failed, returning nan" << std::endl;
    return status;
  }

  // ensure each row sums to 1 by rescaling by rowsum and by getting rid of numerical error negatives
  double scalingFactor = 0;
  double currVal = 0;
  for(unsigned int row = 0; row < destMat->size1; row++) {
    scalingFactor = 0;
    for(unsigned int col = 0; col < destMat->size2; col++) {
      currVal = gsl_matrix_get(destMat, row, col);
      // if an entry is very negative, something went wrong
      if(currVal < -1e-13) {
        return GSL_NAN;
      }
      // if small numerical error, set to 0
      else if(currVal < 0) {
        gsl_matrix_set(destMat, row, col, 0);
        currVal = 0;
      }
      scalingFactor += currVal;
    }
    gsl_vector_view currRow = gsl_matrix_row(destMat, row);
    gsl_vector_scale(&currRow.vector, 1.0 / scalingFactor);
  }
  return status;
}

gsl_matrix* HMM::getTransition() const {
  return this->transition;
}
gsl_matrix* HMM::getBaumWelchTransitionMat() const {
  return this->baumWelchTransitionMat;
}
gsl_matrix* HMM::getNumTransitionsMat() const {
  return this->numTransitionsMat;
}
int HMM::getKploidy() const {
  //return (int) this->transition->size2 - 1; // matrix goes from 0 to k-ploid
  return this->MAX_PLOIDY;
}
void HMM::setInitProb(gsl_vector* initProb) {
  //this->initProb = initProb;
  gsl_vector_memcpy(this->initProb, initProb);
}
gsl_vector* HMM::getInitProb() const {
  return this->initProb;
}
void HMM::setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  if(this->meanVarianceCoefVec != nullptr) {
    gsl_vector_free(this->meanVarianceCoefVec);
  }
  this->meanVarianceCoefVec = meanVarianceCoefVec;
}
gsl_vector* HMM::createMeanVarianceCoefVec() {
  // Wed 25 Apr 2018 11:44:01 AM PDT
  // according to Navin_Nature2011/fitLikelihood.R with "diploidBins.q20.adaptTr.qualTr.dusted.sorted.noDups.overallMinLibFALSE.normCov_250kb/diploidBinsVarMeanNormLibSize_250kb.txt", without lm outliers, requiring at least 5 reads (fitNoLmOutliers_poly2_5)
  // variance = 53.62011293244999166063 + mean * 0.62350146443559151255 + mean^2 * 0.01128075188930972687
  /*double intercept = 53.62011293244999166063;
  double slope = 0.62350146443559151255;
  double poly2 = 0.01128075188930972687;*/
  double intercept = 10.46385711652957084539; // Fri 17 Apr 2020 11:02:53 PM PDT /space/s1/sandra/alleleFreqHmm/Navin_Nature2011/fitMeanVarRlnshp_noDownsampling.R
  double slope = 2.42601321762369614987;
  double poly2 = 0.01114518215725581601;

  gsl_vector* meanVarianceCoefVec = gsl_vector_alloc(3);
  gsl_vector_set(meanVarianceCoefVec, 0, intercept);
  gsl_vector_set(meanVarianceCoefVec, 1, slope);
  gsl_vector_set(meanVarianceCoefVec, 2, poly2);

  return meanVarianceCoefVec;
}
double HMM::getMeanVarianceIntercept() const {
  return gsl_vector_get(this->meanVarianceCoefVec, 0);
}
double HMM::getMeanVarianceSlope() const {
  return gsl_vector_get(this->meanVarianceCoefVec, 1);
}
double HMM::getMeanVariancePoly2() const {
  return gsl_vector_get(this->meanVarianceCoefVec, 2);
}
gsl_vector* HMM::getMeanVarianceCoefVec() const {
  return this->meanVarianceCoefVec;
}
//gsl_vector* HMM::getBestImpInitGuess() const {
//  return this->bestOptimLLInitGuess;
//}
double HMM::setParamsToEst(gsl_vector* params) {
  /*if(this->paramsToEst != nullptr && this->paramsToEst != params) {
    gsl_vector_free(this->paramsToEst);
  }
  this->paramsToEst = params;*/
  gsl_vector_memcpy(this->paramsToEst, params); // Thu 30 Jan 2020 10:47:35 AM PST changed to be a memcpy instead of free/alloc cycle
  //std::cout << "IN HMM SETPARAMSTOEST: ";
  //printRowVector(this->paramsToEst);
  return this->setTransition();
}

double HMM::setFixedParams(gsl_vector* fixedParams) {
  //gsl_vector_free(this->fixedParams);
  //this->fixedParams = fixedParams;
  gsl_vector_memcpy(this->fixedParams, fixedParams); // Thu 30 Jan 2020 11:14:49 AM PST changed to be a memcpy

  return this->setTransition();
}
void HMM::setSimParamsToEst(gsl_vector* params) {
  this->simParamsToEst = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simParamsToEst, params);
}
void HMM::setSimFixedParams(gsl_vector* params) {
  this->simFixedParams = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simFixedParams, params);
}

// list of all chromosomes from first DepthPair
std::vector<std::string>* HMM::getChrVec() const {
  // since all DepthPairs share the same vector of chromosomes (assuming second and subsequent DepthPairs are constructed using the first DepthPair as a reference), just pull out the chrVec from the first DepthPair from this->depths vector
  return (*this->depths)[0]->chrVec;
}

/*
 * function to set this's initial probability vector to the steady state distribution
 */
double HMM::setInitProbSteadyState() {
  // if initProb doesn't exist, allocate
  if(this->initProb == nullptr) {
    this->initProb = gsl_vector_alloc(this->states->size());
  }
  // if initProb was set to (wrong) default size earlier, resize
  if(this->initProb->size != this->states->size()) {
    gsl_vector_free(this->initProb);
    this->initProb = gsl_vector_alloc(this->states->size());
  }
  return this->findSteadyStateDist(this->initProb);
}

/*
 * helper function to set the library scaling factors to the total
 * read depth in the tumor cell divided by the total average diploid read depth.
 * This function is only to be called outside of the ctor since
 * it relies on polymorphism for the correct setLibScalingFactor() fn to be called
 */
void HMM::setLibScalingFactorsToTotalRatio() {
  //for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
  for(int i = 0; i < this->NUM_CELLS; i++) { // TODO debugging libRatio
    //this->setLibScalingFactor(i, totalTumorDepth / totalAvgDiploidDepth);
    this->setLibScalingFactor(i, this->calcLibScalingFactorsToTotalRatio(i));
  }
}
double HMM::calcLibScalingFactorsToTotalRatio(int cellIdx) const {
  double totalAvgDiploidDepth = (*this->depths)[0]->getTotalDiploidDepth();
  double totalTumorDepth = (*this->depths)[cellIdx]->getTotalTumorDepth();;
  return totalTumorDepth / totalAvgDiploidDepth;
}
void HMM::setLibScalingFactors(gsl_vector* libScalingFactors) {
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    //gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, gsl_vector_get(libScalingFactors, i));
    this->setLibScalingFactor(i, gsl_vector_get(libScalingFactors, i));
  }
}

/*
 * function to help estimate library scaling factors based on the posterior decoding.
 * Specifically, let
 * - Xij = reads in tumor j in window i
 * - prev Xij ~ NB (lambda_i = lib_j * ploidy_i * diploid_i / 2 + epsilon = dipoid_i / 50)
 * - now, Xij ~ NB (r_ik + epsilon)
 *   where r_ik = k * diploid_i
 *                --------------------------------------------------------- * T_j
 *                n*epsilon + sum_(j=1)^n diploid_j * [most likely ploidy in window j]
 *         k = current ploidy
 *         i = window num
 *         T_j = total num reads in tumor j
 *         //E[ploidy in window j] = expected ploidy in window j from posterior decoding (ie mult probability of each ploidy by ploidy in window j (ie weighted avg))
 *         [most likely ploidy in window j] = most likely ploidy in window j from viterbi decoding
 *         n = # windows
 * - sanity check: from sim, know prev lib scaling. should = 2 * T_j / (denom of r_ik)
 *
 * This function only calculate the denominator (times T_j) of the above, since the numerator
 * is set during calls to getEmissionProb
 */
void HMM::miscFunctions() {
  //std::cerr << "IN HMM miscFunctions()" << std::endl;
  return;
  this->estLibScalingFactorsPosterior();
}
void HMM::estLibScalingFactorsPosterior() {
  //this->setLibScalingFactorsToTotalRatio();
  //std::cerr << "after HMM::setLib" << std::endl;
  //return; // Tue 17 Mar 2020 06:13:46 PM PDT debugging: is the lib size messing things up? try setting lib sizes to true scaling factors from simulation
  this->viterbiDecode(); // only need to call once per pair; needed for viterbi version
  //this->calcMargLikelihoods(); // needed for for/back version
  for(int i = 0; i < this->NUM_CELLS; i++) {
    this->setLibScalingFactor(i, this->estLibScalingFactorsPosterior(i));
  }

}
double HMM::estLibScalingFactorsPosterior(int cellNum) {
  DepthPair* currDepths = (*this->depths)[cellNum];

  // viterbi version
  std::vector<double>* diploidDepthVec = nullptr;
  std::string currChr;
  double sum = 0;
  double expPloidy = 0;
  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < currDepths->chrVec->size(); chrIdx++) {
    currChr = (*currDepths->chrVec)[chrIdx];
    diploidDepthVec = (*currDepths->chrToDiploidDepthMap)[currChr];

    // for each window in chr
    for(unsigned int windowIdx = 0; windowIdx < (*currDepths->regions)[currChr]->size(); windowIdx++) {
      expPloidy = (*(*currDepths->chrToViterbiPathMap)[currChr])[windowIdx];
      //std::cout << "expPloidy: " << expPloidy << ",\tdiploid reads: " << (*diploidDepthVec)[windowIdx] << std::endl;
      //sum += (*diploidDepthVec)[windowIdx] * expPloidy;
      sum += (*diploidDepthVec)[windowIdx]/2.0 * expPloidy;
    }
  }

  /*// avg viterbi version; same as viterbi += 3 decimal places
  std::string currChr;
  double avgPloidy = 0;
  double totPloidy = 0;
  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < currDepths->chrVec->size(); chrIdx++) {
    currChr = (*currDepths->chrVec)[chrIdx];
    //diploidDepthVec = (*currDepths->chrToDiploidDepthMap)[currChr];

    // for each window in chr
    for(unsigned int windowIdx = 0; windowIdx < (*currDepths->regions)[currChr]->size(); windowIdx++) {
      totPloidy += (*(*currDepths->chrToViterbiPathMap)[currChr])[windowIdx];
      //std::cout << "totPloidy: " << totPloidy << ",\tnumWindows " << currDepths->numWindows << std::endl;
      //sum += (*diploidDepthVec)[windowIdx] * expPloidy;
      //sum += (*diploidDepthVec)[windowIdx]/2.0 * expPloidy;
    }
  }
  avgPloidy = totPloidy / currDepths->numWindows;*/

  /*// forward/backward version // Wed 06 May 2020 03:40:51 PM PDT exceedingly slow, not sure why, but same as viterbi += 3 decimal places
  std::vector<double>* diploidDepthVec = nullptr;
  gsl_matrix* currForBackMat = nullptr;
  std::string currChr;
  double sum = 0;
  double expPloidy = 0;
  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < currDepths->chrVec->size(); chrIdx++) {
    currChr = (*currDepths->chrVec)[chrIdx];
    diploidDepthVec = (*currDepths->chrToDiploidDepthMap)[currChr];
    currForBackMat = (*currDepths->forBackMargMatMap)[currChr];

    for(int state = 0; state < this->MAX_PLOIDY+1; state++) {
      // for each window in chr
      for(unsigned int windowIdx = 1; windowIdx < (*currDepths->regions)[currChr]->size(); windowIdx++) {
       // expPloidy = 0;
        expPloidy = state * exp(gsl_matrix_get(currForBackMat, state, windowIdx));
        //std::cout << "expPloidy: " << expPloidy << ",\t" << state << ", " << exp(gsl_matrix_get((*currDepths->forBackMargMatMap)[currChr], state, windowIdx)) << std::endl;
        sum += (*diploidDepthVec)[windowIdx]/2.0 * expPloidy;
      }
    }
  }*/

  /*// stat dist
  // let X_ij ~ NB(ploidy_i * DR_i / 2 + epsilon) * c_j
  // where c_j = T_j / (n * epsilon + sum_windows (DR_j / 2 * weighted avg ploidy (by stationary dist) + epsilon))
  // get weighted avg ploidy
  int currPloidy = -1;
  double weightedAvgPloidy = 0;
  for(unsigned int stateIdx = 0; stateIdx < this->initProb->size; stateIdx++) {
    currPloidy = getCellPloidyFromStateIdx(cellNum, stateIdx);
    weightedAvgPloidy += currPloidy * gsl_vector_get(this->initProb, stateIdx);
  }*/

  // finally, multiply everything together
  //double nEps = currDepths->numWindows * (currDepths->diploidLibrarySize / this->DEPTH_ERROR_SCALING);
  //double nEps = (currDepths->diploidLibrarySize / this->DEPTH_ERROR_SCALING);
  //double nEps = currDepths->numWindows / this->DEPTH_ERROR_SCALING;
  //double nEps = currDepths->diploidLibrarySize / this->DEPTH_ERROR_SCALING;
  //double nEps = currDepths->numWindows * 65.37034 / this->DEPTH_ERROR_SCALING;
  double nEps = currDepths->numWindows * 272.5568 / this->DEPTH_ERROR_SCALING; // unnormDiploid error term
  //double nEps = 0; // Mon 04 May 2020 06:32:51 PM PDT noErr debugging
  double T = currDepths->tumorLibrarySize;

  //double scalingFactor = T / (nEps + sum); // viterbi version
  double scalingFactor = (T - nEps) / sum; // viterbi version with const err Tue 08 Jun 2021 05:26:51 PM PDT
  //double scalingFactor = T / (nEps + avgPloidy * currDepths->diploidLibrarySize / 2.0); // avg viterbi version
  //double scalingFactor = T / (currDepths->diploidLibrarySize / 2.0 * weightedAvgPloidy + nEps); // stat dist

  // adr  cellNum  nEps  T  sum  2T/sum  T/(nEps+sum)
  //fprintf(stderr, "%%%p\t%i\t%.2f\t%.2f\t%.2f\t%.10f\t%.10f\n", this, cellNum, nEps, T, sum, 2.0*T/sum, T/(nEps+sum));
  if(this->gradientDebug) {
    fprintf(stderr, "%%%p\t%i\t%.2f\t%.2f\t%.2f\t%.10f\t%.10f\t%.10f\t%.10f\n", this, cellNum, nEps, T, sum, (T-nEps)/sum, T/(nEps+sum), currDepths->diploidLibrarySize, currDepths->tumorLibrarySize / currDepths->diploidLibrarySize); // viterbi version
    //fprintf(stderr, "%%%p\t%i\t%.2f\t%.2f\t%.2f\t%.10f\n", this, cellNum, nEps, T, avgPloidy, T / (nEps + avgPloidy * currDepths->diploidLibrarySize / 2.0));
    //fprintf(stderr, "%%%p\t%i\t%.2f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", this, cellNum, nEps, weightedAvgPloidy, scalingFactor, T, currDepths->diploidLibrarySize, currDepths->diploidLibrarySize / 2.0, (currDepths->diploidLibrarySize / 2.0 * weightedAvgPloidy + nEps)); // stationary dist version
  }

  //std::cout << "libScale: " << 2.0 * T / sum << "\t";//<< ", T: " << T << ", sum: " << sum << ", nEps: " << nEps << std::endl;

  //return T / (nEps + sum); // viterbi version
  return scalingFactor;
  //return 2.0 *(T - nEps) / sum;
  //return 2.0 * T / (nEps + sum);
}

/*
 * helper function to convert cellIdx and stateIdx to a ploidy, instead of storing and looking up
 * see GetCellPloidyFromStateIdx.java for testing/development code
 * ex. getCellPloidyFromStateIdx(1, 4) if MAX_PLOIDY==3 and NUM_CELLS==2
 *     states = [(0,0), (0,1), (0,2), (0,3), (1,0), ...]
 *     returns (1,0) ==> 0
 */
int HMM::getCellPloidyFromStateIdx(int cellIdx, int stateIdx) const {
  //return stateIdx / ((int) std::pow(this->MAX_PLOIDY + 1, this->NUM_CELLS - (cellIdx + 1))) % (this->MAX_PLOIDY + 1); // for dummy data, std::pow takes 21 sec, looping takes 21 sec, lookup in parsedStates takes 31 sec
  int div = 1; // unclear if std::pow or loop is faster on small dummy data TODO
  for(int i = 0; i < this->NUM_CELLS - (cellIdx + 1); i++) {
    div *= this->MAX_PLOIDY + 1;
  }
  return stateIdx / div % (this->MAX_PLOIDY + 1);
}

/*
 * helper function to convert a stateIdx that corresponds only to a pair of ploidies into the idxInPair's ploidy.
 * assumes there are only 2 ploidies in this stateIdx, and uses code from getCellPloidyFromStateIdx with 2 hard coded.
 * useful for converting a pair of bins into an individual ploidy when calculting rateMatrixQ
 */
int HMM::getIndvPloidyFromStateIdx(int idxInPair, int stateIdx) const {
  int div = 1;
  for(int i = 0; i < 2 - (idxInPair + 1); i++) {
    div *= this->MAX_PLOIDY + 1;
  }
  return stateIdx / div % (this->MAX_PLOIDY + 1);
}

/*
 * helper function to convert a pair of ploidies into a state idx.
 * assumes there are only 2 ploidies (ie must be a pair).
 * useful for getting the stateIdx or timeDepMatrixP index from a pair of ploidies
 */
int HMM::getStateIdxFromPloidyPair(int ploidy0, int ploidy1) const {
  return ploidy0 * (this->MAX_PLOIDY + 1) + ploidy1;
}

double HMM::getLogLikelihood() {
  return this->runForwardAlg();
}

/*
 * returns the logliklihood of the path taken by the viterbi decoding,
 * rather than the full forward loglikelihood
 */
double HMM::getViterbiLogLikelihood() {
  std::vector<std::string>* chrVec = this->getChrVec();
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;

  // first run viterbi decoding
  this->viterbiDecode();

  double vitLl = 0;
  double diploidDepth = 0;
  double tumorDepth = 0;
  int prevVitState = 0;
  int currVitState = 0;
  double emissionProb = 0;
  double transitionProb = 0;
  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    std::string currChr = (*chrVec)[chrIdx];

    for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      diploidDepth = (*(*currDepths->chrToDiploidDepthMap)[currChr])[0];
      tumorDepth = (*(*currDepths->chrToTumorDepthMap)[currChr])[0];
      prevVitState = -1;
      currVitState = (*(*currDepths->chrToViterbiPathMap)[currChr])[0];

      emissionProb = this->getEmissionProb(tumorDepth, diploidDepth, currVitState, cellIdx);
      transitionProb = gsl_vector_get(this->initProb, currVitState); //gsl_matrix_get(this->transition, prevVitState, currVitState); // won't work for anything other than one cell

      //fprintf(stderr, "%.5f\t%.5f\t%i\t%i\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n", tumorDepth, diploidDepth, prevVitState, currVitState, emissionProb, transitionProb, log(emissionProb), 0.0, log(emissionProb) + 0); // NOTE: checkPerWindowLogLik uncomment HERE
      vitLl += log(emissionProb) + 0;
    }

    // for each window in currChr
    for(unsigned int regionIdx = 1; regionIdx < (*firstDepthPair->regions)[currChr]->size(); regionIdx++) {

      for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currDepths = (*this->depths)[cellIdx];
        diploidDepth = (*(*currDepths->chrToDiploidDepthMap)[currChr])[regionIdx];
        tumorDepth = (*(*currDepths->chrToTumorDepthMap)[currChr])[regionIdx];
        prevVitState = (*(*currDepths->chrToViterbiPathMap)[currChr])[regionIdx - 1];
        currVitState = (*(*currDepths->chrToViterbiPathMap)[currChr])[regionIdx];

        emissionProb = this->getEmissionProb(tumorDepth, diploidDepth, currVitState, cellIdx);
        transitionProb = gsl_matrix_get(this->transition, prevVitState, currVitState); // won't work for anything other than one cell

        //fprintf(stderr, "%.5f\t%.5f\t%i\t%i\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n", tumorDepth, diploidDepth, prevVitState, currVitState, emissionProb, transitionProb, log(emissionProb), log(transitionProb), log(emissionProb) + log(transitionProb)); // NOTE: checkPerWindowLogLik uncomment HERE
        vitLl += log(emissionProb) + log(transitionProb);
      }
    }
  }
  return vitLl;
}

/*
 * returns forward likelihood
 */
double HMM::runForwardAlg() {
  std::string currChr;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  gsl_vector* scalingVec = nullptr;
  gsl_matrix* forwardMat = nullptr;

  //int currPloidy = -1;
  int stateIdx = -1;
  int cellIdx = -1;
  //int windowIdx = -1;
  double scalingFactor = 0;
  //double currTumorDepth = -1;
  //double currDiploidDepth = -1;
  //double currRes = 0;
  //double initProb = 0;
  double emissionProb = -1;

  double currResPloidy = 0;

  double sumLogScalingFactors = 0;
  double sumLastCol = 0;
  double logScaledSumLastCol = 0;
  double totalLogLikelihood = 0;

  unsigned int i = 0;
  unsigned int j = 0;
  int depthIdx = 0;
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    scalingVec = (*this->scalingVecVec)[chrIdx];
    forwardMat = (*this->forwardMatVec)[chrIdx];
    gsl_vector_set_all(scalingVec, 1.0);
    gsl_matrix_set_all(forwardMat, 0);
    scalingFactor = 0;

    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    //std::cerr << "before first col" << std::endl;
    // set first col to initProb * emission prob for first sym
   /* currTumorDepth = -1;
    //currDiploidDepth = -1;
    //currRes = 0;
    //initProb = 0;
    emissionProb = 1;
    currPloidy = -1;
    currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[0];*/
    //windowIdx++;
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      /*//emissionProb = 1;
      emissionProb = 0;

      // for each cell, calc the emission prob and multiply by it
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
        currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[0];
        //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
        //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx); // USE THIS ONE; Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
        //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx)); // USE THIS ONE; Wed 24 Jun 2020 03:43:13 PM PDT precision debugging
        emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, 0, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT

        //// Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
        //int simState = std::atoi((*(*(*this->depths)[cellIdx]->chrToTumorSimStateMap)[currChr])[0].c_str()); // get sim state
        //if(simState == currPloidy) {
        //  emissionProb *= 1;
        //} else {
        //  emissionProb *= 0;
        //}



      //std::cout << "emmissionProb: " << emissionProb << ", " << currTumorDepth << std::endl;
      }*/
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);
      //std::cout << "stateIdx " << stateIdx << ", chrIdx: " << chrIdx << ", depth " << 0 << ", emission: " << emissionProb << std::endl;
      //std::cout << "tot calc: " <<  << std::endl;
      //std::cout << "mult emmissionProb: " << emissionProb << std::endl;
      //initProb = gsl_vector_get(this->initProb, stateIdx);
      //currResPloidy = emissionProb * initProb;
      //currResPloidy = emissionProb * gsl_vector_get(this->initProb, stateIdx);
      //gsl_matrix_set(forwardMat, stateIdx, 0, currResPloidy);
      currResPloidy = emissionProb + log(gsl_vector_get(this->initProb, stateIdx));

      gsl_matrix_set(forwardMat, stateIdx, 0, exp(currResPloidy));
    }
    //gsl_matrix_set_col(forwardMat, 0, this->initProb);

    // rescale first col by max
    gsl_vector_view firstCol = gsl_matrix_column(forwardMat, 0);
    //printRowVector(stderr, &firstCol.vector);
    scalingFactor = gsl_vector_max(&firstCol.vector);
    //scalingFactor = 1.0;
    gsl_vector_scale(&firstCol.vector, 1.0 / scalingFactor);
    gsl_vector_set(scalingVec, 0, scalingFactor);
    //printMatrix(forwardMat);

    /*for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      currResPloidy = gsl_matrix_get(forwardMat, stateIdx, 0);
      if(compareDoubles(currResPloidy, 0, 1e-200)) {
        emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);
        std::cout << "CURRRESPLOIDY(" << stateIdx << ") IS 0: " << currResPloidy << ", " << emissionProb << ", " << gsl_vector_get(this->initProb, stateIdx) << ", " << log(gsl_vector_get(this->initProb, stateIdx)) << std::endl;
      }
    }*/

    //std::cout << "after first col" << std::endl;
    // iter through observed sequence
    gsl_vector_set_zero(this->prevForwardCol);
    gsl_vector_set_zero(this->currForwardCol);
    currResPloidy = 0;
  //for(int tmp = 0; tmp < this->diploidDepthPloidyPreCalc->size(); tmp++) {
  //  std::cout << (*this->diploidDepthPloidyPreCalc)[tmp] << ", ";
  //}
  //std::cout << std::endl;
    //for(i = 1; i < (*firstDepthPair->regions)[currChr]->size(); i++) {
    //for(i = 1, depthIdx = 0; i < forwardMat->size2; i++, depthIdx++) {
    for(i = 1, depthIdx = 1; i < forwardMat->size2; i++, depthIdx++) {
      gsl_vector_set_zero(this->currForwardCol);
      gsl_matrix_get_col(this->prevForwardCol, forwardMat, i-1);

      //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[i];
      //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
      //std::cout << "considering depthidx: " << depthIdx << ", forwardMat idx i: " << i << std::endl;
      //windowIdx++;
      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {
        ////emissionProb = 1;
        //emissionProb = 0;

        ////std::chrono::steady_clock::time_point tmpBegin = std::chrono::steady_clock::now();
        //// for each cell, calc the emission prob and multiply by it
        //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        //  currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
        //  //currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[i];
        //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
        //  //std::cout << currTumorDepth << ", " << currDiploidDepth << ", " << currPloidy << ", " << cellIdx << std::endl;
        //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
        //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, windowIdx, cellIdx);
        //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx); // Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
        //  //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx)); // Wed 24 Jun 2020 03:43:13 PM PDT precision debugging
        //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, depthIdx, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT

        //  /*// Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
        //  int simState = std::atoi((*(*(*this->depths)[cellIdx]->chrToTumorSimStateMap)[currChr])[i].c_str()); // get sim state
        //  if(simState == currPloidy) {
        //    emissionProb *= 1;
        //  } else {
        //    emissionProb *= 0;
        //  }*/






        ////std::cout << "emmissionProb: " << emissionProb << ", " <<  currTumorDepth << std::endl;
        //}
        //std::cout << emissionProb - this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx) << std::endl;
        emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);
        //std::cout << "stateIdx " << stateIdx << ", chrIdx: " << chrIdx << ", depth " << depthIdx << ", emission: " << emissionProb << std::endl;
        //std::chrono::steady_clock::time_point tmpEnd = std::chrono::steady_clock::now();
        //double tmpElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(tmpEnd - tmpBegin).count() / 1000000.0;
        //std::cout << tmpElapsedSec << std::endl;
        //fprintf(stdout, "ONE emi\t%.80f\n", tmpElapsedSec);
        //std::cout << "mult emmissionProb: " << emissionProb << std::endl;
        // iter over possible transitions
        for(j = 0; j < this->prevForwardCol->size; j++) {
          //this->currForwardCol[stateIdx] += emissionProb * this->prevForwardCol[j] * this->transition[j][stateIdx];
          currResPloidy = gsl_vector_get(this->currForwardCol, stateIdx);
          //gsl_vector_set(this->currForwardCol, stateIdx, currResPloidy + emissionProb * gsl_vector_get(this->prevForwardCol, j) * gsl_matrix_get(this->transition, j, stateIdx));
          gsl_vector_set(this->currForwardCol, stateIdx, currResPloidy + exp(emissionProb + log(gsl_vector_get(this->prevForwardCol, j)) + log(gsl_matrix_get(this->transition, j, stateIdx))));
        }
      }
      //std::cout << std::endl;

      // find max of this->currForwardCol to scale with
      scalingFactor = gsl_vector_max(this->currForwardCol);
      //scalingFactor = 1.0;

      // scale this->currForwardCol by scalingFactor
      gsl_vector_scale(this->currForwardCol, 1.0 / scalingFactor);

      // save this->currForwardCol and scalingFactor
      //std::cout << "scalingVec[" << i << "]: " << gsl_vector_get(scalingVec, i) << ", size: " << scalingVec->size << std::endl;
      gsl_vector_set(scalingVec, i, scalingFactor);
      gsl_matrix_set_col(forwardMat, i, this->currForwardCol);
      //printMatrix(forwardMat);
    //printRowVector(stdout, this->currForwardCol);
    }

    // return unscaled forward prob (prob of observed seq given model params)
    sumLogScalingFactors = 0;
    for(i = 0; i < scalingVec->size; i++) {
      //std::cout << "forward scalingVec[" << i << "]: " << gsl_vector_get(scalingVec, i) << ", log(): " << log(gsl_vector_get(scalingVec, i)) << std::endl;
      //fprintf(stdout, "%.20f\n", log(gsl_vector_get(scalingVec, i))); // NOTE: checkPerWindowLogLik uncomment HERE
      sumLogScalingFactors = sumLogScalingFactors + log(gsl_vector_get(scalingVec, i));
    }

    // General case: P(x) = (sum(unscaled last col) / last scaling factor) * product(all scaling factors)
    //                    = sum(scaled last col) * product(all scaling factors)
    //     <==> log(P(x)) = log(sum(last col) / last scaling factor) + sum(log(scaling factors))
    //                    = log(sum(scaled last col)) + sum(log(scaling factors))
    // if scaling factor = sum(col), then
    //   log(sum(last col) / last scaling factor) = log(sum(last col) / sum(last col)) = 0, so
    //   log(P(x)) = sum(log(scaling factors)) <==> P(x) = product(scaling factors)
    // if scaling factor = max(col), then
    //   log(P(x)) = log(sum(last col) / max(last col)) + sum_{all cols}(log(max(col)))
    //             = log(sum(scaled last cols)) + sum_{all cols}(log(max(col)))
    gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
    sumLastCol = gsl_blas_dasum(&lastCol.vector); // Double Absolute SUM
    //printColVector((gsl_vector*)&lastCol.vector);
    logScaledSumLastCol = log(sumLastCol);
    totalLogLikelihood += sumLogScalingFactors + logScaledSumLastCol;
    //std::cout << currChr << ", " << sumLogScalingFactors + logScaledSumLastCol << std::endl;
    //std::cout << "logScaledSumLastCol:" << logScaledSumLastCol << ", sumLogScalingFactors: " << sumLogScalingFactors << ", totalLogLikelihood: " << totalLogLikelihood << ", likelihood: " << exp(totalLogLikelihood) << std::endl;
  }

  if(isnan(totalLogLikelihood) || isnan(-totalLogLikelihood)) {
    totalLogLikelihood = GSL_NAN;
  }

  //delete currChrDepthsVec;
  return totalLogLikelihood;
}

/*
 * returns backward likelihood. Assumes forward algorithm has been run first (uses same scaling constants)
 * see https://web.stanford.edu/~jurafsky/slp3/A.pdf and https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm#Backward_probabilities for references
 */
double HMM::runBackwardAlg() {
  std::string currChr;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  gsl_vector* scalingVec = nullptr;
  gsl_matrix* backwardMat = nullptr;

  //int currPloidy = -1;
  int stateIdx = -1;
  int cellIdx = -1;
  double scalingFactor = 0;
  //double currTumorDepth = -1;
  //double currDiploidDepth = -1;
  double currRes = 0;
  double initProb = 0;
  double emissionProb = -1;
  double currResPloidy = 0;

  double sumLogScalingFactors = 0;
  double sumLastCol = 0;
  double logScaledSumLastCol = 0;
  double totalLogLikelihood = 0;

  int i = 0;
  int j = 0;
  int depthIdx = 0;
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    scalingVec = (*this->scalingVecVec)[chrIdx];
    backwardMat = (*this->backwardMatVec)[chrIdx];
    gsl_matrix_set_all(backwardMat, 0);

    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    //std::cerr << "before first col" << std::endl;
    // set last col to 1
    //currTumorDepth = -1;
    currRes = 0;
    emissionProb = 1;
    //currPloidy = -1;
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      gsl_matrix_set(backwardMat, stateIdx, backwardMat->size2 - 1, 1.0); // last col is just 1's
    }

    // iter through observed sequence
    gsl_vector_set_zero(this->nextBackwardCol); // next in sequence, ie to the right/closer to the end
    gsl_vector_set_zero(this->currBackwardCol);
    currResPloidy = 0;
    for(i = backwardMat->size2 - 2, depthIdx = (*firstDepthPair->regions)[currChr]->size() - 1; i >= 0 ; i--, depthIdx--) { // -1 would get last col
    //for(i = backwardMat->size2 - 1, depthIdx = (*firstDepthPair->regions)[currChr]->size() - 1; i >= 0 ; i--, depthIdx--) { // -1 would get last col
      gsl_vector_set_zero(this->currBackwardCol);
      gsl_matrix_get_col(this->nextBackwardCol, backwardMat, i+1);
      //std::cout << "nextBackwardCol idx: " << i+1 << std::endl;

      //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
      //std::cout << "considering depthidx: " << depthIdx << ", backwardMat idx i: " << i << std::endl;
      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {

        // iter over possible transitions
        // this->currBackwardCol[stateIdx] += transition[stateIdx][j] * this->nextBackwardCol[j] * emissionProb(of j'th symbol)
        for(j = 0; j < (int) this->nextBackwardCol->size; j++) {
          currResPloidy = gsl_vector_get(this->currBackwardCol, stateIdx);

          //// for each cell, calc the emission prob and multiply by it
          ////emissionProb = 1;
          //emissionProb = 0;
          //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          //  currPloidy = getCellPloidyFromStateIdx(cellIdx, j);
          //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
          //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, windowIdx, cellIdx);
          //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx); // USE THIS ONE; Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
          //  //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx)); // USE THIS ONE; Wed 24 Jun 2020 03:43:13 PM PDT precision debugging
          //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, depthIdx, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT
          //  //std::cout << cellIdx << ", " << currTumorDepth << ", " << this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, depthIdx, cellIdx) << std::endl;
          //}
          //std::cout << emissionProb - this->getTotalLogEmissionProb(j, currChrDepthsVec, chrIdx, depthIdx) << std::endl;
          emissionProb = this->getTotalLogEmissionProb(j, currChrDepthsVec, chrIdx, depthIdx);
          //gsl_vector_set(this->currBackwardCol, stateIdx, currResPloidy + emissionProb * gsl_vector_get(this->nextBackwardCol, j) * gsl_matrix_get(this->transition, stateIdx, j));
          gsl_vector_set(this->currBackwardCol, stateIdx, currResPloidy + exp(emissionProb + log(gsl_vector_get(this->nextBackwardCol, j)) + log(gsl_matrix_get(this->transition, stateIdx, j))));
        }
      }

      // scale with scaling factors calculated from runForwardAlg
      scalingFactor = gsl_vector_get(scalingVec, i+1);
      //scalingFactor = 1;//gsl_vector_get(scalingVec, i+1);

      // scale this->currBackwardCol by scalingFactor
      gsl_vector_scale(this->currBackwardCol, 1.0 / scalingFactor);

      // save this->currBackwardCol
      gsl_matrix_set_col(backwardMat, i, this->currBackwardCol);
      //printMatrix(backwardMat);
    }

    // for the very first col, need to mult steady state, emission prob, and backwards mat. however, the backwardMat should not actually be set (see https://web.stanford.edu/~jurafsky/slp3/A.pdf page 12)
    //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[0];
    gsl_matrix_get_col(this->currBackwardCol, backwardMat, 0);

    //scalingFactor = gsl_vector_get(scalingVec, 1); // Wed 03 Jun 2020 12:41:44 PM PDT added this scalingFactor, not sure why it wasn't there before?
    scalingFactor = gsl_vector_get(scalingVec, 0); // Wed 17 Jun 2020 03:13:29 PM PDT changed back to 1:1 mapping of times and indices
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      ////emissionProb = 1;
      //emissionProb = 0;

      //// for each cell, calc the emission prob and multiply by it
      //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      //  currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
      //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[0];
      //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, windowIdx, cellIdx);
      //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx);
      //  //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx));
      //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, 0, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT
      //}
      //std::cout << emissionProb - this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0) << std::endl;
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);

      initProb = gsl_vector_get(this->initProb, stateIdx);
      //currRes = emissionProb * initProb * gsl_matrix_get(backwardMat, stateIdx, 1);
      //currRes = emissionProb * initProb * gsl_matrix_get(backwardMat, stateIdx, 1) / scalingFactor; // Wed 03 Jun 2020 12:41:44 PM PDT added this scalingFactor, not sure why it wasn't there before?
      //currRes = emissionProb * initProb * gsl_matrix_get(backwardMat, stateIdx, 0) / scalingFactor; // Wed 17 Jun 2020 03:13:29 PM PDT changed back to 1:1 mapping of times and indices
      currRes = exp(emissionProb + log(initProb) + log(gsl_matrix_get(backwardMat, stateIdx, 0)) - log(scalingFactor)); // Wed 17 Jun 2020 03:13:29 PM PDT changed back to 1:1 mapping of times and indices
      gsl_vector_set(this->currBackwardCol, stateIdx, currRes);
    }
    //gsl_matrix_set_col(backwardMat, 0, this->currBackwardCol); // according to https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Backward_procedure and stanford link earlier, don't save this
    //printMatrix(backwardMat);

    // return unscaled backward prob (prob of observed seq given model params)
    sumLogScalingFactors = 0;
    //for(i = 1; i < (int) scalingVec->size; i++) {
    for(i = 0; i < (int) scalingVec->size; i++) {
      //std::cout << "backward scalingVec[" << i << "]: " << gsl_vector_get(scalingVec, i) << ", log(): " << log(gsl_vector_get(scalingVec, i)) << std::endl;
      sumLogScalingFactors = sumLogScalingFactors + log(gsl_vector_get(scalingVec, i));
    }

    // using same logic as in runForwardAlg, sum up log of the left most col to get the total loglikelihood
    sumLastCol = gsl_blas_dasum(this->currBackwardCol); // Double Absolute SUM
    logScaledSumLastCol = log(sumLastCol);
    totalLogLikelihood += sumLogScalingFactors + logScaledSumLastCol;
  }

  if(isnan(totalLogLikelihood) || isnan(-totalLogLikelihood)) {
    totalLogLikelihood = GSL_NAN;
  }

  return totalLogLikelihood;
}

/*
 * first calls forward alg, then backward alg, then fills entries in forBackMatVec (specifically, stores the normalized log likelihood)
 */
int HMM::runForBackAlg() {
  //std::cout << "###" << this->runForwardAlg() << std::endl;
  //std::cout << "###" << this->runBackwardAlg()<< std::endl;
  double forLogLik = this->runForwardAlg();
  double backLogLik = this->runBackwardAlg();
  if(isnan(forLogLik) || isnan(backLogLik)) {
    fprintf(stderr, "ERROR: forward loglikelihood (%.10f) or backward loglikelihood (%.10f) is nan. Stopping.\n", forLogLik, backLogLik);
    return -1;
  }
  if(!compareDoubles(forLogLik, backLogLik)) {
    fprintf(stderr, "ERROR: forward loglikelihood (%.10f) does not match backward loglikelihood (%.10f). Stopping.\n", forLogLik, backLogLik);
    std::cerr << "this->transition" << std::endl;
    printMatrix(stderr, this->transition);
    return -1;
    //exit(1);
  }
  //std::cout << "###" << forLogLik << std::endl;
  //std::cout << "###" << backLogLik << std::endl;

  gsl_vector* scalingVec = nullptr;
  gsl_matrix* forwardMat = nullptr;
  gsl_matrix* backwardMat = nullptr;
  gsl_matrix* forBackMat = nullptr;

  unsigned int stateIdx = -1;
  //double forScalingFactor = 1;
  //double backScalingFactor = 1;
  double scalingFactor = 0;

  double forLik = -1;
  double backLik = -1;
  double forBackLik = -1;
  double chrTotLik = -1;

  unsigned int i = 0;
  std::vector<std::string>* chrVec = this->getChrVec();
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    scalingVec = (*this->scalingVecVec)[chrIdx];
    forwardMat = (*this->forwardMatVec)[chrIdx];
    backwardMat = (*this->backwardMatVec)[chrIdx];
    forBackMat = (*this->forBackMatVec)[chrIdx];
    gsl_matrix_set_all(forBackMat, 0);
    //std::cout << "forwardMat, chrIdx: " << chrIdx << std::endl;
    //printMatrix(forwardMat);
    //std::cout << "backwardMat, chrIdx: : " << chrIdx << std::endl;
    //printMatrix(backwardMat);

    //forScalingFactor = 1;
    //backScalingFactor = 1;
    //forScalingFactor = 0;
    //backScalingFactor = 0;
    // recalc scaling factors for backwardMat, based on scalingVec
    scalingFactor = 0;
    for(i = 0; i < forBackMat->size2; i++) {
      //backScalingFactor *= gsl_vector_get(scalingVec, i);
      //backScalingFactor += log(gsl_vector_get(scalingVec, i));
      scalingFactor += log(gsl_vector_get(scalingVec, i));
    }
    //std::cout << "starting backScalingFactor: " << backScalingFactor << std::endl;

    // recalc per chr likelihood by summing the last col of forwardMat, added to the summed and logged scaling factors
    gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
    //chrTotLik = backScalingFactor + log(gsl_blas_dasum(&lastCol.vector));
    chrTotLik = scalingFactor + log(gsl_blas_dasum(&lastCol.vector));
    //std::cout << "chrTotLike: " << chrTotLik << std::endl;

    // iter through observed sequence
    for(i = 0; i < forBackMat->size2; i++) {
      //forScalingFactor *= gsl_vector_get(scalingVec, i);
      //backScalingFactor /= gsl_vector_get(scalingVec, i);
      //forScalingFactor += log(gsl_vector_get(scalingVec, i));
      //backScalingFactor -= log(gsl_vector_get(scalingVec, i));
      //std::cout << "forScaling: " << forScalingFactor << ", backScaling: " << backScalingFactor << std::endl;

      // iter over states in col i
      for(stateIdx = 0; stateIdx < forBackMat->size1; stateIdx++) {
        //forLik = gsl_matrix_get(forwardMat, stateIdx, i) * forScalingFactor;
        //backLik = gsl_matrix_get(backwardMat, stateIdx, i) * backScalingFactor;
        //forLik = log(gsl_matrix_get(forwardMat, stateIdx, i)) + forScalingFactor;
        //backLik = log(gsl_matrix_get(backwardMat, stateIdx, i)) + backScalingFactor;
        forLik = log(gsl_matrix_get(forwardMat, stateIdx, i));
        backLik = log(gsl_matrix_get(backwardMat, stateIdx, i));
        //std::cout << "forLike: " << forLik << ", backLike: " << backLik << ", prod: " << forLik * backLik << std::endl;
        //std::cout << "logfor: " << log(gsl_matrix_get(forwardMat, stateIdx, i)) << ", forScaling: " << forScalingFactor << ", logback: " << log(gsl_matrix_get(backwardMat, stateIdx, i)) << ", backScaling: " << backScalingFactor << ", scalingProd: " << scalingProd << std::endl;
        //forBackLik = forLik * backLik / exp(totalLogLike);
        /*if(isinf(forLik)) {
          forLik = std::numeric_limits<double>::min();
        }
        if(isinf(backLik)) {
          backLik = std::numeric_limits<double>::min();
        }*/
        //std::cout << "for: " << forLik << ", mat: " << gsl_matrix_get(forwardMat, stateIdx, i) << ", back: " << backLik << ", chrT: " << chrTotLik << ", scal: " << scalingFactor << std::endl;
        forBackLik = forLik + backLik - chrTotLik + scalingFactor;
        //gsl_matrix_set(forBackMat, stateIdx, i, exp(forBackLike));
        gsl_matrix_set(forBackMat, stateIdx, i, forBackLik);
      }
    }
    //std::cout << "forBackMat:" << std::endl;
    //printMatrix(forBackMat);
  }
  return 0;
}

/*
 * run Baum Welch EM algorithm, using notation from wikipedia: https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Update
 */
void HMM::runBaumWelch(int numBWIters) {
  // emission prob related variables
  int numStates = this->states->size();
  //int currPloidy = -1;
  int cellIdx = -1;
  int stateIdx = -1;
  //double currDiploidDepth = -1;
  //double currTumorDepth = -1;
  double emissionProb = -1;
  gsl_vector* colEmissionProbs = gsl_vector_alloc(numStates); ;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  std::string currChr;
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  std::vector<std::string>* chrVec = this->getChrVec();
  //double nEps = firstDepthPair->numWindows * 272.5568 / this->DEPTH_ERROR_SCALING; // unnormDiploid error term; for statDist update of lib sizes
  //double weightedAvgPloidy = 0;

  gsl_vector* scalingVec = nullptr;
  gsl_matrix* forwardMat = nullptr; // scaled (by scalingVec) forward probs
  gsl_matrix* backwardMat = nullptr; // scaled (by scalingVec) backward probs
  gsl_matrix* forBackMat = nullptr; // log of matrix gamma, prob of being in state i at time t
  gsl_matrix* edgeMat = nullptr; // matrix xi, prob of being in state i at t and j at t+1
  std::vector<gsl_matrix*>* edgeMatVec = new std::vector<gsl_matrix*>(chrVec->size());

  // intermediates
  unsigned int t = 0;
  unsigned int chrIdx = 0;
  int i = -1;
  int j = -1;
  double alpha_it = -1;
  double a_ij = -1;
  double beta_jt1 = -1;
  double b_jt1 = -1;
  double currScalingVecTot = 0;
  double forScalingFactor = 0;
  double backScalingFactor = 0;
  double chrTotalLik = 0;
  double currTotalLik = 0;
  double prevTotalLik = this->getLogLikelihood(); //-std::numeric_limits<double>::max();
  double initTotalLik = prevTotalLik;
  int forBackStatus = 0;

  // update calculation related variables
  gsl_matrix* updatedTransition = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_vector* updatedInitProb = gsl_vector_alloc(numStates);
  gsl_matrix* summedEdgeMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_vector* summedForBackVec = gsl_vector_alloc(numStates); // from t=1 to T-1, for updating transition probs
  //gsl_vector* estLibScaling = gsl_vector_alloc(this->NUM_CELLS); // one entry per cell
  double currEdge = 0;
  double currForBack = 0;
  double currTransitionProb = 0;
  double currTransitionRowSum = 0;
  //double currSA = 0;
  bool transitionUpdateErr = false;

  // variables for comparing against simulation transition matrix, if have it
  gsl_matrix* simTrMat = nullptr;
  gsl_vector* simTrParams = nullptr;
  gsl_matrix* prevTrMat = nullptr;
  double currChiSq = 0;
  if(this->simParamsToEst != nullptr) {
    simTrMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
    simTrParams = gsl_vector_alloc(this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
    int simTrIdx = 0;
    for(int simParamIdx = 0; simParamIdx < this->NUM_TRANSITION_PARAMS_TO_EST; simParamIdx++, simTrIdx++) {
      gsl_vector_set(simTrParams, simTrIdx, gsl_vector_get(this->simParamsToEst, this->TRANSITION_PROB_START_IDX + simParamIdx));
    }
    for(int simParamIdx = 0; simParamIdx < this->NUM_BRANCH_LENGTHS_TO_EST; simParamIdx++, simTrIdx++) {
      gsl_vector_set(simTrParams, simTrIdx, gsl_vector_get(this->simParamsToEst, this->BRANCH_LENGTH_START_IDX + simParamIdx));
    }
    this->setTransition(simTrMat, simTrParams);
    prevTrMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
    gsl_matrix_memcpy(prevTrMat, this->transition);
  }

  // timing variables
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;
  double elapsedSec = 0;
  double totalTime = 0;
  int countTooClose = 0;

  printf("BAUM WELCH INITIAL LIKELIHOOD: %.40f\n", initTotalLik);
  for(int bwIters = 0; bwIters < numBWIters; bwIters++) {
    //std::cout << "ON bwIters: " << bwIters << std::endl;
    begin = std::chrono::steady_clock::now();
    // #### calculate normalized loglikelihoods into this->forBackMatVec #####
    forBackStatus = this->runForBackAlg();
    if(forBackStatus != 0) {
      fprintf(stderr, "ERROR: problem wtih forward/backward algorithm. Breaking out of baum welch\n");
      break;
    }

    // #### calculate intermediates #####
    // for each chr
    for(chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
      currChr = (*chrVec)[chrIdx];
      scalingVec = (*this->scalingVecVec)[chrIdx];
      forwardMat = (*this->forwardMatVec)[chrIdx];
      backwardMat = (*this->backwardMatVec)[chrIdx];

      // recalc per chr likelihood by summing the last col of forwardMat, added to the summed and logged scaling factors
      currScalingVecTot = 0;
      forScalingFactor = 0;
      backScalingFactor = 0;
      for(t = 0; t < scalingVec->size; t++) {
        backScalingFactor += log(gsl_vector_get(scalingVec, t));
        currScalingVecTot += log(gsl_vector_get(scalingVec, t));
      }
      backScalingFactor -= log(gsl_vector_get(scalingVec, 0)); // backwardMat's only useful from t+1

      gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
      chrTotalLik = currScalingVecTot + log(gsl_blas_dasum(&lastCol.vector));
 
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currDepths = (*this->depths)[cellIdx];
        (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
      }
      (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];


      if((*edgeMatVec)[chrIdx] == nullptr) {
        edgeMat = gsl_matrix_alloc(numStates * numStates, forwardMat->size2 - 1); // -1 because counting edges between cols
        (*edgeMatVec)[chrIdx] = edgeMat;
      }
      else {
        edgeMat = (*edgeMatVec)[chrIdx];
      }
      gsl_matrix_set_zero(edgeMat);

      // iter over observed seq
      for(t = 0; t < edgeMat->size2; t++) {
        // pre calc emission prob entries b_j(y_(t+1)) for a col t+1
        gsl_vector_set_zero(colEmissionProbs);
        //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[t+1];
        for(stateIdx = 0; stateIdx < numStates; stateIdx++) {
          //emissionProb = 1;
          //emissionProb = 0;
          //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          //  currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
          //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[t+1];
          //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx);
          //  //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx));
          //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, 0, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT
          //  //std::cout << currPloidy << ", " << currTumorDepth << ", " << currDiploidDepth << ", " << cellIdx << ", " << t+1 << std::endl;
          //}
          //std::cout << emissionProb - this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, t+1) << std::endl;
          emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, t+1);
          gsl_vector_set(colEmissionProbs, stateIdx, emissionProb);
          //gsl_vector_set(colEmissionProbs, stateIdx, exp(emissionProb));
        }
        //std::cout << "colEmissionProbs" << std::endl;
        //printColVector(colEmissionProbs);
        forScalingFactor += log(gsl_vector_get(scalingVec, t));
        backScalingFactor -= log(gsl_vector_get(scalingVec, t+1));

        // calc unscaled entries of edgeMat
        for(i = 0; i < numStates; i++) { // go down a col
          for(j = 0; j < numStates; j++) {
            ////alpha_it = gsl_matrix_get(forwardMat, i, t); // scaled alpha_i(t)
            //alpha_it = gsl_matrix_get(forwardMat, i, t) * exp(forScalingFactor);
            //a_ij = gsl_matrix_get(this->transition, i, j);
            ////beta_jt1 = gsl_matrix_get(backwardMat, j, t+1); // scaled beta_j(t+1)
            //beta_jt1 = gsl_matrix_get(backwardMat, j, t+1) * exp(backScalingFactor);
            //b_jt1 = gsl_vector_get(colEmissionProbs, j);

            a_ij = log(gsl_matrix_get(this->transition, i, j));
            if(isinf(a_ij)) {
              gsl_matrix_set(edgeMat, i * numStates + j, t, 0);
              continue;
            }
            alpha_it = log(gsl_matrix_get(forwardMat, i, t)) + forScalingFactor;
            beta_jt1 = log(gsl_matrix_get(backwardMat, j, t+1)) + backScalingFactor;
            //b_jt1 = log(gsl_vector_get(colEmissionProbs, j));
            b_jt1 = gsl_vector_get(colEmissionProbs, j);
            //std::cout << alpha_it << ", " << a_ij << ", " << beta_jt1 << ", " << b_jt1 << std::endl;
            //std::cout << alpha_it << ", " << a_ij << ", " << beta_jt1 << ", " << b_jt1 << "; resacled alpha: " << alpha_it * exp(forScalingFactor) << ", scaled beta: " << beta_jt1 * exp(backScalingFactor) << std::endl;
            //gsl_matrix_set(edgeMat, i * numStates + j, t, alpha_it * a_ij * beta_jt1 * b_jt1);
            //gsl_matrix_set(edgeMat, i * numStates + j, t, alpha_it * a_ij * beta_jt1 * b_jt1 * exp(forScalingFactor) * exp(backScalingFactor));
            //gsl_matrix_set(edgeMat, i * numStates + j, t, log(alpha_it) + log(a_ij) + log(beta_jt1) + log(b_jt1) + forScalingFactor + backScalingFactor);
            //gsl_matrix_set(edgeMat, i * numStates + j, t, exp(alpha_it + a_ij + beta_jt1 + b_jt1));
            //gsl_matrix_set(edgeMat, i * numStates + j, t, alpha_it + a_ij + beta_jt1 + b_jt1);
            gsl_matrix_set(edgeMat, i * numStates + j, t, exp(alpha_it + a_ij + beta_jt1 + b_jt1 - chrTotalLik));
            //gsl_matrix_set(edgeMat, i * numStates + j, t, gsl_matrix_get(forwardMat, i, t) * gsl_matrix_get(this->transition, i, j) * gsl_matrix_get(backwardMat, j, t+1) * exp(forScalingFactor + backScalingFactor + gsl_vector_get(colEmissionProbs, j) - chrTotalLik));
            /*if(compareDoubles(gsl_matrix_get(edgeMat, i * numStates + j, t), 1)) {
              std::cout << "EDGEMAT ENTRY 1 (" << chrIdx << ":" << i*numStates + j<< "=" << i << "->" << j << "," << t << "): " << alpha_it << ", " << a_ij << ", " << beta_jt1 << ", " << b_jt1 << ", " << chrTotalLik << ", sum: " << alpha_it + a_ij + beta_jt1 + b_jt1 - chrTotalLik << ", exp: " << exp(alpha_it + a_ij + beta_jt1 + b_jt1 - chrTotalLik)<< std::endl;
            }*/
          }
        }
        //std::cout << "unscaled edgeMat:" << std::endl;
        //printMatrix(edgeMat);

        /*// Tue 23 Jun 2020 12:40:51 PM PDT TODO is this normalization wrong??? scale by chrTotLik, as in runForBackAlg???
        // normalize each col so it sums to 1
        gsl_vector_view currCol = gsl_matrix_column(edgeMat, t);
        sumCurrCol = gsl_blas_dasum(&currCol.vector);
        //gsl_vector_scale(&currCol.vector, 1.0 / sumCurrCol);
        gsl_vector_scale(&currCol.vector, 1.0 / exp(chrTotalLik)); // chrTotLik is in log space
        //std::cout << "exp(chrTotalLik): " << exp(chrTotalLik) << ", log: " << chrTotalLik << std::endl;
        //std::cout << "SUMCURRCOL from scaling edgeMat: " << sumCurrCol << std::endl;*/
      //std::cout << "scaled edgeMat:" << std::endl;
      //printMatrix(edgeMat);
      }
      //std::cout << "chr edgeMat: " << chrIdx << std::endl;
      //printMatrix(edgeMat);
    } // chr loop

    // ##### update step #####
    // calculate updated transition matrix and initProb vectors (set transition matrix directly as setTransition(mat) will set initProb to stationary dist)
    // summed update equations from https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Multiple_sequences (each chr is a "sequence")
    gsl_matrix_set_zero(updatedTransition);
    gsl_vector_set_zero(updatedInitProb);
    gsl_matrix_set_zero(summedEdgeMat);
    gsl_vector_set_zero(summedForBackVec);
    //gsl_vector_set_zero(estLibScaling);

    // for each chr
    for(chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
      currChr = (*chrVec)[chrIdx];
      edgeMat = (*edgeMatVec)[chrIdx];
      //printMatrix(edgeMat);
      forBackMat = (*this->forBackMatVec)[chrIdx];

      // calculate numerator of a_ij* (ie sum_t=1^T-1 xi_ij(t))
      // sum up the edge mat
      for(i = 0; i < numStates; i++) {
        for(j = 0; j < numStates; j++) {
          currEdge = gsl_matrix_get(summedEdgeMat, i, j);

          // iter over observed seq
          for(t = 0; t < edgeMat->size2; t++) {
            currEdge += gsl_matrix_get(edgeMat, i * numStates + j, t);
            //std::cout << gsl_matrix_get(edgeMat, i * numStates + j, t) << ", " << currEdge << std::endl;
          }
          gsl_matrix_set(summedEdgeMat, i, j, currEdge);
        }
      }
      //std::cout << "summedEdgeMat, chrIdx: " << chrIdx << std::endl;
      //printMatrix(summedEdgeMat);

      //std::cout << "forBackMat, chrIdx: " << chrIdx << std::endl;
      //printMatrix(forBackMat);

      // calculate denominator of a_ij* (ie sum_t=1^T-1 gamma_i(t))
      // sum up the forBack mat (ie total exp time in state i)
      for(i = 0; i < numStates; i++) {
        currForBack = gsl_vector_get(summedForBackVec, i);
        // iter over observed seq (to t-1 to match edgeMat)
        for(t = 0; t < forBackMat->size2-1; t++) {
          currForBack += exp(gsl_matrix_get(forBackMat, i, t)); // forBackMat stores loglikelihoods
        }
        gsl_vector_set(summedForBackVec, i, currForBack);

        // also sum first col into updatedInitProb
        currForBack = gsl_vector_get(updatedInitProb, i);
        currForBack += exp(gsl_matrix_get(forBackMat, i, 0)); // forBackMat stores loglikelihoods
        //std::cout << "initProb adding: " << gsl_matrix_get(forBackMat, i, 0) << ", exp(): " << exp(gsl_matrix_get(forBackMat, i, 0)) << std::endl;
        gsl_vector_set(updatedInitProb, i, currForBack);
      }
      //std::cout << "summedForBackVec, chrIdx: " << chrIdx << std::endl;
      //printRowVector(summedForBackVec);

      /*// calcualte numerator of b_i* (ie sum_t=1^T 1_emitSym * gamma_i(t), but 1_emitSym indicator is instead s_A = # reads in window / (currPloidy *diploid/2 + epsilon))
      // i = ploidy, A = cell, t = window
      // s_Ait = # reads in window t, cell A / (currPloidy_t(i) * diploid_t /2 + epsilon))
      // estLibScaling_A = sum_i sum_t s_Ait * gamma_i(t) / (sum_i sum_t gamma_i(t))
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currDepths = (*this->depths)[cellIdx];
        (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
      }
      (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

      // iter over observed seq
      for(t = 0; t < forBackMat->size2; t++) {
        currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[t];
        for(stateIdx = 0; stateIdx < numStates; stateIdx++) {
          for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
            currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
            currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[t];
            //double lambda_i = (ploidy * diploidDepth/2.0 + 272.5568 / this->DEPTH_ERROR_SCALING) * currLibSizeScalingFactor; // withErr, unnormDiploid
            currSA = currTumorDepth / (currPloidy * currDiploidDepth / 2.0 + 272.5568 / this->DEPTH_ERROR_SCALING); // TODO copied over from TwoCell3TrParam2DegPolyHMM getEmissionProb
            //std::cout << cellIdx << ", " << currSA << ", " << currTumorDepth << ", " << t << ", " << currPloidy << ", " << currDiploidDepth << ", " << exp(gsl_matrix_get(forBackMat, stateIdx, t)) << ", " << currSA * exp(gsl_matrix_get(forBackMat, stateIdx, t))<< std::endl;
            gsl_vector_set(estLibScaling, cellIdx, gsl_vector_get(estLibScaling, cellIdx) + currSA * exp(gsl_matrix_get(forBackMat, stateIdx, t))); // forBackMat stores loglikelihoods
          }
          //printRowVector(estLibScaling);
        }
        printRowVector(estLibScaling);
        std::cout << std::endl;
      }*/
      //printColVector(estLibScaling);
    } // chr loop
    //std::cout << "summedEdgeMat" << std::endl;
    //printMatrix(summedEdgeMat);

    //std::cout << "summedForBackVec" << std::endl;
    //printColVector(summedForBackVec);

    // update transition matrix
    for(i = 0; i < numStates; i++) {
      currForBack = gsl_vector_get(summedForBackVec, i);
      currTransitionRowSum = 0;
      for(j = 0; j < numStates; j++) {
        currEdge = gsl_matrix_get(summedEdgeMat, i, j);
        currTransitionProb = currEdge / currForBack;
        //printf("currTransitionProb: %.10f = %.10f / %.10f\n", currTransitionProb, currEdge, currForBack);
        if(isnan(currTransitionProb)) {
          std::cerr << "ERROR in Baum Welch: nan in transition matrix detected (" << i << ", " << j  << ") ==> " << currEdge << " / " << currForBack << ", breaking. Curr transition matrix is:" << std::endl;
          printMatrix(stderr, this->transition);
          std::cout << "ERROR in Baum Welch: nan in transition matrix detected (" << i << ", " << j  << ") ==> " << currEdge << " / " << currForBack << ", breaking. Curr transition matrix is:" << std::endl;
          printMatrix(stdout, this->transition);
          transitionUpdateErr = true;
          break;
        }
        gsl_matrix_set(updatedTransition, i, j, currTransitionProb);
        currTransitionRowSum += currTransitionProb;
      }
      //printf("currTransitionRowSum: %.10f\n", currTransitionRowSum);

      // safety check: divide by rowsum to ensure transition matrix entries stay in [0,1], but only if currTransitionRowSum is large enough
      //if(currTransitionRowSum > 1e-6) {
        for(j = 0; j < numStates; j++) {
          currTransitionProb = gsl_matrix_get(updatedTransition, i, j);
          gsl_matrix_set(updatedTransition, i, j, currTransitionProb / currTransitionRowSum);
        }
      //}
    }
    if(transitionUpdateErr) {
      break;
    }
    gsl_matrix_memcpy(this->transition, updatedTransition);
    //std::cout << "updated transition:" << std::endl;
    //printMatrix(this->transition);

    // update initProb
    gsl_vector_scale(updatedInitProb, 1.0 / chrVec->size());
    //this->findSteadyStateDist(updatedInitProb);
    gsl_vector_memcpy(this->initProb, updatedInitProb);
    //std::cout << "updated initProb" << std::endl;
    //printColVector(this->initProb);

    /*// update emission probs
    std::cout << "updated libs" << std::endl;
    // scale each cell's entry by summed forBack across all states and genomic positions
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currSA = gsl_vector_get(estLibScaling, cellIdx) / firstDepthPair->numWindows; // sum across all genomic positions for all states of the forBack (gamma) matrix is the same as the sum of all cols in that matrix. Each col sums to 1 (by definition) and so this is just the number of all genomic positions.
      this->setLibScalingFactor(cellIdx, currSA);
      std::cout << this->getLibScalingFactor(cellIdx) << std::endl;
    }
    std::cout << std::endl;*/

    /*// TODO Mon 22 Jun 2020 11:11:04 AM PDT check if the stat dist from BW predicts the correct number of total reads
    // stat dist libs. Wed 21 Oct 2020 07:51:39 PM PDT this doesn't work if have absorbing states
    // let X_ij ~ NB(ploidy_i * DR_i / 2 + epsilon) * c_j
    // where c_j = T_j / (n * epsilon + sum_windows (DR_j / 2 * weighted avg ploidy (by stationary dist) + epsilon))
    // get weighted avg ploidy
    this->findSteadyStateDist(updatedInitProb);
    //std::cout << "updated steady state dist:" << std::endl;
    //printColVector(updatedInitProb);
    //std::cout << "updated libs" << std::endl;
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      weightedAvgPloidy = 0;
      for(stateIdx = 0; stateIdx < (int) updatedInitProb->size; stateIdx++) {
        currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
        weightedAvgPloidy += currPloidy * gsl_vector_get(updatedInitProb, stateIdx);
      }
      currSA = currDepths->tumorLibrarySize / (currDepths->diploidLibrarySize / 2.0 * weightedAvgPloidy + nEps);

      //this->setLibScalingFactor(cellIdx, currSA);
      std::cerr << bwIters << "\tupdatedLib_" << cellIdx << "_stat\t" << currSA << std::endl;

      //double estT = this->getLibScalingFactor(cellIdx) * (currDepths->diploidLibrarySize / 2.0 * weightedAvgPloidy + nEps);
      //std::cout << "######### storedLib: " << this->getLibScalingFactor(cellIdx) << ", estT: " << estT << ", trueT: " << currDepths->tumorLibrarySize << std::endl;
    }*/
    // viterbi libs
    this->estLibScalingFactorsPosterior();
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      std::cerr << bwIters << "\tupdatedLib_" << cellIdx << "_vit\t" << this->getLibScalingFactor(cellIdx) << std::endl;
    }

    // eval how well this step worked
    // compare updated transition mat to orig mat
    if(simTrMat != nullptr) {
      currChiSq = calcChiSqOfMatrix(simTrMat, this->transition);
      //std::cout << "Curr chiSq between simulation matrix and current transition matrix: " << currChiSq << std::endl;
      //std::cerr << bwIters << "\tCurr chiSq between simulation matrix and current transition matrix:\t" << currChiSq << std::endl;
      std::cerr << bwIters << "\tcurrChiSq_sim_curr\t" << currChiSq << std::endl;

      currChiSq = calcChiSqOfMatrix(prevTrMat, this->transition);
      //std::cout << "Curr chiSq between previous transition matrix and current transition matrix: " << currChiSq << std::endl;
      //std::cerr << bwIters << "\tCurr chiSq between previous transition matrix and current transition matrix:\t" << currChiSq << std::endl;
      std::cerr << bwIters << "\tcurrChiSq_prev_curr\t" << currChiSq << std::endl;
      gsl_matrix_memcpy(prevTrMat, this->transition);
    }

    // check change in loglikelihood
    currTotalLik = this->getLogLikelihood();
    //std::cout << "currTotalLik: " << currTotalLik << std::endl;
    //std::cout << "diffTotalLik: " << currTotalLik - prevTotalLik << std::endl;
    //std::cout << bwIters << "\tcurrTotalLik:\t" << currTotalLik << std::endl;
    //std::cout << bwIters << "\tdiffTotalLik:\t" << currTotalLik - prevTotalLik << std::endl;
    //std::cout << "End of BW iter " << bwIters << std::endl;

    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    totalTime += elapsedSec;
    printf("ON BW ITER %d, likelihood %.20f, time elapsed (sec) %.5f, iter change in loglikelihood %.5f, total change in loglikelihood %.5f\n", bwIters, currTotalLik, elapsedSec, currTotalLik - prevTotalLik, currTotalLik - initTotalLik);

    if(std::abs(currTotalLik - prevTotalLik) < 1e-4) {
      countTooClose++;
      if(countTooClose >= 2) {
        std::cout << "STATUS IS: consecutive small changes in loglikelihood detected, breaking" << std::endl;
        break;
      }
    } else {
      countTooClose = 0;
    }
    prevTotalLik = currTotalLik;
  } // baum welch loop

  printf("BAUM WELCH INDV HMM LOGLIKELIHOOD FOUND: %.40f\n", currTotalLik);
  printf("CHANGE IN BAUM WELCH INDV HMM LIKELIHOOD: %.40f\n", currTotalLik - initTotalLik);
  printf("BAUM WELCH INDV HMM TIME (sec): %.10f\n", totalTime);
  printf("DONE WITH BAUM WELCH INDV HMM.\n\n");

  // clean up
  gsl_vector_free(colEmissionProbs);
  gsl_matrix_free(updatedTransition);
  gsl_vector_free(updatedInitProb);
  gsl_matrix_free(summedEdgeMat);
  gsl_vector_free(summedForBackVec);
  gsl_matrix_free(simTrMat);
  gsl_vector_free(simTrParams);
  gsl_matrix_free(prevTrMat);
  for(std::vector<gsl_matrix*>::iterator itr = edgeMatVec->begin(); itr != edgeMatVec->end(); ++itr) {
    gsl_matrix_free(*itr);
  }
  delete edgeMatVec;
  //gsl_vector_free(estLibScaling);
}

/*
 * Viterbi algorithm for decoding all cells
 * sets chrToViterbiPathMap for DepthPairs
 *
 * //returns number of state transitions decoded
 */
void HMM::viterbiDecode() {
  gsl_matrix_set_zero(this->numTransitionsMat);
  std::string currChr;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;

  gsl_matrix* backTrace = nullptr;
  gsl_matrix* forwardMat = nullptr;
  double scalingFactor = 0;

  //double currTumorDepth = -1;
  //double currDiploidDepth = -1;
  double emissionProb = -1;
  double currMax = -1;
  int currPloidy = -1;
  int currMaxIdx = -1;
  int stateIdx = -1;
  int cellIdx = -1;
  double currResPloidy = 0;

  unsigned int i = 0;
  unsigned int j = 0;
  int depthIdx = 0;
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    backTrace = (*this->backTraceVec)[chrIdx];
    forwardMat = (*this->forwardMatVec)[chrIdx];
    gsl_matrix_set_all(backTrace, 0);
    gsl_matrix_set_all(forwardMat, 0);

    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    //currTumorDepth = -1;
    emissionProb = 1;
    //currPloidy = -1;
    //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[0];
    // for each state in first col
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      //emissionProb = 0;

      //// for each cell, calc the emission prob and multiply by it
      //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      //  currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
      //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[0];
      //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx);
      //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, 0, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT
      //}
      //currResPloidy = emissionProb * (gsl_vector_get(this->initProb, stateIdx));
      //gsl_matrix_set(forwardMat, stateIdx, 0, currResPloidy);
      //std::cout << emissionProb - this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0) << std::endl;
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);
      currResPloidy = emissionProb + log(gsl_vector_get(this->initProb, stateIdx));
      gsl_matrix_set(forwardMat, stateIdx, 0, exp(currResPloidy));
      /*if(currRes > currMax) {
        currMax = currRes;
        currMaxIdx = i;
      }*/
    }

    /*// set first col to initProb
    gsl_matrix_set_col(forwardMat, 0, this->initProb);*/

    //std::cerr << "rescale first col by max" << std::endl;
    // rescale first col by max
    gsl_vector_view firstCol = gsl_matrix_column(forwardMat, 0);
    scalingFactor = gsl_vector_max(&firstCol.vector);
    gsl_vector_scale(&firstCol.vector, 1.0 / scalingFactor);

    // iter through observed sequence
    gsl_vector_set_zero(this->prevForwardCol);
    gsl_vector_set_zero(this->currForwardCol);
    currResPloidy = 0;
    //std::cerr << "iter through obs seq" << std::endl;
    //for(i = 1, depthIdx = 0; i < forwardMat->size2; i++, depthIdx++) {
    for(i = 1, depthIdx = 1; i < forwardMat->size2; i++, depthIdx++) {
      gsl_vector_set_zero(this->currForwardCol);
      gsl_matrix_get_col(this->prevForwardCol, forwardMat, i-1);

      //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[i];
      //currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {
        //emissionProb = 0;

        //// for each cell, calc the emission prob and multiply by it
        //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        //  currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
        //  currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
        //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
        //  //emissionProb *= this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0, cellIdx);
        //  emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, chrIdx, 0, cellIdx)); // USE THIS ONE; Wed 01 Jul 2020 02:16:18 PM PDT
        //}
        //std::cout << emissionProb - this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx) << std::endl;
        emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);

        currMax = -1;
        currMaxIdx = -1;
        // iter over possible transitions
        for(j = 0; j < this->prevForwardCol->size; j++) {
          //this->currForwardCol[stateIdx] = max(emissionProb * this->prevForwardCol[j] * this->transition[j][stateIdx]);
          //currResPloidy = gsl_vector_get(this->currForwardCol, stateIdx) + emissionProb * gsl_vector_get(this->prevForwardCol, j) * gsl_matrix_get(this->transition, j, stateIdx);
          currResPloidy = gsl_vector_get(this->currForwardCol, stateIdx) + exp(emissionProb) * gsl_vector_get(this->prevForwardCol, j) * gsl_matrix_get(this->transition, j, stateIdx);
          //std::cout << "currResPloidy: " << currResPloidy << ", emissionProb: " << emissionProb << ", prevCoL: " << gsl_vector_get(this->prevForwardCol, j) << ", transition: " << gsl_matrix_get(this->transition, j, stateIdx) << std::endl;;
          if(currResPloidy > currMax) {
            currMax = currResPloidy;
            currMaxIdx = j;
          }
        }
        gsl_vector_set(this->currForwardCol, stateIdx, currMax);
        gsl_matrix_set(backTrace, stateIdx, i, currMaxIdx);
      }

      // find max of this->currForwardCol to scale with
      scalingFactor = gsl_vector_max(this->currForwardCol);

      // scale this->currForwardCol by scalingFactor
      gsl_vector_scale(this->currForwardCol, 1.0 / scalingFactor);

      // save this->currForwardCol and scalingFactor
      gsl_matrix_set_col(forwardMat, i, this->currForwardCol);
      //printMatrix(backTrace);
    }

    //std::cerr << "recover seq" << std::endl;
    // recover sequence
    // find best state at end of sequence (ie last col)
    gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
    currMax = gsl_vector_max(&lastCol.vector);
    currMaxIdx = gsl_vector_max_index(&lastCol.vector);

    // backtrack
    //std::cout << "currMaxIdx: " << currMaxIdx << std::endl;
    std::vector<std::forward_list<int>> paths(this->NUM_CELLS, std::forward_list<int>());
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currPloidy = getCellPloidyFromStateIdx(cellIdx, currMaxIdx);
      paths[cellIdx].push_front(currPloidy);
    }
    int prevMaxIdx = currMaxIdx;
    //for(int k = (int) backTrace->size2-1; k > 1; k--) {
    for(int k = (int) backTrace->size2-1; k > 0; k--) {
      currMaxIdx = gsl_matrix_get(backTrace, currMaxIdx, k);
      //std::cerr << "currMaxIdx: " << currMaxIdx << ", k: " << k << std::endl;

      // count up number of state changes
      gsl_matrix_set(this->numTransitionsMat, currMaxIdx, prevMaxIdx, gsl_matrix_get(this->numTransitionsMat, currMaxIdx, prevMaxIdx) + 1); // going backwards, so indices are flipped
      prevMaxIdx = currMaxIdx;
      //std::cout << "next currMaxIdx: " << currMaxIdx << ", k: " << k << std::endl;
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currPloidy = getCellPloidyFromStateIdx(cellIdx, currMaxIdx);
        paths[cellIdx].push_front(currPloidy);
      }
    }
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      if(currDepths->chrToViterbiPathMap->count(currChr) > 0) {
        currDepths->chrToViterbiPathMap->erase(currChr);
      }
      (*currDepths->chrToViterbiPathMap)[currChr] = new std::vector<int>(std::begin(paths[cellIdx]), std::end(paths[cellIdx]));
    }
  }

  //delete currChrDepthsVec;
}


/*
 * function to calculate the unweighted average ploidy from the current viterbi
 * decoding
 * returns a vector of average ploidies with one entry per cell
 */
gsl_vector* HMM::getAveragePloidy() {
  this->viterbiDecode();

  gsl_vector* averagePloidies = gsl_vector_alloc(this->NUM_CELLS);
  std::vector<std::string>* chrVec = this->getChrVec();
  DepthPair* currDepths = nullptr;
  std::string currChr;
  int totalNumWindows = 0;
  int totalPloidy = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currDepths = (*this->getDepths())[cellIdx];
    totalNumWindows = 0;
    totalPloidy = 0;
    for(unsigned int i = 0; i < chrVec->size(); i++) {
      currChr = (*chrVec)[i];
      for(unsigned int regionIdx = 0; regionIdx < (*currDepths->regions)[currChr]->size(); regionIdx++) {
        totalPloidy += (*(*currDepths->chrToViterbiPathMap)[currChr])[regionIdx];
        totalNumWindows++;
      }
    }
    double averagePloidy = (double) totalPloidy / (double) totalNumWindows;
    gsl_vector_set(averagePloidies, cellIdx, averagePloidy);
  }
  return averagePloidies;
}


/*
 * function to get marginal likelihoods (in order to get most likely state at a given point in time across pairs), using the forward/backward matrix
 * First calls forward/backward alg, then for each cell, sets each DepthPair's forBackMargMatMap to the log summed likelihood
 * Of note, this unlogs the forBackMat entries, sums them up by state, then stores them as logs (ie if you unlog each cell's individual forBackMargMat, the columns will sum to 1)
 */
void HMM::calcMargLikelihoods() {
  // forward probability of ending up in any particular state given the first t observations P(data | seq_1:t)
  // backward probability of observing the remaining observations given any starting point P(data | seq_t+1:end)
  // contents of forBackMatVec is probability of being in state j at time t.// Noteabily, index 0 is the stat dist (initProb), so the 0'th observation depth occurs at index 1 in forBackMatVec // Wed 17 Jun 2020 03:13:29 PM PDT changed back to 1:1 mapping of times and indices
  this->runForBackAlg();

  // set up DepthPair forBackMargMatMap fields
  std::string currChr;
  std::vector<std::string>* chrVec = this->getChrVec();
  for(int cellNum = 0; cellNum < this->NUM_CELLS; cellNum++) {
    for(unsigned int i = 0; i < chrVec->size(); i++) {
      currChr = (*chrVec)[i];
      if((*((*this->depths)[cellNum])->forBackMargMatMap)[currChr] != nullptr) {
        gsl_matrix_free((*((*this->depths)[cellNum])->forBackMargMatMap)[currChr]);
      }
      int currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
      //gsl_matrix* currMargMat = gsl_matrix_alloc(this->MAX_PLOIDY+1, currNumWindows+1); // because this is marginalized, only need k x numWindows matrix
      gsl_matrix* currMargMat = gsl_matrix_alloc(this->MAX_PLOIDY+1, currNumWindows); // because this is marginalized, only need k x numWindows matrix
      gsl_matrix_set_zero(currMargMat);
      (*((*this->depths)[cellNum])->forBackMargMatMap)[currChr] = currMargMat;
    }
  }

  gsl_matrix* forBackMat = nullptr;
  gsl_matrix* currForBackMargMat = nullptr;
  unsigned int i = 0;
  unsigned int stateIdx = -1;

  int cellIdx = -1;
  int currPloidy = -1;

  double currForBackProb = 0;
  double currForBackMargProb = 0;

  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    forBackMat = (*this->forBackMatVec)[chrIdx];
    currChr = (*chrVec)[chrIdx];

    // for each chr position //(indexed by forBackMat, so start at 1) // Wed 17 Jun 2020 03:13:29 PM PDT changed back to 1:1 mapping of times and indices
    //for(i = 1; i < forBackMat->size2; i++) {
    for(i = 0; i < forBackMat->size2; i++) {

      // for each state
      for(stateIdx = 0; stateIdx < forBackMat->size1; stateIdx++) {
        currForBackProb = gsl_matrix_get(forBackMat, stateIdx, i);

        // for each cell
        for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
          // get cell ploidy (from state)
          currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);

          // add the curr forBackMat entry to cell's DepthPair forBackMargMatMap entry
          currForBackMargMat = (*((*this->depths)[cellIdx])->forBackMargMatMap)[currChr];
          currForBackMargProb = gsl_matrix_get(currForBackMargMat, currPloidy, i);
          //currForBackMargProb += currForBackProb; // use this one if you don't care about the marginal matrices summing to 1 in probability space, and then remove the re-logging below
          currForBackMargProb += exp(currForBackProb);
          gsl_matrix_set(currForBackMargMat, currPloidy, i, currForBackMargProb);
        }
      }
      // re-log for consistency (all probs stored as logs)
      // for each cell
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currForBackMargMat = (*((*this->depths)[cellIdx])->forBackMargMatMap)[currChr];
        // for each ploidy in that cell
        for(currPloidy = 0; currPloidy < (int) currForBackMargMat->size1; currPloidy++) {
          currForBackMargProb = gsl_matrix_get(currForBackMargMat, currPloidy, i);
          currForBackMargProb = log(currForBackMargProb);
          gsl_matrix_set(currForBackMargMat, currPloidy, i, currForBackMargProb);
        }
      }
    }

    //for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    //  std::cout << "cellIdx: " << cellIdx << std::endl;
    //  currForBackMargMat = (*((*this->depths)[cellIdx])->forBackMargMatMap)[currChr];
    //  printMatrix(currForBackMargMat);
    //}
  }
}

/*
 * function to find the steady state distribution of this HMM.
 * Assumes transtion matrix is already set, stores steady state dist
 * in passed steadStateVec
 */
double HMM::findSteadyStateDist(gsl_vector* steadyStateVec) const {
  /*raiseMatrixToPower(this->transitionTranspose, 500000, this->transition);
  gsl_vector_view firstRowView = gsl_matrix_row(this->transitionTranspose, 1);
  std::cout << "transition matrix raised to 500000'th power" << std::endl;
  printMatrix(this->transitionTranspose);
  gsl_vector_memcpy(steadyStateVec, &firstRowView.vector);
  return 0;*/


  // solve for eigenvalues and eigenvectors pi*P = pi
  // steady state dist is pi such that pi * P = pi
  // pi is the eigenvector corresponding to the eigenvalue of 1 for P^T (ie P transpose)

  // turn off error handler to get rid of "gsl: francis.c:209: ERROR: maximum iterations reached without finding all eigenvalues" message that stops program execution
  // see https://www.gnu.org/software/gsl/doc/html/err.html#c.gsl_set_error_handler_off
  gsl_error_handler_t* errHandler = gsl_set_error_handler_off();

  // get transpose of transition matrix
  gsl_matrix_transpose_memcpy(this->transitionTranspose, this->transition);

  // solve for eigenvalues and eigenvectors. transitionTranspose is nonsymmetric
  // see https://www.gnu.org/software/gsl/manual/html_node/Eigenvalue-and-Eigenvector-Examples.html#Eigenvalue-and-Eigenvector-Examples
  double status = gsl_eigen_nonsymmv(this->transitionTranspose, this->ssEval, this->ssEvec, this->ssEigenWorkspace);

  if(status != GSL_SUCCESS) {
    //return false;
    //printMatrix(this->transition);
    //printColVector(this->paramsToEst);
    //exit(-1);
    //GSL_ERROR("could not find eigenvalues for steady state dist", GSL_EMAXITER);

    std::cerr << "ERROR: could not find eigenvalues to find steady state distribution" << std::endl;
    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // get eigenvalue == 1 and corresponding eigenvector. store in steadyStateVec
  bool foundSteadyState = false;
  /*std::cout << "curr eval_i to find steady state: " << std::endl;
  for (unsigned int i = 0; i < this->ssEval->size; i++) {
    gsl_complex eval_i = gsl_vector_complex_get(this->ssEval, i);
    std::cout << GSL_REAL(eval_i) << std::endl;
  }*/
  for (unsigned int i = 0; i < this->ssEval->size; i++) {
    gsl_complex eval_i = gsl_vector_complex_get(this->ssEval, i);

    if(compareDoubles(GSL_REAL(eval_i), 1.0)) {
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column(this->ssEvec, i);
      gsl_vector_view evecReal = gsl_vector_complex_real(&evec_i.vector);

      // get absolute value to make all positive
      vectorAbsoluteValue(&evecReal.vector);

      // normalize so sum = 1
      double sum = gsl_blas_dasum(&evecReal.vector);
      gsl_vector_scale(&evecReal.vector, 1.0 / sum);
      gsl_vector_memcpy(steadyStateVec, &evecReal.vector);
      foundSteadyState = true;
      break;
    }
  }
  if(!foundSteadyState) {
    std::cerr << "ERROR: could not find steady state distribution" << std::endl;
    //exit(-1);
  }

  // restore error handler
  gsl_set_error_handler(errHandler);
  return GSL_SUCCESS;
  //return foundSteadyState;
}

/*
 * function to print elements that all HMMs have
 * ex hmm->print(stdout);
 */
void HMM::print(FILE* stream) {
  /*// depths
  int cellCtr = 0;
  for(std::vector<DepthPair*>::iterator it = this->depths->begin(); it != this->depths->end(); ++it, cellCtr++) {
    fprintf(stream, "cell pair %d:\n", cellCtr);
    (*it)->print(stream);
  }*/

  /*// states
  fprintf(stream, "states:\n");
  for(std::vector<std::string>::iterator it = this->states->begin(); it != this->states->end(); ++it) {
    fprintf(stream, "%s\t", (*it).c_str());
  }
  fprintf(stream, "\n\n");*/

  /*// alphabet
  fprintf(stream, "alphabet:\n");
  for(std::set<int>::iterator it = this->alphabet->begin(); it != this->alphabet->end(); ++it) {
    fprintf(stream, "%d\t", *it);
  }
  fprintf(stream, "\n\n");*/

  // transition matrix
  fprintf(stream, "transition matrix:\n");
  fprintf(stream, "\t");
  for(unsigned int i = 0; i < this->states->size(); i++) {
    fprintf(stream, "%s \t", ((*this->states)[i]).c_str());
  }
  fprintf(stream, "\n");
  for(unsigned int row = 0; row < this->states->size(); row++) {
    fprintf(stream, "%s \t", ((*this->states)[row]).c_str());
    for(unsigned int col = 0; col < this->states->size(); col++) {
      fprintf(stream, "%0.10g \t", gsl_matrix_get(this->transition, row, col));
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");

  /*// transition matrix string representation
  if(this->transitionStrings != nullptr) {
    fprintf(stream, "\t");
    for(unsigned int i = 0; i < this->states->size(); i++) {
      fprintf(stream, "%s\t", ((*this->states)[i]).c_str());
    }
    fprintf(stream, "\n");
    for(unsigned int row = 0; row < this->states->size(); row++) {
      fprintf(stream, "%s\t", ((*this->states)[row]).c_str());
      for(unsigned int col = 0; col < this->states->size(); col++) {
        fprintf(stream, "%s\t", (*(*this->transitionStrings)[row])[col]->c_str());
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }*/

  /*// rate matrix string representation
  if(this->rateMatrixStrings != nullptr) {
    fprintf(stream, "rate matrix symbols:\n");
    fprintf(stream, "\t");
    for(unsigned int i = 0; i < this->rateMatrixQ->size2; i++) {
      fprintf(stream, "%s \t", ((*this->adjBinLabels)[i]).c_str());
    }
    fprintf(stream, "\n");
    for(unsigned int row = 0; row < this->rateMatrixQ->size1; row++) {
      fprintf(stream, "%s \t", ((*this->adjBinLabels)[row]).c_str());
      for(unsigned int col = 0; col < this->rateMatrixQ->size2; col++) {
        fprintf(stream, "%s \t", (*(*this->rateMatrixStrings)[row])[col]->c_str());
      }
      fprintf(stream, "\n");
    }
    fprintf(stream, "\n");
  }*/
  // TODO add timeDepMatrixP printing?

  // initProb vector
  fprintf(stream, "initProb vector:\n");
  for(unsigned int i = 0; i < this->initProb->size; i++) {
    fprintf(stream, "%s\t", ((*this->states)[i]).c_str());
    fprintf(stream, "%.40f\n", gsl_vector_get(this->initProb, i));
  }
  fprintf(stream, "\n");

  /*// numTransitionsMat
  fprintf(stream, "viterbi decoded numTransitionsMat:\n");
  fprintf(stream, "\t");
  for(unsigned int i = 0; i < this->states->size(); i++) {
    fprintf(stream, "%s\t", ((*this->states)[i]).c_str());
  }
  fprintf(stream, "\n");
  for(unsigned int row = 0; row < this->states->size(); row++) {
    fprintf(stream, "%s\t", ((*this->states)[row]).c_str());
    for(unsigned int col = 0; col < this->states->size(); col++) {
      fprintf(stream, "%.0f\t", gsl_matrix_get(this->numTransitionsMat, row, col));
    }
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
  */

  // paramsToEst
  fprintf(stream, "paramsToEst:\n");
  printColVector(stream, this->paramsToEst);

  // fixedParams
  fprintf(stream, "fixedParams:\n");
  printColVector(stream, this->fixedParams);

  // likelihood
  fprintf(stream, "loglikelihood: %.10f\n\n", this->runForwardAlg());
}
/*
 * helper function to save an HMM to a file, internally calls print
 */
void HMM::saveHMMToFile(std::string filename) {
  FILE* outFile = fopen(filename.c_str(), "w");
  this->print(outFile);
  fclose(outFile);
}

/*
 * function to check if the state of this HMMs is currently ok.
 * That is, returns nan if has a completely uniform transition matrix, 0 o/w
 */
double HMM::checkStateValidity(double epsilon) const {
  return this->checkStateValidity(this->transition, epsilon);
}
double HMM::checkStateValidity(gsl_matrix* mat, double epsilon) const {
  bool isUnifMatrix = checkUniformMatrix(mat, epsilon);
  if(isUnifMatrix) {
    return GSL_NAN;
  }
  return 0;
}
/*
 * Function to check for transient states
 * (ex state 1 is only seen in 0->1->2 or 2->1->0 transitions. That is, 1 is the stepping stone between 0 and 2)
 * transient states are defined as (% time seen as stepping stone state > 80%) && (total num times in state >= 2).
 * Returns nan if transient states are found, 0 o/w
 * As of Mon 07 Jun 2021 11:04:26 AM PDT, these cutoffs are arbitrary and made up
 */
double HMM::checkForTransientStates() {
  double status = 0;
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<int>* currVitPath = nullptr;
  DepthPair* currDepths = nullptr;
  int currVitState = 0;
  int stateCount = 0;
  int pathSize = 0;
  gsl_vector* countTransientTimes = gsl_vector_alloc(this->MAX_PLOIDY + 1);
  gsl_vector* countTimesInPloidy = gsl_vector_alloc(this->MAX_PLOIDY + 1);
  //gsl_vector* countTransitionsIntoPloidy = gsl_vector_alloc(this->MAX_PLOIDY + 1);
  gsl_vector* countTransitionsOutOfPloidy = gsl_vector_alloc(this->MAX_PLOIDY + 1);
  gsl_vector* fracTransient = gsl_vector_alloc(this->MAX_PLOIDY + 1);

  this->viterbiDecode();

  // for each cell
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currDepths = (*this->depths)[cellIdx];
    gsl_vector_set_zero(countTransientTimes);
    gsl_vector_set_zero(countTimesInPloidy);
    //gsl_vector_set_zero(countTransitionsIntoPloidy);
    gsl_vector_set_zero(countTransitionsOutOfPloidy);
    gsl_vector_set_zero(fracTransient);

    // for each chr
    for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
      std::string currChr = (*chrVec)[chrIdx];
      currVitPath = (*currDepths->chrToViterbiPathMap)[currChr];
      /*currVitPath = new std::vector<int>();
      if(chrIdx == 0) {
      currVitPath->push_back(0);
      currVitPath->push_back(1);
      currVitPath->push_back(2);
      currVitPath->push_back(2);
      currVitPath->push_back(5);
      currVitPath->push_back(0);
      }
      if(chrIdx == 1) {
      currVitPath->push_back(2);
      currVitPath->push_back(1);
      currVitPath->push_back(0);
      currVitPath->push_back(2);
      currVitPath->push_back(0);
      currVitPath->push_back(0);
      }*/

      // for each window in currChr, get run length encoding
      pathSize = currVitPath->size();
      for(int regionIdx = 0; regionIdx < pathSize; regionIdx++) {
        stateCount = 1;
        currVitState = (*currVitPath)[regionIdx];
        //gsl_vector_set(countTransitionsIntoPloidy, currVitState, gsl_vector_get(countTransitionsIntoPloidy, currVitState) + 1);
        //std::cout << "regionidx: " << regionIdx << ", curr: " << (*currVitPath)[regionIdx] << ", next: " << (*currVitPath)[regionIdx + 1] << std::endl;
        while(regionIdx < pathSize - 1 && (*currVitPath)[regionIdx] == (*currVitPath)[regionIdx + 1]) {
          regionIdx++;
          stateCount++;
        }
        if(stateCount == 1 && (regionIdx - 1) >= 0 && (regionIdx + 1) < pathSize &&
          (((*currVitPath)[regionIdx - 1] == currVitState - 1 && (*currVitPath)[regionIdx + 1] == currVitState + 1) ||
          ((*currVitPath)[regionIdx - 1] == currVitState + 1 && (*currVitPath)[regionIdx + 1] == currVitState - 1))) {
          gsl_vector_set(countTransientTimes, currVitState, gsl_vector_get(countTransientTimes, currVitState) + 1);
        }
        gsl_vector_set(countTimesInPloidy, currVitState, gsl_vector_get(countTimesInPloidy, currVitState) + stateCount);
        gsl_vector_set(countTransitionsOutOfPloidy, currVitState, gsl_vector_get(countTransitionsOutOfPloidy, currVitState) + 1);
      }
      //std::cout << "next chr" << std::endl;
      //printRowVector(countTransientTimes);
      //printRowVector(countTimesInPloidy);
    }
    //std::cout << "done itering over chrs" << std::endl;

    // check if this cell has any transient states. if so, set status to GSL_NAN and break. else, continue to next cell
    gsl_vector_memcpy(fracTransient, countTransientTimes);
    //gsl_vector_div(fracTransient, countTimesInPloidy);
    gsl_vector_div(fracTransient, countTransitionsOutOfPloidy);
    std::cout << "HMM::checkForTransientStates countTransientTimes, countTimesInPloidy, countTransitionsOutOfPloidy, fracTransient" << std::endl;
    printColVector(countTransientTimes);
    printColVector(countTimesInPloidy);
    //printColVector(countTransitionsIntoPloidy);
    printColVector(countTransitionsOutOfPloidy);
    printColVector(fracTransient);
    /*int countTotalTransient = 0;
    double frac = 0;
    for(unsigned int i = 0; i < fracTransient->size; i++) {
      frac = gsl_vector_get(fracTransient, i);
      //if(!gsl_isnan(frac) && frac > 0.8 && gsl_vector_get(countTimesInPloidy, i) >= 2) {
      //if(!gsl_isnan(frac) && frac >= 0.49 && gsl_vector_get(countTimesInPloidy, i) >= 1) {
      //if(!gsl_isnan(frac) && frac > 0.5 && gsl_vector_get(countTransitionsOutOfPloidy, i) >= 2) {
      if(!gsl_isnan(frac) && frac >= 0.3 && gsl_vector_get(countTransitionsOutOfPloidy, i) >= 1) {
        //std::cout << "detected transient state at idx " << i << std::endl;
        countTotalTransient++;
        //status = GSL_NAN;
        //break;
      }
    }
    if(countTotalTransient >= 1) {
      status = GSL_NAN;
    }*/
    double sumTransient = gsl_blas_dasum(countTransientTimes); // Double Absolute SUM
    //if(sumTransient >= 5.9) { // transient4 transient threshold
    if(sumTransient >= 4.9) { // transient5
      status = GSL_NAN;
      break;
    }
    if(gsl_isnan(status)) { // break out of cell loop
      //std::cout << "breaking cell loop" << std::endl;
      break;
    }
  }
  gsl_vector_free(countTransientTimes);
  gsl_vector_free(countTimesInPloidy);
  gsl_vector_free(fracTransient);
  return status;
}

/*
 * helper function to write the table needed to plot the viterbi decoded CNA path for each tumor cell
 * coord | diploid_mean | diploid_var | tumor0 | viterbiDecoded0_0-MAX_PLOIDY | tumor1 | viterbiDecoded1_0-MAX_PLOIDY |...
 */
void HMM::saveViterbiDecodedCNA(std::string filename) {
  std::ofstream outFile(filename);
  std::string sep = "\t";
  std::string currChr;
  std::vector<std::string>* chrVec = this->getChrVec();
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;

  // write header
  outFile << "coord\tdiploid_mean\tdiploid_var";
  if(firstDepthPair->chrToDiploidSimStateMap != nullptr && firstDepthPair->chrToDiploidSimStateMap->size() > 0) {
    outFile << "\tdiploid_simState";
  }
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    outFile << "\ttumor" << cellIdx << "\tviterbiDecoded" << cellIdx << "_0-" << this->MAX_PLOIDY;
    if(firstDepthPair->chrToTumorSimStateMap != nullptr && firstDepthPair->chrToTumorSimStateMap->size() > 0) {
      outFile << "\tsimulated" << cellIdx << "_0-" << this->MAX_PLOIDY;
    }
  }
  outFile << std::endl;

  // for each chr
  for(unsigned int i = 0; i < chrVec->size(); i++) {
    currChr = (*chrVec)[i];

    // for each window in currChr
    for(unsigned int regionIdx = 0; regionIdx < (*firstDepthPair->regions)[currChr]->size(); regionIdx++) {
      // coord
      outFile << (*(*firstDepthPair->regions)[currChr])[regionIdx] << sep;

      // diploid mean
      outFile << (*(*firstDepthPair->chrToDiploidDepthMap)[currChr])[regionIdx] << sep;

      // diploid var
      outFile << (*(*firstDepthPair->chrToDiploidVarMap)[currChr])[regionIdx];

      // diploid simulated state (if set)
      if(firstDepthPair->chrToDiploidSimStateMap->size() > 0) {
        outFile << sep << (*(*firstDepthPair->chrToDiploidSimStateMap)[currChr])[regionIdx];
      }

      // for each tumor cell
      for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        currDepths = (*this->depths)[cellIdx];
        // tumor depth
        outFile << sep << (*(*currDepths->chrToTumorDepthMap)[currChr])[regionIdx];

        // viterbi decoded CNA
        outFile << sep << (*(*currDepths->chrToViterbiPathMap)[currChr])[regionIdx];

        // tumor simulated state (if set)
        if(currDepths->chrToTumorSimStateMap->size() > 0) {
          outFile << sep << (*(*currDepths->chrToTumorSimStateMap)[currChr])[regionIdx];
        }
      }

      // new line
      outFile << std::endl;
    }
  }
  outFile.close();
}

/*
 * helper method to get a state idx according to prob p in transition matrix
 * from fromStateIdx
 */
int HMM::getRandStateIdx(double p, int fromStateIdx) const {
  gsl_vector_view fromRow = gsl_matrix_row(this->transition, fromStateIdx);
  //std::cout << "\nfrom " << fromStateIdx << std::endl;
  return getRandStateIdx(p, &fromRow.vector);
}

/*
 * helper method to get a state idx according to prob p in probVec
 */
int HMM::getRandStateIdx(double p, gsl_vector* probVec) const {
  double cummSum = 0;
  int toState = 0;
  //std::cout << "p: " << p << std::endl;
  for(unsigned int i = 0; i < this->states->size(); i++) {
    if(cummSum > p) {
      break;
    }
    cummSum += gsl_vector_get(probVec, i);
    //std::cout << cummSum << std::endl;
    toState = i;
  }
  //std::cout << "to " << toState << std::endl;
  return toState;
}

