#include "HMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
HMM::HMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int numTransitionParamsToEst, int numCells, int numBranches, int maxPloidy, int numFixedTrParams, int numFixedLibs) :
                                        MAX_PLOIDY(maxPloidy),
                                        NUM_CELLS(numCells),
                                        NUM_LIBS_TO_EST(numCells - numFixedLibs),
                                        NUM_TRANSITION_PARAMS_TO_EST(numTransitionParamsToEst),
                                        NUM_BRANCH_LENGTHS_TO_EST(numBranches),
                                        NUM_FIXED_LIBS(numFixedLibs),
                                        NUM_FIXED_TRANSITION_PARAMS(numFixedTrParams),
                                        LIB_SIZE_SCALING_FACTOR_START_IDX(0),
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
  this->rateMatrixStrings = nullptr;

  this->transition = gsl_matrix_alloc(this->states->size(), this->states->size());
  gsl_matrix_set_zero(this->transition);

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
  this->meanVarianceCoefVec = nullptr;

  // params to estimate with BFGS: a library scaling factor for each cell and the transition params
  this->paramsToEst = gsl_vector_alloc(this->NUM_CELLS - numFixedLibs + this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
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

  this->numTransitionsMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_matrix_set_zero(this->numTransitionsMat);

}

// destructor
HMM::~HMM() {
  for(std::vector<DepthPair*>::iterator it = this->depths->begin(); it != this->depths->end(); ++it) {
    delete *it;
  }
  delete this->alphabet;
  gsl_matrix_free(this->transition);
  gsl_matrix_free(this->rateMatrixQ);
  gsl_matrix_free(this->timeDepMatrixP);
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
  for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
    gsl_vector_set(tmpTrProbs, tmpTrProbsIdx, gsl_vector_get(this->paramsToEst, TRANSITION_PROB_START_IDX + i));
    tmpTrProbsIdx++;
  }
  for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
    gsl_vector_set(tmpTrProbs, tmpTrProbsIdx, gsl_vector_get(this->paramsToEst, BRANCH_LENGTH_START_IDX + i));
    tmpTrProbsIdx++;
  }
  double status = this->setTransition(tmpTrProbs); // call subclass overridden version
  gsl_vector_free(tmpTrProbs);
  return status;
}

// always saves into this->transition; delegates to subclass
double HMM::setTransition(gsl_vector* transitionParams) {
  double status = this->setTransition(this->transition, transitionParams);
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
void HMM::setRateMatrixQ(double alpha, double beta, double lambda) {
  // if haven't saved the string represetation of the rate matrix yet, save it now
  bool saveStrs = false;
  std::string* currTrStr = nullptr;
  if(this->rateMatrixStrings == nullptr) {
    saveStrs = true;
    this->rateMatrixStrings = new std::vector<std::vector<std::string*>*>(this->rateMatrixQ->size1);
    for(unsigned int i = 0; i < this->rateMatrixQ->size1; i++) {
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

      // if move in tandem (include lambda term)
      // ex 0,1 -> 1,2  (+1)
      //    1,0 -> 2,1  (+1)
      //    2,3 -> 0,1  (-2)
      //    3,2 -> 1,0  (-2)
      if((frI - toI) == (frJ - toJ)){// || (toI - frI) == (toJ - frJ)) {
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
      }
      // else only one bin moves
      else if(frI == toI || frJ == toJ) {
        // move only one (include alpha term)
        if(std::abs(frI - toI) == 1 || std::abs(frJ - toJ) == 1) {
          currRate += alpha;
          if(saveStrs) {
            *currTrStr += "a+";
          }
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
  return this->MAX_PLOIDY;
}
void HMM::setInitProb(gsl_vector* initProb) {
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
  double intercept = 10.46385711652957084539; // Navin_Nature2011
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
double HMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);
  return this->setTransition();
}

double HMM::setFixedParams(gsl_vector* fixedParams) {
  gsl_vector_memcpy(this->fixedParams, fixedParams);

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
  for(int i = 0; i < this->NUM_CELLS; i++) {
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
  return;
  this->estLibScalingFactorsPosterior();
}
void HMM::estLibScalingFactorsPosterior() {
  this->viterbiDecode(); // only need to call once per pair; needed for viterbi version
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
      sum += (*diploidDepthVec)[windowIdx]/2.0 * expPloidy;
    }
  }

  // finally, multiply everything together
  double nEps = currDepths->numWindows * 272.5568 / this->DEPTH_ERROR_SCALING;
  double T = currDepths->tumorLibrarySize;
  double scalingFactor = (T - nEps) / sum; // viterbi version with const err

  if(this->gradientDebug) {
    fprintf(stderr, "%%%p\t%i\t%.2f\t%.2f\t%.2f\t%.10f\t%.10f\t%.10f\t%.10f\n", this, cellNum, nEps, T, sum, (T-nEps)/sum, T/(nEps+sum), currDepths->diploidLibrarySize, currDepths->tumorLibrarySize / currDepths->diploidLibrarySize); // viterbi version
  }
  return scalingFactor;
}

/*
 * helper function to convert cellIdx and stateIdx to a ploidy, instead of storing and looking up
 * see GetCellPloidyFromStateIdx.java for testing/development code
 * ex. getCellPloidyFromStateIdx(1, 4) if MAX_PLOIDY==3 and NUM_CELLS==2
 *     states = [(0,0), (0,1), (0,2), (0,3), (1,0), ...]
 *     returns (1,0) ==> 0
 */
int HMM::getCellPloidyFromStateIdx(int cellIdx, int stateIdx) const {
  int div = 1;
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
      transitionProb = gsl_vector_get(this->initProb, currVitState);

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

  int stateIdx = -1;
  int cellIdx = -1;
  double scalingFactor = 0;
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

    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);
      currResPloidy = emissionProb + log(gsl_vector_get(this->initProb, stateIdx));

      gsl_matrix_set(forwardMat, stateIdx, 0, exp(currResPloidy));
    }

    // rescale first col by max
    gsl_vector_view firstCol = gsl_matrix_column(forwardMat, 0);
    scalingFactor = gsl_vector_max(&firstCol.vector);
    gsl_vector_scale(&firstCol.vector, 1.0 / scalingFactor);
    gsl_vector_set(scalingVec, 0, scalingFactor);

    // iter through observed sequence
    gsl_vector_set_zero(this->prevForwardCol);
    gsl_vector_set_zero(this->currForwardCol);
    currResPloidy = 0;
    for(i = 1, depthIdx = 1; i < forwardMat->size2; i++, depthIdx++) {
      gsl_vector_set_zero(this->currForwardCol);
      gsl_matrix_get_col(this->prevForwardCol, forwardMat, i-1);

      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {
        emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);

        // iter over possible transitions
        for(j = 0; j < this->prevForwardCol->size; j++) {
          currResPloidy = gsl_vector_get(this->currForwardCol, stateIdx);
          gsl_vector_set(this->currForwardCol, stateIdx, currResPloidy + exp(emissionProb + log(gsl_vector_get(this->prevForwardCol, j)) + log(gsl_matrix_get(this->transition, j, stateIdx))));
        }
      }

      // find max of this->currForwardCol to scale with
      scalingFactor = gsl_vector_max(this->currForwardCol);

      // scale this->currForwardCol by scalingFactor
      gsl_vector_scale(this->currForwardCol, 1.0 / scalingFactor);

      // save this->currForwardCol and scalingFactor
      gsl_vector_set(scalingVec, i, scalingFactor);
      gsl_matrix_set_col(forwardMat, i, this->currForwardCol);
    }

    // return unscaled forward prob (prob of observed seq given model params)
    sumLogScalingFactors = 0;
    for(i = 0; i < scalingVec->size; i++) {
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
    logScaledSumLastCol = log(sumLastCol);
    totalLogLikelihood += sumLogScalingFactors + logScaledSumLastCol;
  }

  if(isnan(totalLogLikelihood) || isnan(-totalLogLikelihood)) {
    totalLogLikelihood = GSL_NAN;
  }

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

  int stateIdx = -1;
  int cellIdx = -1;
  double scalingFactor = 0;
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

    // set last col to 1
    currRes = 0;
    emissionProb = 1;
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      gsl_matrix_set(backwardMat, stateIdx, backwardMat->size2 - 1, 1.0); // last col is just 1's
    }

    // iter through observed sequence
    gsl_vector_set_zero(this->nextBackwardCol); // next in sequence, ie to the right/closer to the end
    gsl_vector_set_zero(this->currBackwardCol);
    currResPloidy = 0;
    for(i = backwardMat->size2 - 2, depthIdx = (*firstDepthPair->regions)[currChr]->size() - 1; i >= 0 ; i--, depthIdx--) { // -1 would get last col
      gsl_vector_set_zero(this->currBackwardCol);
      gsl_matrix_get_col(this->nextBackwardCol, backwardMat, i+1);

      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {

        // iter over possible transitions
        for(j = 0; j < (int) this->nextBackwardCol->size; j++) {
          currResPloidy = gsl_vector_get(this->currBackwardCol, stateIdx);

          emissionProb = this->getTotalLogEmissionProb(j, currChrDepthsVec, chrIdx, depthIdx);
          gsl_vector_set(this->currBackwardCol, stateIdx, currResPloidy + exp(emissionProb + log(gsl_vector_get(this->nextBackwardCol, j)) + log(gsl_matrix_get(this->transition, stateIdx, j))));
        }
      }

      // scale with scaling factors calculated from runForwardAlg
      scalingFactor = gsl_vector_get(scalingVec, i+1);

      // scale this->currBackwardCol by scalingFactor
      gsl_vector_scale(this->currBackwardCol, 1.0 / scalingFactor);

      // save this->currBackwardCol
      gsl_matrix_set_col(backwardMat, i, this->currBackwardCol);
    }

    // for the very first col, need to mult steady state, emission prob, and backwards mat. however, the backwardMat should not actually be set (see https://web.stanford.edu/~jurafsky/slp3/A.pdf page 12)
    gsl_matrix_get_col(this->currBackwardCol, backwardMat, 0);

    scalingFactor = gsl_vector_get(scalingVec, 0);
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);

      initProb = gsl_vector_get(this->initProb, stateIdx);
      currRes = exp(emissionProb + log(initProb) + log(gsl_matrix_get(backwardMat, stateIdx, 0)) - log(scalingFactor));
      gsl_vector_set(this->currBackwardCol, stateIdx, currRes);
    }

    // return unscaled backward prob (prob of observed seq given model params)
    sumLogScalingFactors = 0;
    for(i = 0; i < (int) scalingVec->size; i++) {
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
  }

  gsl_vector* scalingVec = nullptr;
  gsl_matrix* forwardMat = nullptr;
  gsl_matrix* backwardMat = nullptr;
  gsl_matrix* forBackMat = nullptr;

  unsigned int stateIdx = -1;
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

    // recalc scaling factors for backwardMat, based on scalingVec
    scalingFactor = 0;
    for(i = 0; i < forBackMat->size2; i++) {
      scalingFactor += log(gsl_vector_get(scalingVec, i));
    }

    // recalc per chr likelihood by summing the last col of forwardMat, added to the summed and logged scaling factors
    gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
    chrTotLik = scalingFactor + log(gsl_blas_dasum(&lastCol.vector));

    // iter through observed sequence
    for(i = 0; i < forBackMat->size2; i++) {

      // iter over states in col i
      for(stateIdx = 0; stateIdx < forBackMat->size1; stateIdx++) {
        forLik = log(gsl_matrix_get(forwardMat, stateIdx, i));
        backLik = log(gsl_matrix_get(backwardMat, stateIdx, i));
        forBackLik = forLik + backLik - chrTotLik + scalingFactor;

        gsl_matrix_set(forBackMat, stateIdx, i, forBackLik);
      }
    }
  }
  return 0;
}

/*
 * run Baum Welch EM algorithm, using notation from wikipedia: https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Update
 */
void HMM::runBaumWelch(int numBWIters) {
  // emission prob related variables
  int numStates = this->states->size();
  int cellIdx = -1;
  int stateIdx = -1;
  double emissionProb = -1;
  gsl_vector* colEmissionProbs = gsl_vector_alloc(numStates); ;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  std::string currChr;
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  std::vector<std::string>* chrVec = this->getChrVec();

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
  double prevTotalLik = this->getLogLikelihood();
  double initTotalLik = prevTotalLik;
  int forBackStatus = 0;

  // update calculation related variables
  gsl_matrix* updatedTransition = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_vector* updatedInitProb = gsl_vector_alloc(numStates);
  gsl_matrix* summedEdgeMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  gsl_vector* summedForBackVec = gsl_vector_alloc(numStates); // from t=1 to T-1, for updating transition probs
  double currEdge = 0;
  double currForBack = 0;
  double currTransitionProb = 0;
  double currTransitionRowSum = 0;
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
        for(stateIdx = 0; stateIdx < numStates; stateIdx++) {
          emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, t+1);
          gsl_vector_set(colEmissionProbs, stateIdx, emissionProb);
        }
        forScalingFactor += log(gsl_vector_get(scalingVec, t));
        backScalingFactor -= log(gsl_vector_get(scalingVec, t+1));

        // calc unscaled entries of edgeMat
        for(i = 0; i < numStates; i++) { // go down a col
          for(j = 0; j < numStates; j++) {
            a_ij = log(gsl_matrix_get(this->transition, i, j));
            if(isinf(a_ij)) {
              gsl_matrix_set(edgeMat, i * numStates + j, t, 0);
              continue;
            }
            alpha_it = log(gsl_matrix_get(forwardMat, i, t)) + forScalingFactor;
            beta_jt1 = log(gsl_matrix_get(backwardMat, j, t+1)) + backScalingFactor;
            b_jt1 = gsl_vector_get(colEmissionProbs, j);
            gsl_matrix_set(edgeMat, i * numStates + j, t, exp(alpha_it + a_ij + beta_jt1 + b_jt1 - chrTotalLik));
          }
        }
      }
    } // chr loop

    // ##### update step #####
    // calculate updated transition matrix and initProb vectors (set transition matrix directly as setTransition(mat) will set initProb to stationary dist)
    // summed update equations from https://en.wikipedia.org/wiki/Baum%E2%80%93Welch_algorithm#Multiple_sequences (each chr is a "sequence")
    gsl_matrix_set_zero(updatedTransition);
    gsl_vector_set_zero(updatedInitProb);
    gsl_matrix_set_zero(summedEdgeMat);
    gsl_vector_set_zero(summedForBackVec);

    // for each chr
    for(chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
      currChr = (*chrVec)[chrIdx];
      edgeMat = (*edgeMatVec)[chrIdx];
      forBackMat = (*this->forBackMatVec)[chrIdx];

      // calculate numerator of a_ij* (ie sum_t=1^T-1 xi_ij(t))
      // sum up the edge mat
      for(i = 0; i < numStates; i++) {
        for(j = 0; j < numStates; j++) {
          currEdge = gsl_matrix_get(summedEdgeMat, i, j);

          // iter over observed seq
          for(t = 0; t < edgeMat->size2; t++) {
            currEdge += gsl_matrix_get(edgeMat, i * numStates + j, t);
          }
          gsl_matrix_set(summedEdgeMat, i, j, currEdge);
        }
      }

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
        gsl_vector_set(updatedInitProb, i, currForBack);
      }
    } // chr loop

    // update transition matrix
    for(i = 0; i < numStates; i++) {
      currForBack = gsl_vector_get(summedForBackVec, i);
      currTransitionRowSum = 0;
      for(j = 0; j < numStates; j++) {
        currEdge = gsl_matrix_get(summedEdgeMat, i, j);
        currTransitionProb = currEdge / currForBack;
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

      // safety check: divide by rowsum to ensure transition matrix entries stay in [0,1]
      for(j = 0; j < numStates; j++) {
        currTransitionProb = gsl_matrix_get(updatedTransition, i, j);
        gsl_matrix_set(updatedTransition, i, j, currTransitionProb / currTransitionRowSum);
      }
    }
    if(transitionUpdateErr) {
      break;
    }
    gsl_matrix_memcpy(this->transition, updatedTransition);

    // update initProb
    gsl_vector_scale(updatedInitProb, 1.0 / chrVec->size());
    gsl_vector_memcpy(this->initProb, updatedInitProb);

    // viterbi libs
    this->estLibScalingFactorsPosterior();
    if(this->gradientDebug) {
      for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
        std::cerr << bwIters << "\tupdatedLib_" << cellIdx << "_vit\t" << this->getLibScalingFactor(cellIdx) << std::endl;
      }
    }

    // eval how well this step worked
    // compare updated transition mat to orig mat
    if(simTrMat != nullptr) {
      currChiSq = calcChiSqOfMatrix(simTrMat, this->transition);
      if(this->gradientDebug) {
        std::cerr << bwIters << "\tcurrChiSq_sim_curr\t" << currChiSq << std::endl;
      }

      currChiSq = calcChiSqOfMatrix(prevTrMat, this->transition);
      if(this->gradientDebug) {
        std::cerr << bwIters << "\tcurrChiSq_prev_curr\t" << currChiSq << std::endl;
      }
      gsl_matrix_memcpy(prevTrMat, this->transition);
    }

    // check change in loglikelihood
    currTotalLik = this->getLogLikelihood();

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
}

/*
 * Viterbi algorithm for decoding all cells
 * sets chrToViterbiPathMap for DepthPairs
 */
void HMM::viterbiDecode() {
  gsl_matrix_set_zero(this->numTransitionsMat);
  std::string currChr;
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;

  gsl_matrix* backTrace = nullptr;
  gsl_matrix* forwardMat = nullptr;
  double scalingFactor = 0;

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

    emissionProb = 1;
    // for each state in first col
    for(stateIdx = 0; stateIdx < (int) this->initProb->size; stateIdx++) {
      emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, 0);
      currResPloidy = emissionProb + log(gsl_vector_get(this->initProb, stateIdx));
      gsl_matrix_set(forwardMat, stateIdx, 0, exp(currResPloidy));
    }

    // rescale first col by max
    gsl_vector_view firstCol = gsl_matrix_column(forwardMat, 0);
    scalingFactor = gsl_vector_max(&firstCol.vector);
    gsl_vector_scale(&firstCol.vector, 1.0 / scalingFactor);

    // iter through observed sequence
    gsl_vector_set_zero(this->prevForwardCol);
    gsl_vector_set_zero(this->currForwardCol);
    currResPloidy = 0;
    for(i = 1, depthIdx = 1; i < forwardMat->size2; i++, depthIdx++) {
      gsl_vector_set_zero(this->currForwardCol);
      gsl_matrix_get_col(this->prevForwardCol, forwardMat, i-1);

      // iter over states in col i
      for(stateIdx = 0; stateIdx < (int) this->states->size(); stateIdx++) {
        emissionProb = this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);
        currMax = -1;
        currMaxIdx = -1;

        // iter over possible transitions
        for(j = 0; j < this->prevForwardCol->size; j++) {
          currResPloidy = gsl_vector_get(this->currForwardCol, stateIdx) + exp(emissionProb) * gsl_vector_get(this->prevForwardCol, j) * gsl_matrix_get(this->transition, j, stateIdx);
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
    }

    // recover sequence
    // find best state at end of sequence (ie last col)
    gsl_vector_view lastCol = gsl_matrix_column(forwardMat, forwardMat->size2-1);
    currMax = gsl_vector_max(&lastCol.vector);
    currMaxIdx = gsl_vector_max_index(&lastCol.vector);

    // backtrack
    std::vector<std::forward_list<int>> paths(this->NUM_CELLS, std::forward_list<int>());
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currPloidy = getCellPloidyFromStateIdx(cellIdx, currMaxIdx);
      paths[cellIdx].push_front(currPloidy);
    }
    int prevMaxIdx = currMaxIdx;
    for(int k = (int) backTrace->size2-1; k > 0; k--) {
      currMaxIdx = gsl_matrix_get(backTrace, currMaxIdx, k);

      // count up number of state changes
      gsl_matrix_set(this->numTransitionsMat, currMaxIdx, prevMaxIdx, gsl_matrix_get(this->numTransitionsMat, currMaxIdx, prevMaxIdx) + 1); // going backwards, so indices are flipped
      prevMaxIdx = currMaxIdx;
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

    // for each chr position
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
  }
}

/*
 * function to find the steady state distribution of this HMM.
 * Assumes transtion matrix is already set, stores steady state dist
 * in passed steadStateVec
 */
double HMM::findSteadyStateDist(gsl_vector* steadyStateVec) const {
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
    std::cerr << "ERROR: could not find eigenvalues to find steady state distribution" << std::endl;
    // restore error handler
    gsl_set_error_handler(errHandler);
    return GSL_NAN;
  }

  // get eigenvalue == 1 and corresponding eigenvector. store in steadyStateVec
  bool foundSteadyState = false;
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
  }

  // restore error handler
  gsl_set_error_handler(errHandler);
  return GSL_SUCCESS;
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
 * These cutoffs work well in practice
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
  gsl_vector* countTransitionsOutOfPloidy = gsl_vector_alloc(this->MAX_PLOIDY + 1);
  gsl_vector* fracTransient = gsl_vector_alloc(this->MAX_PLOIDY + 1);

  this->viterbiDecode();

  // for each cell
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currDepths = (*this->depths)[cellIdx];
    gsl_vector_set_zero(countTransientTimes);
    gsl_vector_set_zero(countTimesInPloidy);
    gsl_vector_set_zero(countTransitionsOutOfPloidy);
    gsl_vector_set_zero(fracTransient);

    // for each chr
    for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
      std::string currChr = (*chrVec)[chrIdx];
      currVitPath = (*currDepths->chrToViterbiPathMap)[currChr];

      // for each window in currChr, get run length encoding
      pathSize = currVitPath->size();
      for(int regionIdx = 0; regionIdx < pathSize; regionIdx++) {
        stateCount = 1;
        currVitState = (*currVitPath)[regionIdx];

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
    }

    // check if this cell has any transient states. if so, set status to GSL_NAN and break. else, continue to next cell
    gsl_vector_memcpy(fracTransient, countTransientTimes);
    gsl_vector_div(fracTransient, countTransitionsOutOfPloidy);
    if(this->gradientDebug) {
      std::cout << "HMM::checkForTransientStates countTransientTimes, countTimesInPloidy, countTransitionsOutOfPloidy, fracTransient" << std::endl;
      printColVector(countTransientTimes);
      printColVector(countTimesInPloidy);
      printColVector(countTransitionsOutOfPloidy);
      printColVector(fracTransient);
    }
    double sumTransient = gsl_blas_dasum(countTransientTimes); // Double Absolute SUM
    if(sumTransient >= 4.9) {
      status = GSL_NAN;
      break;
    }
    if(gsl_isnan(status)) { // break out of cell loop
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
  for(unsigned int i = 0; i < this->states->size(); i++) {
    if(cummSum > p) {
      break;
    }
    cummSum += gsl_vector_get(probVec, i);
    toState = i;
  }
  return toState;
}

