#include "OneCell3TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
OneCell3TrParam2DegPolyHMM::OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy) : OneCell3TrParam2DegPolyHMM(depths, nullptr, maxPloidy, 2, 1, 0, 1) { // depths, nullptr fixed params, maxPloidy, 2 transition params to est (beta+gamma), 1 fixedTrParams (alpha), 0 numFixedLibs (est all 1 lib), 1 branch to est (t)
}
OneCell3TrParam2DegPolyHMM::OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : HMM (depths, fixedParams, numTrParamsToEst, 1, numBranches, maxPloidy, numFixedTrParams, numFixedLibs) {
  this->logFacKVec = nullptr;
  this->maxNumBFGSStarts = std::numeric_limits<int>::max();
  this->setAlpha(0.1); // arbitrary starting value

  /*// params to estimate with BFGS: the transition params and a library scaling factor for each cell
  // start libSizeScalingFactor at avgTumorDepth / avgDiploidDepth
  // [(t1 + t2 + ... + tY)/24] / [(d1 + d2 + ... + dY)/24] = (t1 + t2 + ... + tY) / (d1 + d2 + ... + dY)
  // removed Tue 14 Jan 2020 03:02:22 PM PST; call HMM's setLibScalingFactorsToAvg() instead
  double sumAvgTumorDepth = std::accumulate((*this->depths)[0]->avgTumorDepth->begin(), (*this->depths)[0]->avgTumorDepth->end(), 0.0); // TODO maybe not hardcode taking the 0th DepthPair?
  double sumAvgDiploidDepth = std::accumulate((*this->depths)[0]->avgDiploidDepth->begin(), (*this->depths)[0]->avgDiploidDepth->end(), 0.0);
  //gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX, sumAvgTumorDepth / sumAvgDiploidDepth);
  this->setLibScalingFactor(0, sumAvgTumorDepth / sumAvgDiploidDepth);*/
}

/*OneCell3TrParam2DegPolyHMM::OneCell3TrParam2DegPolyHMM(const OneCell3TrParam2DegPolyHMM& otherHMM) : HMM(otherHMM) {
  this->logFacKVec = nullptr;
  this->maxNumBFGSStarts = otherHMM.maxNumBFGSStarts;
}*/
OneCell3TrParam2DegPolyHMM::~OneCell3TrParam2DegPolyHMM() {
  delete this->logFacKVec;
}

//int OneCell3TrParam2DegPolyHMM::getMaxNumBFGSStarts() const {
//  return this->maxNumBFGSStarts;
//}

// ex hmm->print(stdout);
// ex hmm->print(stderr);
void OneCell3TrParam2DegPolyHMM::print(FILE* stream) {
  // print standard info first
  HMM::print(stream);

  // emission likelihood function
  fprintf(stream, "lambda_ij = ploidy_ij * %f * mu_i/2 + mu_i/%i\n", this->getLibScalingFactor(0), this->DEPTH_ERROR_SCALING);
  fprintf(stream, "X_ij ~ NegBinom(mean=lambda_ij, variance=%f + %f*lambda_ij + %f*lambda_ij^2)\n", this->getMeanVarianceIntercept(), this->getMeanVarianceSlope(), this->getMeanVariancePoly2());
}

void OneCell3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
}
double OneCell3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}
// TODO does it make sense to move this into HMM? this is copied from TwoCell3TrParam2DegPolyHMM
double OneCell3TrParam2DegPolyHMM::getAlpha() const {
  return gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX);
}
void OneCell3TrParam2DegPolyHMM::setAlpha(double alpha) {
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX, alpha);
}


/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
/*
 * //transitionParams = [beta, gamma, t]
 * transitionParams = [beta, lambda, t]
 */
double OneCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  //double beta  = gsl_vector_get(transitionParams, 0);
  //double gamma = gsl_vector_get(transitionParams, 1);
  //double t     = gsl_vector_get(transitionParams, 2);
  double beta   = gsl_vector_get(transitionParams, 0);
  double lambda = gsl_vector_get(transitionParams, 1);
  double t      = gsl_vector_get(transitionParams, 2);

  //std::cout << "alpha: " << alpha << ", beta: " << beta << ", lambda: " << lambda << ", t: " << t << std::endl;
  // save into paramsToEst if saving into this->transition
  if(dest == this->transition) {
    int transitionParamsIdx = 0;
    for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->TRANSITION_PROB_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
    for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
  }

  return this->setTransition(dest, this->getAlpha(), beta, lambda, t);
}

//double OneCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double gamma, double t) {
double OneCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t) {
  // first set rate matrix
  this->setRateMatrixQ(alpha, beta, lambda);
  //std::cout << "rate matrix q" << std::endl;
  //printMatrix(this->rateMatrixQ);

  // then set time dependent matrix P for time t
  double status = this->setTimeDepMatrixP(t);
  //std::cout << "OneCell3TrParam2DegPolyHMM::setTransition this->setTimeDepMatrixP(t) STATUS: " << status << std::endl;
  if(status != GSL_SUCCESS) {
    return status;
  }
  //std::cout << "time dep matrix P" << std::endl;
  //printMatrix(this->timeDepMatrixP);

  // then calculate matrix M(t) = {m_ij(t)}, where m_ij(t) = P_(2,2),(i,j)(t) / sum_v=0^k P_(2,2),(i,v)(t)
  int ancDiploidStateIdx = getStateIdxFromPloidyPair(2, 2); // index of ancestral diploid state (2,2)
  double currVal = 0;
  double rowsum = 0;
  for(unsigned int i = 0; i < this->states->size(); i++) {
    for(unsigned int j = 0; j < this->states->size(); j++) {
      currVal = gsl_matrix_get(this->timeDepMatrixP, ancDiploidStateIdx, getStateIdxFromPloidyPair(i, j));
      gsl_matrix_set(dest, i, j, currVal);
    }

    // rescale by rowsum
    gsl_vector_view currRow = gsl_matrix_row(dest, i);

    //std::cout << "unscaled row " << i << " in transition mat" << std::endl;
    //printRowVector(&currRow.vector);

    rowsum = gsl_blas_dasum(&currRow.vector); // Double Absolute SUM
    gsl_vector_scale(&currRow.vector, 1.0 / rowsum);

    //std::cout << "scaled row " << i << " in transition mat" << std::endl;
    //printRowVector(&currRow.vector);
  }
  //std::cout << "scaled dest" << std::endl;
  //printMatrix(dest);

  //std::cout << "OneCell3TrParam2DegPolyHMM::setTransition return STATUS: " << status << std::endl;
  return status;

  //// alpha: P(one CNA)
  //// beta: P(any CNA)
  //// gamma: P(return to diploid)
  //// t: branch length (scale everything by this)
  ////      ploidy x ploidy
  ////     0         1          2                      ... k
  //// 0 [ 1-rowSum  (a+b)*t    (b+g)*t                ... b*t]
  //// 1 [ (a+b)*t   1-rowSum   (a+b+g)*t  b*t     b*t ... b*t]
  //// 2 [ b*t       (a+b)*t    1-rowSum  (a+b)*t  b*t ... b*t]
  //// 3 [ ...                                             b*t]
  //// . [ ...                                             b*t]
  //// k [ b*t                                    ... 1-rowSum]

  //// init all values to beta*t
  //gsl_matrix_set_all(dest, beta*t);

  //// if haven't saved the string represetation of the transition matrix yet, save it now
  //bool saveStrs = false;
  //std::string* currTrStr = nullptr;
  //if(this->transitionStrings == nullptr) {
  //  saveStrs = true;
  //  this->transitionStrings = new std::vector<std::vector<std::string*>*>(this->states->size());
  //  for(unsigned int i = 0; i < this->states->size(); i++) {
  //    (*this->transitionStrings)[i] = new std::vector<std::string*>(this->states->size());
  //    for(unsigned int j = 0; j < this->states->size(); j++) {
  //      currTrStr = new std::string();
  //      *currTrStr = "b*t";
  //      (*(*this->transitionStrings)[i])[j] = currTrStr;
  //    }
  //  }
  //}

  //// set subdiagonal and superdiagonal elements to (a+b)*t
  //gsl_vector_view subdiag = gsl_matrix_subdiagonal(dest, 1);
  //gsl_vector_set_all(&subdiag.vector, (alpha+beta)*t);

  //gsl_vector_view superdiag = gsl_matrix_superdiagonal(dest, 1);
  //gsl_vector_set_all(&superdiag.vector, (alpha+beta)*t);
  //if(saveStrs) {
  //  for(unsigned int row = 0; row < this->states->size(); row++) {
  //    for(unsigned int col = 0; col < this->states->size(); col++) {
  //      if(col == row + 1 || col == row - 1) {
  //        currTrStr = (*(*this->transitionStrings)[row])[col];
  //        *currTrStr = "(a+b)*t";
  //      }
  //    }
  //  }
  //}

  //// set diploid col to (b+g)*t (assumes there are enough rows and columns, no safety checks here)
  //gsl_vector_view diploidCol = gsl_matrix_column(dest, 2);
  //for(unsigned int row = 0; row < this->states->size(); row++) {
  //  // set super and sub diag in diploid col special case to be (a+b+g)*t
  //  if(row == 1 || row == 3) {
  //    gsl_vector_set(&diploidCol.vector, row, (alpha+beta+gamma)*t);
  //    if(saveStrs) {
  //      currTrStr = (*(*this->transitionStrings)[row])[2];
  //      *currTrStr = "(a+b+g)*t";
  //    }
  //  }
  //  else {
  //    gsl_vector_set(&diploidCol.vector, row, (beta+gamma)*t);
  //    if(saveStrs) {
  //      currTrStr = (*(*this->transitionStrings)[row])[2];
  //      *currTrStr = "(b+g)*t";
  //    }
  //  }
  //}

  //// reset diagonal elements to 1-rowSum
  //gsl_vector_view diag = gsl_matrix_diagonal(dest);
  //gsl_vector_set_all(&diag.vector, 0); // first clear everything out
  //for(unsigned int row = 0; row < this->states->size(); row++) {
  //  gsl_vector_view currRow = gsl_matrix_row(dest, row);
  //  double rowSum = gsl_blas_dasum(&currRow.vector); // Double Absolute SUM
  //  gsl_matrix_set(dest, row, row, 1.0 - rowSum);
  //  if(saveStrs) {
  //    currTrStr = (*(*this->transitionStrings)[row])[row];
  //    *currTrStr = "(r)";
  //  }
  //}

  //// check if anything went negative
  //// if somehow the transition matrix went negative, return NAN
  //double transitionMatMin = gsl_matrix_min(dest);
  //if(transitionMatMin < 0 || gsl_isnan(transitionMatMin)) {
  //  return GSL_NAN;
  //}

  return 0;
}

/*
 * Because BFGS optimizes unrestrained values, we need to transform
 * probabilities and parameters at different steps.
 *
 *   BFGS will optimize params x,y \in (-inf, inf)
 *   Probs a,b must be \in [0,1]
 *   and for the transition matrix, the constraint is
 *   //0 <= 1-2a-kb-g <= 1 <==> 2a+kb+g <= 1
 *   0 <= 1-t(2a+kb+g) <= 1 <==> 0 <= 2a+kb+g <= 1 && 0 <= t <= 1
 *
 * The following transformation converts params x,y,z to probabilities a,b,g
 *   2a = e^x / (1 + e^x + e^y + e^z) ==> a = e^x / [2 (1 + e^x + e^y + e^z)]
 *   kb = e^y / (1 + e^x + e^y + e^z) ==> b = e^y / [k (1 + e^x + e^y + e^z)]
 *   g = e^z / (1 + e^x + e^y + e^z) ==> g = e^z / [(1 + e^x + e^y + e^z)]
 * 
 * The following transformation converts probs a,b,g to params x,y,z
 *   Let c = 1+ e^x + e^y + e^z
 *   2a + kb + g + (1-2a-kb-g) = 1 ==> e^x / c + e^y / c + e^z / c + 1 / c = 1
 *     ==> (1-2a-kb-g) = 1 / c
 *     ==> c = 1 / (1-2a-kb-g)
 *   a = e^x / (2c) ==> a * 2c = e^x
 *   x = ln(a * 2c)
 *   y = ln(b * kc)
 *   z = ln(g * c)
 *
 * The following transformation converts branch length t (copied/reduced from TwoCell3TrParam2DegPolyHMM)
 *   d*t = e^(T) / (1 + exp(T))
 *   T = ln(-(d*t) / (d*t - 1))
 *
 * The library scaling factor (s) must be strictly positive. BFGS will optimize w \in (-inf, inf)
 *   s = exp(w) > 0
 *   w = ln(s)
 *
 * The following two functions convert the elements of src into the elements of dest
 * Also, copies over any other elements into dest (ie libSizeScalingFactor)
 */
void OneCell3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  // Wed 09 Sep 2020 04:09:01 PM PDT all vars must be > 0, but that's it. just exp/log everything for continuous time model
  //double d = (double) (*this->depths)[0]->maxWindowSize;
  double beta = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  //double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX) / d;
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);

  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);

  //double a = this->getAlpha();
  //double c = (b * this->getKploidy() + g - 1);

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z

  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(-(d * t) / (d * t - 1))); // set T

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z

  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t)); // set T

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(s));
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCell3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);

  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);

  //double a = this->getAlpha();
  //double c = 1 - 2.0*a + exp(y) + exp(z);

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // set beta
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z) / c); // set gamma
  //c = 1.0 / (1 + exp(T));

  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T) * c); // set t

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y)); // set beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z)); // set lambda
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T)); // set t

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
}

/*
 ********
 * functions
 ********
 */
// Thu 23 Jul 2020 05:05:33 PM PDT copied over from TwoCell3TrParam2DegPolyHMM
double OneCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  // if diploidDepth == 0, assume this is a hard to map region (ie any tumor reads mapping here are spurious)
  // P(any number reads | no data) = 1
  if(diploidDepth < 1e-4) {
    return 1;
  }

  double currLibSizeScalingFactor = this->getLibScalingFactor(cellIdx);
  // Fri 07 Feb 2020 09:50:37 PM PST trying out new X_ij ~ NB(r_ik + epsilon)
  //double lambda_i = (ploidy * diploidDepth/2.0 + 272.5568 / this->DEPTH_ERROR_SCALING) * currLibSizeScalingFactor; // withErr, unnormDiploid
  //double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor; // no error
  double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor + 272.5568 / this->DEPTH_ERROR_SCALING; // constant error
  //double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor + 1; // constant error of 1

  // tumorDepth ~ NegBinom(lambda_i, var=intercept + slope * lambda_i + poly * lambda_i + poly2 * lambda_i * lambda_i
  double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i;

  // in the boost library, p = P(success), r = # successes
  double p = lambda_i / var;
  double r = lambda_i * lambda_i / (var - lambda_i);

  //std::cerr << "currLibSizeScalingFactor: " << currLibSizeScalingFactor << ", lambda_i: " << lambda_i << ", var: " << var << ", p: " << p << ", r: " << r << std::endl;
  // if diploidDepth is so low that p or r goes negative, assume this is a hard to map region
  if((p < 0 || r < 0) && diploidDepth < 1) {
    //std::cerr << "P OR R TOO SMALL AND DIPLOIDDEPTH < 1; RETURNING 1" << std::endl;
    return 1;
  }

  // Mon 30 Mar 2020 03:55:48 PM PDT try matching the simulation fix for small r
  //std::cerr << r << std::endl;
  if(r < 1) {
    //std::cerr << "R TOO SMALL; ADJUSTING" << std::endl;
    r = 1;
    var = lambda_i * lambda_i + lambda_i;
    p = lambda_i / var;
  }

  double likelihood = -1;
  int k = (int) tumorDepth;
  likelihood = exp(lgamma(r+k) - lgamma(r) - this->getLogFacK(k) + r * log(p) + k * log(1-p));
  return likelihood;
}

// copied from TwoCell3TrParam2DegPolyHMM
double OneCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  double emissionProb = 0;
  int currPloidy = -1;
  double currTumorDepth = -1;
  double currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
  // for each cell, calc the emission prob and multiply by it
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
    currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
    emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx));
  }
  return emissionProb;
}

/*
 * helper function to return log(k!) part of the negative binomial coefficient
 * given k (some read depth).
 * to be consistent with prior boost implementation:
 * f(k; r, p) = gamma(r+k)/(k! gamma(r)) * p^r * (1-p)^k
 * https://www.boost.org/doc/libs/1_70_0/libs/math/doc/html/math_toolkit/dist_ref/dists/negative_binomial_dist.html
 * The underlying member variable vector (this->logFacKVec)
 * has structure k:[log(k!) == lgamma(k+1)]
 */
double OneCell3TrParam2DegPolyHMM::getLogFacK(int k) {
  // if haven't called getLogFacK yet, alloc the lookup vector
  if(this->logFacKVec == nullptr) {
    this->setLogFacK();
  }
  return (*this->logFacKVec)[k];
}
/*
 * helper method to set this->logFacKVec
 */
void OneCell3TrParam2DegPolyHMM::setLogFacK() {
  // k is 0..maxReadDepth
  int maxDepth = *(this->alphabet->rbegin()); // sets are stored in sorted order
  if(this->logFacKVec == nullptr) {
    this->logFacKVec = new std::vector<double>();
  }
  for(int i = this->logFacKVec->size(); i < maxDepth+1; i++) {
    this->logFacKVec->push_back(lgamma(i+1));
  }
}

/*
 * calls bfgs, returns a bestGuessHMM
 */
OneCell3TrParam2DegPolyHMM* OneCell3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, bool verbose) {
  // create new HMM with the best guess parameters and return it
  //OneCell3TrParam2DegPolyHMM* bestGuessHMM = new OneCell3TrParam2DegPolyHMM(*this);
  OneCell3TrParam2DegPolyHMM* bestGuessHMM = this;//new OneCell3TrParam2DegPolyHMM(*this);
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  //this->convertProbToParam(initGuessAsParams, initGuess, this->getKploidy());
  this->convertProbToParam(initGuessAsParams, initGuess);
  //HMM::bfgs(initGuess, bestGuessHMM, verbose);
  //HMM::bfgs(initGuessAsParams, bestGuessHMM, verbose);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, verbose);
  return bestGuessHMM;
}
  
/*
 * method to simulate a sequence of states, then a sequence of read depths
 * given those states.
 * Tue 11 Jun 2019 04:03:05 PM PDT seed is only used once (first time simulate is called)
 */
void OneCell3TrParam2DegPolyHMM::simulate() {
  this->simulate(43);
}
void OneCell3TrParam2DegPolyHMM::simulate(int seed) {
  this->simulate(seed, true); // by default, simulate diploid depths
}
void OneCell3TrParam2DegPolyHMM::simulate(int seed, bool simDiploid, double diploid_lambda_i, int numDiploidCells) {
  if(this->generator == nullptr) {
    this->generator = new base_generator_type(seed);
  }
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->generator, uni_dist);

  // simulate states
  std::vector<std::string>* diploidStates = nullptr;
  std::vector<std::string>* tumorStates = nullptr;
  gsl_matrix_set_zero(this->numTransitionsMat);

  // clear out DepthPair variables. This is important if multiple simulations are run
  for(unsigned int i = 0; i < this->depths->size(); i++) {
    if(simDiploid) {
      (*this->depths)[i]->diploidLibrarySize = 0;
      (*this->depths)[i]->maxDiploidDepth = 0;
    }
    (*this->depths)[i]->tumorLibrarySize = 0;
    (*this->depths)[i]->maxTumorDepth = 0;
  }

  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  int fromStateIdx = -1;
  int toStateIdx = -1;
  // for each chr
  for(unsigned int i = 0; i < chrVec->size(); i++) {
    // init this chr's vectors
    diploidStates = new std::vector<std::string>();
    tumorStates = new std::vector<std::string>();

    currChr = (*chrVec)[i];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();

    // diploid always in state diploid
    diploidStates->push_back("2");

    // start tumor in steady state
    double stateProb = uni();
    fromStateIdx = getRandStateIdx(stateProb, this->initProb);
    tumorStates->push_back((*this->states)[fromStateIdx]);

    // for each window in this chr
    for(int winIdx = 1; winIdx < currNumWindows; winIdx++) {
      // diploid always in state diploid
      diploidStates->push_back("2");

      // tumor state changes according to transition matrix
      stateProb = uni();
      toStateIdx = getRandStateIdx(stateProb, fromStateIdx);
      tumorStates->push_back((*this->states)[toStateIdx]);
      gsl_matrix_set(this->numTransitionsMat, fromStateIdx, toStateIdx, gsl_matrix_get(this->numTransitionsMat, fromStateIdx, toStateIdx) + 1);
      fromStateIdx = toStateIdx;
    }

    // save
    (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr] = diploidStates;
    (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr] = tumorStates;
  }

  double diploid_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * diploid_lambda_i + this->getMeanVariancePoly2() * diploid_lambda_i * diploid_lambda_i; // Mon 20 Apr 2020 08:10:12 PM PDT now all diploid variance is calculated using the polynomial mean/var relationship
  double diploid_p = diploid_lambda_i / diploid_var;
  double diploid_r = diploid_lambda_i * diploid_lambda_i / (diploid_var - diploid_lambda_i);
  //std::cout << diploid_r << ", " << diploid_p  << ", " << std::endl;
  boost::random::negative_binomial_distribution<> diploid_negBinom_dist(diploid_r, diploid_p);
  boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->generator, diploid_negBinom_dist);

  // given states, simulate coverage according to neg binom model
  double tumor_lambda_i = -1;
  double tumor_var = -1;
  double tumor_p = -1;
  double tumor_r = -1;
  int currTumorPloidy = -1;
  double maxSim = -1;
  double currDiploidSim = -1;
  double currTumorSim = -1;
  double currLibSizeScalingFactor = -1;

  if(simDiploid) {
    // for each chr
    for(unsigned int j = 0; j < chrVec->size(); j++) {
      currChr = (*chrVec)[j];
      diploidStates = (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr];
      for(unsigned int i = 0; i < diploidStates->size(); i++) {
        double dipSimVec[numDiploidCells];
        for(int dipIdx = 0; dipIdx < numDiploidCells; dipIdx++) {
          dipSimVec[dipIdx] = diploid_negBinom();
          std::cerr << dipSimVec[dipIdx] << std::endl;
        }
        currDiploidSim = gsl_stats_mean(dipSimVec, 1, numDiploidCells);
        (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i] = currDiploidSim;
        if(currDiploidSim > maxSim) {
          maxSim = currDiploidSim;
          (*this->depths)[0]->maxDiploidDepth = currDiploidSim;
        }
        (*this->depths)[0]->diploidLibrarySize += currDiploidSim;
        (*this->depths)[1]->diploidLibrarySize += currDiploidSim;
        //(*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = diploid_var;
        (*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = gsl_stats_variance(dipSimVec, 1, numDiploidCells); // Fri 17 Apr 2020 10:42:38 PM PDT debugging simulating multiple diploid cells, then averaging
      }
    }
  }

  // sim tumor
  // for each chr
  for(unsigned int j = 0; j < chrVec->size(); j++) {
    currChr = (*chrVec)[j];
    tumorStates = (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr];

    // for each simulated state
    for(unsigned int i = 0; i < tumorStates->size(); i++) {
      diploid_lambda_i = (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i];

      // parse tumor states
      currTumorPloidy = atoi((*tumorStates)[i].c_str());
      currLibSizeScalingFactor = this->getLibScalingFactor(0);
      tumor_lambda_i = currLibSizeScalingFactor * (currTumorPloidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      if(tumor_r < 1) {
        // fix r=1, adj var, recalc p
        // r = lambda^2 / (var - lambda); v = lambda^2 / r + lambda
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;
      }

      boost::random::negative_binomial_distribution<> tumor_negBinom_dist(tumor_r,tumor_p);
      boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor_negBinom(*this->generator, tumor_negBinom_dist);
      currTumorSim = tumor_negBinom();

      // save simulated values, as well as max for alphabet
      (*(*(*this->depths)[0]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
      if(currTumorSim > maxSim) {
        maxSim = currTumorSim;
        (*this->depths)[0]->maxTumorDepth = currTumorSim;
      }
      (*this->depths)[0]->tumorLibrarySize += currTumorSim;
    }
  }

  // save alphabet
  std::vector<int> fullAlphabet((int)maxSim + 1); // add one to make inclusive
  std::iota(fullAlphabet.begin(), fullAlphabet.end(), 0);
  this->setAlphabet(new std::set<int>(fullAlphabet.begin(), fullAlphabet.end()));

  // save this->logFacKVec
  this->setLogFacK();

  // set library sizes according to final tumor size
  this->estLibScalingFactorsPosterior();
}

void OneCell3TrParam2DegPolyHMM::setUpBaumWelchLeastSquares() {
  this->baumWelchTransitionMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
//  /*
//   * transition matrix entries are of the form:
//   * (a+b+g)*t
//   * (a+b)*t
//   * (b+g)*t
//   * b*t
//   * which can be rewritten as f(.) = (xa + yb + zg) * ut
//   * where u,x,y,z are all indicator variables (and u is always 1)
//   *
//   * We want to solve the overdetermined system by BFGS on least squares
//   * S = sum_i r_i^2
//   * r_i = a_ij* - f(.)
//   */
//
//  // one residual per coef entry, which has one entry per transition matrix element
//  this->residuals = gsl_vector_alloc(this->states->size() * this->states->size());
//  gsl_vector_set_zero(this->residuals);
//
//  /*
//   * structure of coefs:
//   *   x  y  z  u
//   * [            ] // a*[0,0]
//   * [            ] // a*[0,1]
//   *      ...
//   * [            ] // a*[k,k]
//   */
//  int xIdx = 0; // alpha
//  int yIdx = 1; // beta
//  int zIdx = 2; // gamma
//  int uIdx = 3; // t
//  this->coefs = gsl_matrix_alloc(this->states->size() * this->states->size(), 1 + this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST); // +1 for alpha coef
//  gsl_matrix_set_zero(coefs);
//
//  // fill in coefs. logic adapted from TwoCell3TrParam2DegPolyHMM::setTransition
//  int coefsIdx = 0;
//  gsl_vector* zerosVec = gsl_vector_alloc(this->coefs->size2);
//  gsl_vector_set_zero(zerosVec);
//  gsl_vector* onesVec = gsl_vector_alloc(this->coefs->size1);
//  gsl_vector_set_all(onesVec, 1.0);
//  gsl_matrix_set_col(coefs, yIdx, onesVec); // beta appears everywhere
//  gsl_matrix_set_col(coefs, uIdx, onesVec); // t appears everywhere
//  gsl_vector_free(onesVec);
//  for(int from = 0; from < (int) this->states->size(); from++) {
//    for(int to = 0; to < (int) this->states->size(); to++, coefsIdx++) {
//
//      // if from == to (no move)
//      if(from == to) {
//        // diagonal; reset row to 0
//        gsl_matrix_set_row(coefs, coefsIdx, zerosVec);
//      }
//      else {
//        // if abs(from - to) == 1
//        if(std::abs(from - to) == 1) { // adj CNA (alpha)
//          gsl_matrix_set(coefs, coefsIdx, xIdx, 1);
//        }
//        // if to == 2 (return to diploid (gamma))
//        if(to == 2) {
//          gsl_matrix_set(coefs, coefsIdx, zIdx, 1);
//        }
//      }
//    }
//  }
//  gsl_vector_free(zerosVec);
//
//  this->a_ij = gsl_vector_alloc(this->states->size() * this->states->size());
//  int a_ijIdx = 0;
//  for(unsigned int row = 0; row < this->states->size(); row++) {
//    for(unsigned int col = 0; col < this->states->size(); col++, a_ijIdx++) {
//      gsl_vector_set(a_ij, a_ijIdx, gsl_matrix_get(this->transition, row, col));
//    }
//  }
}

/*
 * param probs should be in prob space (ie already converted away from BFGS space). Order should be
 * [b]
 * [g]
 * [t]
 */
double OneCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs) {
  double sumSqResid = 0;
  double resid = 0;

  // recalc what the transition matrix would be under these probs
  double status = this->setTransition(this->baumWelchTransitionMat, probs);
  if(gsl_isnan(status)) {
    return status;
  }

  // compare each entry in the hypothetical proposed transition matrix to the baum welch calculated transition matrix (ie sum up the squared difference)
  for(unsigned int row = 0; row < this->transition->size1; row++) {
    for(unsigned int col = 0; col < this->transition->size2; col++) {
      resid = (gsl_matrix_get(this->baumWelchTransitionMat, row, col) - gsl_matrix_get(this->transition, row, col));
      sumSqResid += resid * resid;
    }
  }
  //std::cout << "OneCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f: " << sumSqResid << std::endl;
  //printRowVector(probs);
  //printMatrix(this->transition);
  //printMatrix(this->baumWelchTransitionMat);
  return sumSqResid;


  //double sumSqResid = 0;
  //double resid = 0;
  //double transitionProbs = 0;
  //double branchLengths = 0;
  //double coefVarProd = 0;
  //for(unsigned int residIdx = 0; residIdx < residuals->size; residIdx++) {
  //  // skip over the diagonals
  //  if(residIdx % (this->states->size() + 1) == 0) {
  //    resid = 0;
  //  }
  //  else {
  //    transitionProbs  = gsl_matrix_get(coefs, residIdx, 0) * this->getAlpha(); // alpha
  //    transitionProbs += gsl_matrix_get(coefs, residIdx, 1) * gsl_vector_get(probs, 0); // beta
  //    transitionProbs += gsl_matrix_get(coefs, residIdx, 2) * gsl_vector_get(probs, 1); // gamma
  //    branchLengths  = gsl_matrix_get(coefs, residIdx, 3) * gsl_vector_get(probs, 2); // t

  //    coefVarProd = transitionProbs * branchLengths;
  //    resid = gsl_vector_get(a_ij, residIdx) - coefVarProd;
  //  }
  //  gsl_vector_set(residuals, residIdx, resid);
  //  sumSqResid += resid * resid;
  //}
  //return sumSqResid;
}

void OneCell3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  // TODO
}

/*
 * method to check validity of parameters proposed by BFGS. Assumes probs
 * is in probability space.
 * returns 0 if probs is valid, GSL_NAN otherwise
 */
double OneCell3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
      //std::cerr << "shortcut for bad lib sizes: ";// << std::endl;
      //printRowVector(stderr, probs);
      return GSL_NAN;
    }
  }

  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
    //std::cerr << "shortcut for any probabilities becoming too small: " ;//<< probMin << std::endl;
    //printRowVector(stderr, probs);
    return GSL_NAN;
  }
  return 0;
}

