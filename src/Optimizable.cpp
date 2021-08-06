#include "Optimizable.hpp"

Optimizable::Optimizable() {
  this->paramsToEst = nullptr;
  this->fixedParams = nullptr;
  this->optimSuccess = false;
  this->generator = nullptr;
  this->initGuessCopyBeforeBFGS = nullptr;
  this->bestOptimLLInitGuess = nullptr;
  this->BFGSParamResults = nullptr;
  this->maxNumBFGSStarts = std::numeric_limits<int>::max();
  this->changeInBFGSLoglikelihood = 0;
  this->probParamConversionVec = nullptr;
  this->gradientDebug = false;
  this->simParamsToEst = nullptr;
  this->simFixedParams = nullptr;
}
Optimizable::Optimizable(const Optimizable& other) {
  this->paramsToEst = gsl_vector_alloc(other.getNumParamsToEst());
  gsl_vector_memcpy(this->paramsToEst, other.paramsToEst);
  this->fixedParams = gsl_vector_alloc(other.getNumFixedParams());
  gsl_vector_memcpy(this->fixedParams, other.fixedParams);
  this->optimSuccess = false;
  this->generator = nullptr;
  this->initGuessCopyBeforeBFGS = nullptr;
  this->bestOptimLLInitGuess = nullptr;
  this->BFGSParamResults = nullptr;
  this->maxNumBFGSStarts = other.maxNumBFGSStarts;
  this->changeInBFGSLoglikelihood = 0;
  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst());
  this->gradientDebug = other.gradientDebug;
  if(other.simParamsToEst != nullptr) {
    this->simParamsToEst = gsl_vector_alloc(other.simParamsToEst->size);
    gsl_vector_memcpy(this->simParamsToEst, other.simParamsToEst);
  }
  if(other.simFixedParams != nullptr) {
    this->simFixedParams = gsl_vector_alloc(other.simFixedParams->size);
    gsl_vector_memcpy(this->simFixedParams, other.simFixedParams);
  }
}

gsl_vector* Optimizable::getParamsToEst() const {
  return this->paramsToEst;
}
gsl_vector* Optimizable::getFixedParams() const {
  return this->fixedParams;
}
int Optimizable::getNumParamsToEst() const {
  return this->paramsToEst->size;
}
int Optimizable::getNumFixedParams() const {
  return this->fixedParams->size;
}
gsl_vector* Optimizable::getBestLLInitGuess() const {
  return this->bestOptimLLInitGuess;
}
std::vector<gsl_vector*>* Optimizable::getBFGSParamResults() const {
  return this->BFGSParamResults;
}
int Optimizable::getMaxNumBFGSStarts() const {
  return this->maxNumBFGSStarts;
}
/*
 * evalLikelihoodAtPoint function assumes
 * v is the current point to evaluate at in BFGS space, with the same format as paramsToEst
 */
double Optimizable::evalLikelihoodAtPoint(const gsl_vector* v, void* params) {
  Optimizable* optimObj = (Optimizable*) params;

  // convert vector v (bfgs params point to evaluate) back to probs
  gsl_vector* probs = optimObj->probParamConversionVec;
  optimObj->convertParamToProb(probs, v);

  // if can't find steady state dist using this param set, return NAN
  double status = optimObj->checkOptimProbValidity(probs);
  if(gsl_isnan(status)) {
    return status;
  }
  status = optimObj->setParamsToEst(probs); // includes a call to this->setTransition()

  if(gsl_isnan(status)) {
    return status;
  }

  // call forward alg
  double loglikelihood = -optimObj->getLogLikelihood(); // negate because gsl provides a minimizer, and we want to maximize the likelihood
  return loglikelihood;
}

/* 
 * evalGradientAtPoint evaluates the gradient at point v using the finite difference method, and stores
 * the gradient in df
 */
void Optimizable::evalGradientAtPoint(const gsl_vector* v, void* params, gsl_vector* df) {
  // approx using finite difference method:
  // f'(x) = [f(x + h) - f(x - h)] / (2h)
  double h = 1e-4; // step size
  double f_x = evalLikelihoodAtPoint(v, params);
  if(gsl_isnan(f_x)) {
    gsl_vector_set_all(df, GSL_NAN);

    gsl_vector* probVec = gsl_vector_alloc(v->size); // added for printGradientPerIter
    ((Optimizable*)params)->convertParamToProb(probVec, v);
    if(((Optimizable*)params)->gradientDebug) {
      fprintf(stderr, "%.20f\t", f_x);
      fprintf(stderr, "|\t");
      for(unsigned int j = 0; j < probVec->size; j++) {
        fprintf(stderr, "%.20f\t", gsl_vector_get(probVec, j));
      }
      fprintf(stderr, "|\t");
      for(unsigned int j = 0; j < v->size; j++) {
        fprintf(stderr, "%.20f\t", gsl_vector_get(v, j));
      }
      fprintf(stderr, "|\t");
      for(unsigned int j = 0; j < df->size; j++) {
        fprintf(stderr, "%.20f\t", gsl_vector_get(df, j));
      }
      fprintf(stderr, "\n");
    }

    return;
  }
  double f_x_ph = 0; // f(x+h)
  double f_x_mh = 0; // f(x-h)
  double x = 0;
  double derivApprox = 0;
  gsl_vector_set_zero(df);
  gsl_vector* currVec = gsl_vector_alloc(v->size);
  gsl_vector_set_zero(currVec);
  gsl_vector* probVec = gsl_vector_alloc(v->size); // added for printGradientPerIter
  gsl_vector_set_zero(probVec); // added for printGradientPerIter
  for(unsigned int i = 0; i < v->size; i++) {
    gsl_vector_memcpy(currVec, v);
    x = gsl_vector_get(v, i);
    gsl_vector_set(currVec, i, x + h);
    f_x_ph = evalLikelihoodAtPoint(currVec, params);

    if(gsl_isnan(f_x_ph)) {
      derivApprox = GSL_NAN;
    } else {
      gsl_vector_set(currVec, i, x - h);
      f_x_mh = evalLikelihoodAtPoint(currVec, params);
      derivApprox = (f_x_ph - f_x_mh) / (2*h); // central difference
    }
    gsl_vector_set(df, i, derivApprox);

  }

    // added for a more granular version of printGradientPerIter
  if(((Optimizable*)params)->gradientDebug) {
    ((Optimizable*)params)->convertParamToProb(probVec, currVec);
    fprintf(stderr, "-%.20f\t", f_x);
    fprintf(stderr, "|\t");
    for(unsigned int j = 0; j < probVec->size; j++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(probVec, j));
    }
    fprintf(stderr, "|\t");
    for(unsigned int j = 0; j < currVec->size; j++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(currVec, j));
    }
    fprintf(stderr, "|\t");
    for(unsigned int j = 0; j < df->size; j++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(df, j));
    }
    fprintf(stderr, "\n");
  }
    


  gsl_vector_free(currVec);
  gsl_vector_free(probVec); // added for printGradientPerIter
}
void Optimizable::evalLikelihoodGradAtPoint(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
  *f = evalLikelihoodAtPoint(v, params); 
  evalGradientAtPoint(v, params, df);
}

/*
 * function to estimate parameters using bfgs and forward algorithm.
 * bestGuessOptim is the dest HMM to copy into; should be specific to the subclass.
 * initGuess is in BFGS space
 */
void Optimizable::bfgs(gsl_vector* initGuess, Optimizable* bestGuessOptim, bool verbose) const {
  // set up GSL optimizer
  gsl_multimin_function_fdf my_func; // which function to minimize
  my_func.df = &evalGradientAtPoint; // gradient of function
  my_func.fdf = &evalLikelihoodGradAtPoint; // how to set both f and df
  my_func.f = &evalLikelihoodAtPoint; // function itself
  my_func.params = bestGuessOptim; // this optimizable should be passed around

  // starting point/initial guess, in BFGS space
  gsl_vector* x = nullptr;
  gsl_vector* probs = nullptr;
  int status = GSL_CONTINUE;
  double initLikelihood = 0;

  if(initGuess != nullptr) {
    // before doing anything, check if probs are valid
    probs = gsl_vector_alloc(initGuess->size);
    bestGuessOptim->convertParamToProb(probs, initGuess);

    status = bestGuessOptim->checkOptimProbValidity(probs);
    if(isnan(status)) {
      std::cerr << "ERROR: initial param set is invalid, not running optimization" << std::endl;
      std::cerr << "##################################################################" << std::endl;
      bestGuessOptim->optimSuccess = false;
      return;
    }

    x = initGuess;
    my_func.n = initGuess->size;

    initLikelihood = Optimizable::evalLikelihoodAtPoint(x, bestGuessOptim);
    bestGuessOptim->initLl = initLikelihood;
    if(verbose) {
      printf("BFGS INITIAL LIKELIHOOD: %.40f\n", -initLikelihood);
    }
    if(isnan(initLikelihood)) {
      std::cerr << "ERROR: initial param set yields nan loglikelihood, not running optimization" << std::endl;
      std::cerr << "##################################################################" << std::endl;
      bestGuessOptim->optimSuccess = false;
      return;
    }
    else if(gsl_isnan(this->checkStateValidity(1e-4))) {
      std::cout << "ERROR: at least one HMM had an invalid transition matrix from this initial BFGS parameter set, not running optimization" << std::endl;
      std::cerr << "##################################################################" << std::endl;
      bestGuessOptim->optimSuccess = false;
      return;
    }
  }
  else {
    my_func.n = bestGuessOptim->getNumParamsToEst();
    x = gsl_vector_alloc(my_func.n);
    probs = gsl_vector_alloc(my_func.n);
    convertProbToParam(x, bestGuessOptim->paramsToEst);
  }

  const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2; // which algorithm to use
  gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, my_func.n);
  gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-4, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)

  int iter = 0;
  double oldLik = initLikelihood;
  double currLik = 0;
  double deltaLik = 0;
  int countTooClose = 0;
  gsl_vector* prevBestParams = gsl_vector_alloc(my_func.n); // store the previously best params in case currLik turns into NaN
  gsl_vector_memcpy(prevBestParams, x);
  status = GSL_CONTINUE; // reset status flag

  // timing from https://stackoverflow.com/a/27739925
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;
  double elapsedSec = 0;
  double totalTime = 0;

  while(status == GSL_CONTINUE) {
    begin = std::chrono::steady_clock::now();
    status = gsl_multimin_fdfminimizer_iterate(s); // do one iter
    end = std::chrono::steady_clock::now();

    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    iter++;
    totalTime += elapsedSec;
    currLik = gsl_multimin_fdfminimizer_minimum(s);
    deltaLik = std::abs(currLik - oldLik);

    // print update message
    if(verbose) {
      printf("ON BFGS ITER %d, likelihood %.20f, time elapsed (sec) %.5f, iter change in loglikelihood %.5f, total change in loglikelihood %.5f\n", iter, -currLik, elapsedSec, deltaLik, initLikelihood - currLik);
    }

    // if any error has occurred, break
    if(status) {
      std::cout << "STATUS IS: " << gsl_strerror(status) << std::endl;
      break;
    }

    // if currLik is NaN, stop
    if(gsl_isnan(currLik)) {
      printf("STATUS IS: ERROR: current likelihood is %f\n", currLik);
      break;
    }

    // check if have converged (is the gradient close enough to 0?)
    status = gsl_multimin_test_gradient(s->gradient, 1e-4);
    if(status == GSL_SUCCESS) {
      printf("STATUS IS: SUCCESS: Minimum found.\n");
    }

    // check if done too many iters
    if(iter >= 500) {
      std::cout << "STATUS IS: minimum not found in 500 iterations" << std::endl;
      break;
    }

    // check if have converged based on change in likelihood
    if(deltaLik < 1e-6) {
      countTooClose++;
      if(countTooClose >= 3) {
        std::cout << "STATUS IS: converged by consecutive small change in likelihood" << std::endl;
        break;
      }
    } else {
      countTooClose = 0;
    }
    oldLik = currLik;
    gsl_vector_memcpy(prevBestParams, gsl_multimin_fdfminimizer_x(s));
  }

  if(bestGuessOptim->gradientDebug) {
    std::cerr << "##################################################################" << std::endl;
  }

  // save results
  if(gsl_isnan(currLik)) {
    bestGuessOptim->optimSuccess = false;

    // use the previously best results
    if(verbose) {
      std::cout << "Using previously best params" << std::endl;
      std::cerr << "ERROR: Optimzation failed: likelihood is nan, using previously best params" << std::endl;
    }
    bestGuessOptim->convertParamToProb(probs, prevBestParams);
  }
  // check if final parameters are actually valid
  else if(gsl_isnan(this->checkStateValidity())) {
    std::cout << "WARNING: at least one HMM had an invalid transition matrix from this BFGS parameter set" << std::endl;
    bestGuessOptim->optimSuccess = false;
  }
  else {
    bestGuessOptim->optimSuccess = true;
    bestGuessOptim->convertParamToProb(probs, gsl_multimin_fdfminimizer_x(s));
  }

  bestGuessOptim->setParamsToEst(probs);
  printRowVector(probs);
  bestGuessOptim->changeInBFGSLoglikelihood = initLikelihood - currLik;
  if(compareDoubles(0.0, bestGuessOptim->changeInBFGSLoglikelihood)) {
    std::cout << "WARNING: bestGuessOptim->changeInBFGSLoglikelihood = " << bestGuessOptim->changeInBFGSLoglikelihood << std::endl;
    bestGuessOptim->optimSuccess = false;
  }

  if(verbose) {
    printf("BFGS LOGLIKELIHOOD FOUND: %.40f\n", -currLik);
    printf("CHANGE IN BFGS LIKELIHOOD: %.40f\n", initLikelihood - currLik);
    printf("BFGS TOTAL TIME (sec): %.10f\n", totalTime);
    printf("DONE WITH BFGS.\n\n");
  }
  // clean up
  gsl_multimin_fdfminimizer_free(s);
  if(initGuess == nullptr) {
    gsl_vector_free(x);
  }
}

/*
 * function to call BFGS n times on the same HMM
 * from various starting points (mostly copied from the HMM.cpp implementation)
 * returns HMM with highest likelihood
 */
Optimizable* Optimizable::callBFGSNTimes(int numRuns, bool verbose, int seed) {
  if(this->generator == nullptr) {
    this->generator = new base_generator_type(seed);
  }
  double bestLL = GSL_NEGINF;
  double currLL = GSL_NEGINF;
  this->BFGSParamResults = new std::vector<gsl_vector*>();
  gsl_vector* currParams = nullptr;
  gsl_vector* bestParams = gsl_vector_alloc(this->getNumParamsToEst());
  gsl_vector* initGuess = gsl_vector_alloc(this->getNumParamsToEst());
  int bestRun = -1;
  int numSuccessfulRuns = 0;
  int numRetries = 0;
  bool needToRetry = false;
  bool ranOutOfRetries = false;
  for(int i = 1; i <= numRuns || needToRetry; i++) {
    // set initGuess
    if(needToRetry) {
      std::cout << "RETRYING RUN " << i << ", USING DEFAULT INITGUESS: ";
      this->setInitGuessNthTime(initGuess, i - 1 - numRetries, -numRetries);
    }
    else {
      std::cout << "CALLING BFGS, RUN " << i << " OF " << numRuns << " PLANNED, USING INITGUESS:" << std::endl;
      this->setInitGuessNthTime(initGuess, i - 1, numRuns);
    }
    printRowVector(initGuess);

    if(this->initGuessCopyBeforeBFGS == nullptr) {
      this->initGuessCopyBeforeBFGS = gsl_vector_alloc(initGuess->size);
    }
    gsl_vector_memcpy(this->initGuessCopyBeforeBFGS, initGuess);

    // run BFGS
    this->bfgs(initGuess, verbose);

    // if optim succeeded (according to flag) then add to count of successful runs
    if(this->optimSuccess) {
      numSuccessfulRuns++;
      currParams = gsl_vector_alloc(this->getNumParamsToEst());
      gsl_vector_memcpy(currParams, this->getParamsToEst());
      this->BFGSParamResults->push_back(currParams);

      // if have previously retried (ie bw starting point failed), run BFGS several times in case one of the predetermined starting points got stuck in a local max
      if(numRetries >= 1 && numSuccessfulRuns < 3) { // if 3 or more successful runs, break
        needToRetry = true;
      }
      else {
        needToRetry = false;
      }

      // if previously retried, keep tallying upwards
      if(numRetries > 0) {
        numRetries++;
      }
    }
    else {
      if(verbose) {
        std::cout << "WARNING: BFGS optimization failed for run " << i << std::endl;
      }
      // try again
      needToRetry = true;
      numRetries++;
      if(numRetries >= 16) { // one more than neg cases in setInitGuessNthTime
        std::cout << "WARNING: out of BFGS retries, exiting BFGS" << std::endl;
        ranOutOfRetries = true;
        break;
      }

      if(this->gradientDebug) {
        std::cerr << "##################################################################" << std::endl; // added for more granular printGradientPerIter
      }
    }
    currLL = this->getLogLikelihood();
    if(this->gradientDebug) {
      std::cerr << "##################################################################" << std::endl; // added for more granular printGradientPerIter
    }

    // save if currLL is better than what's stored before, and it was a successful run or we're out of runs
    if(currLL > bestLL && (this->optimSuccess || ranOutOfRetries)) {
      gsl_vector_memcpy(bestParams, this->getParamsToEst());
      bestLL = currLL;
      bestRun = i;

      // save the initGuess with the best optimized loglikelihood so far
      if(this->bestOptimLLInitGuess == nullptr) {
        this->bestOptimLLInitGuess = gsl_vector_alloc(initGuess->size);
      }
      gsl_vector_memcpy(this->bestOptimLLInitGuess, initGuess);

    }
    if(verbose) {
      printf("\n\n");
    }
  }
  if(numSuccessfulRuns == 0) {
    std::cerr << "ERROR: all starting param sets failed BFGS optimization. Saving last params" << std::endl;
    currParams = gsl_vector_alloc(this->getNumParamsToEst());
    gsl_vector_memcpy(currParams, this->getParamsToEst());
    this->BFGSParamResults->push_back(currParams);
    gsl_vector_memcpy(bestParams, this->getParamsToEst());
  }
  this->setParamsToEst(bestParams); // save best params
  printf("BEST OVERALL BFGS LIKELIHOOD (FROM RUN %i): %.40f\n", bestRun, this->getLogLikelihood());
  std::cout << "BEST BFGS PARAMS:" << std::endl;
  printColVector(bestParams);
  return this;
}

