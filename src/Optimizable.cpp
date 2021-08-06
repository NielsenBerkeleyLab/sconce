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
  this->gradientDebug = true; // TODO this should be false by default. need to set in HMM and base classes
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

  //gsl_vector* oldParamsToEst = gsl_vector_alloc(optimObj->getNumParamsToEst());
  //gsl_vector_memcpy(oldParamsToEst, optimObj->getParamsToEst());

  //std::cerr << "IN EVALLIKELIHOODATPOINT" << std::endl;
  //return GSL_NAN;
  //HMM* hmm = (HMM*) params; // this hmm instance can (and should) be changed/updated

  // convert vector v (bfgs params point to evaluate) back to probs
  //gsl_vector* probs = gsl_vector_alloc(v->size);
  gsl_vector* probs = optimObj->probParamConversionVec;
  //hmm->convertParamToProb(probs, v);
  optimObj->convertParamToProb(probs, v);
  //std::cerr << "in evalLL at point, probs: " << std::endl;
  //printColVector(probs);

  // if can't find steady state dist using this param set, return NAN
  //double status = hmm->setParamsToEst(probs); // includes a call to this->setTransition()
  double status = optimObj->checkOptimProbValidity(probs);
  if(gsl_isnan(status)) {
  /*if(std::abs(status - 0) > 1e-4) { // debugging penalty
    //std::cout << "optimObj->checkOptimProbValidity(probs) RETURNED NAN" << std::endl;
    //printColVector(probs);
    fprintf(stderr, "~");
    for(unsigned int i = 0; i < probs->size; i++) {
      fprintf(stderr, "\t%.15f", gsl_vector_get(probs, i));
    }
    fprintf(stderr, "\n");
    return optimObj->initLl * pow(10, status); // debugging penalty*/
    //printColVector((gsl_vector*) v);
    //return status;
    //return std::numeric_limits<double>::max(); // this seems to fix the line search problem where BFGS never progesses, still testing Wed 23 Oct 2019 03:42:51 PM PDT

    // reset to old params if this failed Wed 19 May 2021 03:49:30 PM PDT don't reset, it's a waste of time
    //std::cout << "WOULD HAVE RESET TO OLD PARAMS, NEW PARAMS FAILED checkOptimProbValidity, returning nan" << std::endl;
    //optimObj->setParamsToEst(oldParamsToEst);
    //gsl_vector_free(oldParamsToEst);
    //return std::numeric_limits<double>::max() / 1e50; // max() is ~1.7e310. This goes to inf sometimes in gradient calculations if you go too small (ie add f(x) in finite difference method) Tue 19 Nov 2019 05:27:49 PM PST
    return status;
  }
  //std::cout << "checkOptimProbValidity passed" << std::endl;
  status = optimObj->setParamsToEst(probs); // includes a call to this->setTransition()


  // Fri 28 Feb 2020 02:19:09 PM PST
  // try adjusting the lib size scaling factors according to the stationary dist
  //optimObj->miscFunctions(); // statDist

  //std::cerr << "v, probs, paramsToEst, trans" << std::endl;
  //printRowVector((gsl_vector*)v);
  //printRowVector(probs);
  //printRowVector(optimObj->getParamsToEst());
  //printMatrix(hmm->getTransition());

  //std::cout << "#####################" << std::endl;
  //optimObj->print(stdout);

  //std::cout << "status: " << status << std::endl;

  if(gsl_isnan(status)) {
    // reset to old params if this failed Wed 19 May 2021 03:49:30 PM PDT don't reset, it's a waste of time
    //std::cout << "WOULD HAVE RESET TO OLD PARAMS, NEW PARAMS RESULTED IN NAN STATUS, returning nan" << std::endl;
    //optimObj->setParamsToEst(oldParamsToEst); // Wed 19 May 2021 11:43:01 AM PDT TEST 1
    //gsl_vector_free(oldParamsToEst);
    //std::cerr << "Optimizable::evalLikelihoodAtPoint caught nan, returning nan" << std::endl;
    return status;
  }
  //std::cout << "here" << std::endl;

  // call forward alg
  //double negProbObsSeq = -hmm->runForwardAlg(); // negate because gsl provides a minimizer, and we want to maximize the likelihood
  double loglikelihood = -optimObj->getLogLikelihood(); // negate because gsl provides a minimizer, and we want to maximize the likelihood
  //std::cout << "loglikelihood: " << loglikelihood << std::endl;



  /*// added for more granular printGradientPerIter
  //std::cerr << "optimObj->gradientDebug: " << optimObj->gradientDebug << std::endl;
  if(optimObj->gradientDebug) {
    fprintf(stderr, "-%.20f\t", loglikelihood);
    //fprintf(stderr, "-%.20f\t", f_x_ph);
    //fprintf(stderr, "-%.20f\t", f_x_mh);
    //fprintf(stderr, "%.20f\t", derivApprox);
    fprintf(stderr, "0\t");
    fprintf(stderr, "0\t");
    fprintf(stderr, "0\t");
    for(unsigned int j = 0; j < probs->size; j++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(probs, j));
    }
    for(unsigned int j = 0; j < v->size; j++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(v, j));
    }
    for(unsigned int j = 0; j < v->size; j++) {
      //fprintf(stderr, "%.20f\t", gsl_vector_get(df, j));
      fprintf(stderr, "0\t");
    }
    fprintf(stderr, "\n");
  }*/



  //gsl_vector_free(oldParamsToEst);
  return loglikelihood;
}

/* 
 * evalGradientAtPoint evaluates the gradient at point v using the finite difference method, and stores
 * the gradient in df
 */
void Optimizable::evalGradientAtPoint(const gsl_vector* v, void* params, gsl_vector* df) {
  //std::cerr << "Optimizable::evalGradientAtPoint" << std::endl;
  // approx using finite difference method:
  // f'(x) = [f(x + h) - f(x)] / h
  //double h = 1e-2; // step size
  //double h = 1e-3; // step size
  double h = 1e-4; // step size // oh4
  //double h = 1e-5; // step size // oh4
  //double h = 1e-6; // step size ORIG
  //double h = 1e-8; // step size from testOptim // ORIG Fri 28 May 2021 12:14:08 PM PDT
  double f_x = evalLikelihoodAtPoint(v, params);
  if(gsl_isnan(f_x)) {
    //std::cerr << "##### f_x is nan" << std::endl;
    gsl_vector_set_all(df, GSL_NAN);

    gsl_vector* probVec = gsl_vector_alloc(v->size); // added for printGradientPerIter
    ((Optimizable*)params)->convertParamToProb(probVec, v);
    fprintf(stderr, "%.20f\t", f_x);
    //fprintf(stderr, "-%.20f\t", f_x_ph);
    //fprintf(stderr, "-%.20f\t", f_x_mh);
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
    //x_h = gsl_vector_get(v, i) + h;
    x = gsl_vector_get(v, i);
    gsl_vector_set(currVec, i, x + h);
    f_x_ph = evalLikelihoodAtPoint(currVec, params);

    if(gsl_isnan(f_x_ph)) {
      derivApprox = GSL_NAN;
    } else {
      //derivApprox = (f_x_ph - f_x) / h; // forward difference

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
    //fprintf(stderr, "-%.20f\t", f_x_ph);
    //fprintf(stderr, "-%.20f\t", f_x_mh);
    //fprintf(stderr, "0\t"); // Fri 31 Jan 2020 04:55:35 PM PST don't want the plotting scripts to get messed up
    //fprintf(stderr, "%.20f\t", derivApprox);
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
  //*f = this->evalLikelihoodAtPoint(v); 
  //this->evalGradientAtPoint(v, df);
}

/*
 * function to estimate parameters using bfgs and forward algorithm.
 * bestGuessOptim is the dest HMM to copy into; should be specific to the subclass.
 * initGuess is in BFGS space
 */
void Optimizable::bfgs(gsl_vector* initGuess, Optimizable* bestGuessOptim, bool verbose) const {
  //std::cout << "Optimizable::bfgs:" << std::endl;
  //printColVector(bestGuessOptim->getParamsToEst());
  //printColVector(initGuess);

  // BFGS version
  gsl_multimin_function_fdf my_func; // which function to minimize
  my_func.df = &evalGradientAtPoint; // gradient of function
  my_func.fdf = &evalLikelihoodGradAtPoint; // how to set both f and df
  //*/

  /*
  // simplex version
  gsl_multimin_function my_func; // which function to minimize
  */

  my_func.f = &evalLikelihoodAtPoint; // function itself
  my_func.params = bestGuessOptim; // this optimizable should be passed around


  // Fri 07 Feb 2020 09:50:37 PM PST trying out new X_ij ~ NB(r_ik + epsilon)
  //bestGuessOptim->miscFunctions();



  // Starting point/initial guess, in BFGS space
  gsl_vector* x = nullptr;
  gsl_vector* probs = nullptr;
  int status = GSL_CONTINUE;
  double initLikelihood = 0;

  if(initGuess != nullptr) {
    // before doing anything, check if probs are valid
    probs = gsl_vector_alloc(initGuess->size);
    bestGuessOptim->convertParamToProb(probs, initGuess);

    /*std::cout << "initGuess in BFGS(" << std::endl;
    printColVector(initGuess);
    std::cout << "probs in BFGS(" << std::endl;
    printColVector(probs);*/

    status = bestGuessOptim->checkOptimProbValidity(probs);
    if(isnan(status)) {
      std::cerr << "ERROR: initial param set is invalid, not running optimization" << std::endl;
      std::cerr << "##################################################################" << std::endl;
      bestGuessOptim->optimSuccess = false;
      return;
    }

    x = initGuess;
    my_func.n = initGuess->size;

    //bestGuessOptim->miscFunctions(); // update libs for viterbi and statDist

    initLikelihood = Optimizable::evalLikelihoodAtPoint(x, bestGuessOptim);
    //bestGuessOptim->print(stdout);
    bestGuessOptim->initLl = initLikelihood; // Sat 23 May 2020 12:15:43 PM PDT debugging penalty
    if(verbose) {
      printf("BFGS INITIAL LIKELIHOOD: %.40f\n", -initLikelihood);

      /*// Fri 06 Mar 2020 04:07:05 PM PST
      // check that the same likelihood is reproduced for the same params at each iteration of bfgs (ie make sure background structures are being correctly set)
      gsl_vector* tmpProbs = gsl_vector_alloc(bestGuessOptim->getNumParamsToEst());
      gsl_vector* tmpParams = gsl_vector_alloc(bestGuessOptim->getNumParamsToEst());

      gsl_vector_set(tmpProbs, 0, 0.01); // beta
      gsl_vector_set(tmpProbs, 1, 0.02); // gamma

      // pairwise branch lengths
      for(int pairIdx = 0; pairIdx < 20; pairIdx++) {
        gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 0, 0.1); // set t1
        gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 1, 0.2); // set t2
        gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 2, 0.3); // set t3
      }

      bestGuessOptim->convertProbToParam(tmpParams, tmpProbs);
      double tmpLl = evalLikelihoodAtPoint(tmpParams, bestGuessOptim);
      fprintf(stderr, "#\t-%0.20f\n", tmpLl);
      gsl_vector_free(tmpProbs);
      gsl_vector_free(tmpParams);*/

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
    //convertProbToParam(x, bestGuessOptim->paramsToEst, bestGuessOptim->getKploidy());
    convertProbToParam(x, bestGuessOptim->paramsToEst);
  }

  // BFGS version
  const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2; // which algorithm to use
  //const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr; // which algorithm to use
  //const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_pr; // which algorithm to use
  gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, my_func.n);
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-6, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-4, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance) // ORIG
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-2, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 1, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 10, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 100, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  //gsl_multimin_fdfminimizer_set(s, &my_func, x, 1000, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance)
  //gsl_multimin_fdfminimizer_set (s, &my_func, x, 1e-1, .1); // from testOptim.cpp

  /*
  // simplex version
  const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2; // which algorithm to use
  gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc(T, my_func.n);
  gsl_vector* stepSize = gsl_vector_alloc(my_func.n);
  gsl_vector_set_all(stepSize, 1);
  gsl_multimin_fminimizer_set(s, &my_func, x, stepSize); // minimizer s, function my_func, starting at x, size of initial trial steps
  */

  int iter = 0;
  double oldLik = initLikelihood;//std::numeric_limits<double>::max();
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
  //int countTooLong = 0;
  double totalTime = 0;


  /*// Thu 12 Sep 2019 01:19:21 PM PDT want to print a big table of gradient values at each BFGS iteration
  // used for printGradientPerIter (bfgs iters only)
  // ll | {HMM probs} | {bfgs params} | {gradients}
  // convert vector v (bfgs params point to evaluate) back to probs
  fprintf(stderr, "-%.20f\t", gsl_multimin_fdfminimizer_minimum(s));
  gsl_vector* v = gsl_multimin_fdfminimizer_x(s);
  gsl_vector* df = gsl_multimin_fdfminimizer_gradient(s);
  //bestGuessOptim->convertParamToProb(probs, v, bestGuessOptim->getKploidy());
  bestGuessOptim->convertParamToProb(probs, v);
  for(unsigned int i = 0; i < probs->size; i++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(probs, i));
  }
  for(unsigned int i = 0; i < v->size; i++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(v, i));
  }
  for(unsigned int i = 0; i < df->size; i++) {
    fprintf(stderr, "%.20f\t", gsl_vector_get(df, i));
  }
  //fprintf(stderr, "%.20f\t", gsl_multimin_fminimizer_size((const gsl_multimin_fminimizer*)s)); // simplex
  fprintf(stderr, "\n");*/

  while(status == GSL_CONTINUE) {
    /*// Fri 06 Mar 2020 04:07:05 PM PST
    // check that the same likelihood is reproduced for the same params at each iteration of bfgs (ie make sure background structures are being correctly set)
    gsl_vector* tmpProbs = gsl_vector_alloc(bestGuessOptim->getNumParamsToEst());
    gsl_vector* tmpParams = gsl_vector_alloc(bestGuessOptim->getNumParamsToEst());

    gsl_vector_set(tmpProbs, 0, 0.01); // beta
    gsl_vector_set(tmpProbs, 1, 0.02); // gamma

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < 20; pairIdx++) {
      gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 0, 0.1); // set t1
      gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 1, 0.2); // set t2
      gsl_vector_set(tmpProbs, 2 + 3 * pairIdx + 2, 0.3); // set t3
    }

    bestGuessOptim->convertProbToParam(tmpParams, tmpProbs);
    double tmpLl = evalLikelihoodAtPoint(tmpParams, bestGuessOptim);
    fprintf(stderr, "#\t-%0.20f\n", tmpLl);
    gsl_vector_free(tmpProbs);
    gsl_vector_free(tmpParams);*/




    begin = std::chrono::steady_clock::now();

    // BFGS
    status = gsl_multimin_fdfminimizer_iterate(s); // do one iter
    /*
    // simplex
    status = gsl_multimin_fminimizer_iterate(s); // do one iter
    */

    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    iter++;
    totalTime += elapsedSec;
    // BFGS
    currLik = gsl_multimin_fdfminimizer_minimum(s);
    deltaLik = std::abs(currLik - oldLik);
    //*/
    /*
    // simplex
    currLik = gsl_multimin_fminimizer_minimum(s);
    */

    // BFGS
    //// check gradient
    //if(verbose) {
    //  printRowVector(stderr, gsl_multimin_fdfminimizer_gradient(s));
    //}

    //*/

    // print update message
    if(verbose) {
      printf("ON BFGS ITER %d, likelihood %.20f, time elapsed (sec) %.5f, iter change in loglikelihood %.5f, total change in loglikelihood %.5f\n", iter, -currLik, elapsedSec, deltaLik, initLikelihood - currLik);
      //printf("ON ITER %d, likelihood %.40f, status %i, time elapsed (sec) %.5f\n", iter, -currLik, status, elapsedSec);
    }

    // if any error has occurred, break
    if(status) {
      std::cout << "STATUS IS: " << gsl_strerror(status) << std::endl;
      break;
    }

    // if currLik is NaN, stop
    if(gsl_isnan(currLik)) {
      printf("STATUS IS: ERROR: current likelihood is %f\n", currLik);
      //std::cout << "Using previously best params instead:" << std::endl;
      //printRowVector(prevBestParams);
      break;
    }

    // BFGS
    // check if have converged (is the gradient close enough to 0?)
    //status = gsl_multimin_test_gradient(s->gradient, 1e-8);
    //status = gsl_multimin_test_gradient(s->gradient, 1e-6); ORIG
    status = gsl_multimin_test_gradient(s->gradient, 1e-4);
    //*/
    /*
    // simplex
    // check if have converged (is the size close enough to tolerance?)
    status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), 1e-6); // TODO how to pick epsabs?
    */
    if(status == GSL_SUCCESS) {
      printf("STATUS IS: SUCCESS: Minimum found.\n");
      //printRowVector(gsl_multimin_fminimizer_x(s));
    }

    // check if done too many iters
    //if(iter >= 500000) { // simplex
    //if(iter >= 5000) { // simplex
    if(iter >= 500) { // BFGS
    //if(iter >= 100) {
    //if(iter >= 50) {
    //if(iter >= 15) {
    //if(iter >= 10) {
    //if(iter >= 5) {
    //if(iter >= 2) {
    //if(iter >= 1) {
      std::cout << "STATUS IS: minimum not found in 500 iterations" << std::endl;
      break;
    }

    // check if have converged based on change in likelihood
    //if(deltaLik < 1e-4) {
    //if(deltaLik < 1e-5) {
    if(deltaLik < 1e-6) { // ctc6
    //if(deltaLik < 1e-8) {
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
    //std::cout << "HERE, " << currLik << std::endl;
    //bestGuessOptim->convertParamToProb(probs, prevBestParams);
    //printRowVector(probs);
    //printRowVector(prevBestParams);
    //bestGuessOptim->print(stdout);

    //std::cout << "SETTING PROBS EXPLICITLY" << std::endl;
    //std::cerr << "SETTING PROBS EXPLICITLY" << std::endl;
    //std::cout << "EXPLICIT SET STATUS: " << bestGuessOptim->setParamsToEst(probs) << std::endl;
    //std::cerr << "DONE SETTING PROBS EXPLICITLY" << std::endl;
    //bestGuessOptim->print(stdout);

    // check if too many consecutive iters are taking too long (tends to happen if bfgs optim gets stuck)
    //if(elapsedSec > 250) { // 250 chosen arbitrarily after some empirical evidence
    /*if(elapsedSec > 2500) { // TODO set a longer timing threshold for allPairs. maybe make a member variable?
      countTooLong++;
      if(countTooLong >= 2) {
        std::cout << "STATUS IS: too many consecutive iterations taking too long" << std::endl;
        break;
      }
    } else {
      countTooLong = 0;
    }*/

    
    /*// Thu 12 Sep 2019 01:19:21 PM PDT want to print a big table of gradient values at each BFGS iteration
    // used for printGradientPerIter (bfgs iters only)
    // ll | {HMM probs} | {bfgs params} | {gradients}
    // convert vector v (bfgs params point to evaluate) back to probs
    fprintf(stderr, "-%.20f\t", gsl_multimin_fdfminimizer_minimum(s));
    gsl_vector* v = gsl_multimin_fdfminimizer_x(s);
    gsl_vector* df = gsl_multimin_fdfminimizer_gradient(s);
    //bestGuessOptim->convertParamToProb(probs, v, bestGuessOptim->getKploidy());
    bestGuessOptim->convertParamToProb(probs, v);
    for(unsigned int i = 0; i < probs->size; i++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(probs, i));
    }
    for(unsigned int i = 0; i < v->size; i++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(v, i));
    }
    for(unsigned int i = 0; i < df->size; i++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(df, i));
    }
    //fprintf(stderr, "%.20f\t", gsl_multimin_fminimizer_size((const gsl_multimin_fminimizer*)s)); // simplex
    fprintf(stderr, "\n");*/
    

    // Fri 07 Feb 2020 09:50:37 PM PST trying out new X_ij ~ NB(r_ik + epsilon)
    //bestGuessOptim->print(stdout);
    // try updating the lib scaling factors
    //bestGuessOptim->miscFunctions(); // viterbi version




  }
  /*// used for printGradientPerIter (bfgs iters only)
  std::cerr << "##################################################################" << std::endl;*/
  //std::cout << std::endl;
  if(bestGuessOptim->gradientDebug) {
    std::cerr << "##################################################################" << std::endl; // added for more granular printGradientPerIter; used when bfgs is only called once
  }

  // save results
  if(gsl_isnan(currLik)) {
    //*(bestGuessOptim->optimSuccess) = false;
    bestGuessOptim->optimSuccess = false;

    // use the previously best results
    if(verbose) {
      std::cout << "Using previously best params" << std::endl;
      std::cerr << "ERROR: Optimzation failed: likelihood is nan, using previously best params" << std::endl;
    }
    bestGuessOptim->convertParamToProb(probs, prevBestParams);
    //printRowVector(prevBestParams);
    //printRowVector(probs);
    //bestGuessOptim->print(stdout);
    //if(verbose) {
    //  std::cerr << "ERROR: Optimzation failed: likelihood is nan" << std::endl;
    //}
  }
  // check if final parameters are actually valid
  else if(gsl_isnan(this->checkStateValidity())) {
    std::cout << "WARNING: at least one HMM had an invalid transition matrix from this BFGS parameter set" << std::endl;
    bestGuessOptim->optimSuccess = false;
  }
  else {
    //*(bestGuessOptim->optimSuccess) = true;
    bestGuessOptim->optimSuccess = true;
    // BFGS
    //bestGuessOptim->convertParamToProb(probs, gsl_multimin_fdfminimizer_x(s), bestGuessOptim->getKploidy());
    bestGuessOptim->convertParamToProb(probs, gsl_multimin_fdfminimizer_x(s));
    //*/
    /*
    // simplex
    //bestGuessOptim->convertParamToProb(probs, gsl_multimin_fminimizer_x(s), bestGuessOptim->getKploidy());
    bestGuessOptim->convertParamToProb(probs, gsl_multimin_fminimizer_x(s));
    */
  }

  //std::cout << "BESTGUESSOPTIM AFTER SETTING PROBS BEFORE ENDING" << std::endl;
  bestGuessOptim->setParamsToEst(probs);
  printRowVector(probs);
  bestGuessOptim->changeInBFGSLoglikelihood = initLikelihood - currLik;
  if(compareDoubles(0.0, bestGuessOptim->changeInBFGSLoglikelihood)) {
    std::cout << "WARNING: bestGuessOptim->changeInBFGSLoglikelihood = " << bestGuessOptim->changeInBFGSLoglikelihood << std::endl;
    bestGuessOptim->optimSuccess = false;
  }
  //bestGuessOptim->print(stdout);

  if(verbose) {
    printf("BFGS LOGLIKELIHOOD FOUND: %.40f\n", -currLik);
    printf("CHANGE IN BFGS LIKELIHOOD: %.40f\n", initLikelihood - currLik);
    printf("BFGS TOTAL TIME (sec): %.10f\n", totalTime);
    printf("DONE WITH BFGS.\n\n");
  }
  // clean up
  // BFGS
  gsl_multimin_fdfminimizer_free(s);
  //*/
  /*
  // simplex
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(stepSize);
  */
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
  //Optimizable* bestGuess = nullptr;
  //Optimizable* currGuess = nullptr;
  int bestRun = -1;
  int numSuccessfulRuns = 0;
  int numRetries = 0;
  bool needToRetry = false;
  bool ranOutOfRetries = false;
  //for(int i = 0; numSuccessfulRuns < numRuns && i < this->getMaxNumBFGSStarts(); i++) {
  //for(int i = 0; i < numRuns && i < this->getMaxNumBFGSStarts(); i++) {
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
    //currGuess = this->bfgs(initGuess, verbose);
    this->bfgs(initGuess, verbose);

    // if optim succeeded (according to flag) then add to count of successful runs
    //if(currGuess->optimSuccess) {
    if(this->optimSuccess) {
      //needToRetry = false;
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
      //continue;

      if(this->gradientDebug) {
        std::cerr << "##################################################################" << std::endl; // added for more granular printGradientPerIter
      }

      //// reset initGuess
      //this->setInitGuessNthTime(initGuess, i - 1, -1);
      //if(verbose) {
      //  std::cout << "RETRYING RUN " << i << ", USING DEFAULT INITGUESS: ";
      //  printRowVector(initGuess);
      //}

      //// run BFGS
      ////currGuess = this->bfgs(initGuess, verbose);
      //this->bfgs(initGuess, verbose);
      //gsl_vector_memcpy(currParams, this->getParamsToEst()); // currParams is already stored
      ////this->BFGSParamResults->push_back(currParams);

    }
    //currLL = currGuess->runForwardAlg();
    //currLL = currGuess->getLogLikelihood();
    currLL = this->getLogLikelihood();
    //if(currGuess->gradientDebug) {
    if(this->gradientDebug) {
      std::cerr << "##################################################################" << std::endl; // added for more granular printGradientPerIter
    }

    //if(verbose) {
    //  std::cout << "OPTIMIZED BFGS PARAMS FOR RUN " << i << " FOUND: ";
    //  //printRowVector(currGuess->getParamsToEst());
    //  printRowVector(this->getParamsToEst());
    //}
    //if(currLL > bestLL) {
    // save if currLL is better than what's stored before, and it was a successful run or we're out of runs
    if(currLL > bestLL && (this->optimSuccess || ranOutOfRetries)) {
      //bestGuess = currGuess;
      gsl_vector_memcpy(bestParams, this->getParamsToEst());
      bestLL = currLL;
      bestRun = i;

      // save the initGuess with the best optimized loglikelihood so far
      //bestGuess->bestOptimLLInitGuess = gsl_vector_alloc(initGuess->size);
      //gsl_vector_memcpy(bestGuess->bestOptimLLInitGuess, initGuess);
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
  //gsl_vector_free(initGuess);
  this->setParamsToEst(bestParams); // save best params
  printf("BEST OVERALL BFGS LIKELIHOOD (FROM RUN %i): %.40f\n", bestRun, this->getLogLikelihood());
  std::cout << "BEST BFGS PARAMS:" << std::endl;
  printColVector(bestParams);
  //return bestGuess;
  return this;
}

