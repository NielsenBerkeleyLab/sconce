#ifndef OPTIMIZABLE_HPP
#define OPTIMIZABLE_HPP

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <iostream>
#include <chrono>

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
// https://www.boost.org/doc/libs/1_70_0/libs/random/example/random_demo.cpp
// This is a typedef for a random number generator. Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

#include "util.hpp"

/*
 * generic class that can optimize a function and estimate a set of parameters using BFGS.
 * Functions to be optimzed are implemented in subclasses in
 * evalLikelihoodAtPoint, evalGradientAtPoint, and evalLikelihoodGradAtPoint
 */
class Optimizable {
  protected:
    double initLl = 0; // Sat 23 May 2020 12:13:52 PM PDT debugging, for saving the init likelihood as a penalty term
    gsl_vector* paramsToEst; // params to estimate with BFGS; indicies into vector laid out in subclasses
    gsl_vector* fixedParams; // fixed parameters, such as lib sizes or alpha/beta/gamma, depending on subclass. Ordering is the same as paramsToEst
    bool optimSuccess; // flag to indicate if the bfgs optimization failed or not. TODO consider changing this to an int to save more detailed statuses, rather than a binary flag
    int maxNumBFGSStarts; // maximum number of starting points for optimization, usually set by setInitGuessNthTime
    double changeInBFGSLoglikelihood; // total change in loglikelihood after BFGS
    gsl_vector* initGuessCopyBeforeBFGS; // copy of initGuess, before BFGS. useful for keeping variables from changing too much from initial values
    gsl_vector* bestOptimLLInitGuess; // initGuess with best optimized loglikelihood, found via calling BFGS N times
    std::vector<gsl_vector*>* BFGSParamResults; // vector of paramsToEst from BFGS, starting from different parameter sets, stored in order discovered

    gsl_vector* probParamConversionVec; // intermediate for evalLikelihoodAtPoint. Hold param conversion somewhere instead of constantly allocating new memory

    gsl_vector* simParamsToEst; // copies of paramsToEst and fixedParams, useful for comparing to simualation params if we have them
    gsl_vector* simFixedParams;

  public:
    // random number generator, useful for simulations
    base_generator_type* generator;
    bool gradientDebug; // debugging flag; will print granular gradient information if true

    Optimizable();
    Optimizable(const Optimizable& other);

    // methods that must be implemented by subclasses because they depend on the model
    virtual int getMaxNumBFGSStarts() const;
    virtual double setParamsToEst(gsl_vector* params) = 0;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const = 0;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const = 0;
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const = 0;
    virtual Optimizable* bfgs(gsl_vector* initGuess, bool verbose = true) = 0; // this one should construct a subclass Optimizable, convert initGuess into BFGS space, and call Optimizable::bfgs
    virtual double getLogLikelihood() = 0;
    virtual void print(FILE* stream) = 0;
    virtual double checkOptimProbValidity(gsl_vector* probs) const = 0;
    virtual double checkStateValidity(double epsilon = 1e-8) const = 0; // should check if transition matrices are ok
    virtual void setSimParamsToEst(gsl_vector* params) = 0;
    virtual void setSimFixedParams(gsl_vector* params) = 0;

    static double evalLikelihoodAtPoint(const gsl_vector* v, void* params);
    static void evalGradientAtPoint(const gsl_vector* v, void* params, gsl_vector* df); // TODO evalGradientAtPoint and evalLikelihoodGradAtPoint are the same in HMM and AllPairs_. consider moving into Optimizable? ==> also evalLikelihoodAtPoint turns out to be the same? would need to add a checkParams method. this would also conveniently get rid of the problem with the non-static versions changing the calling object
    static void evalLikelihoodGradAtPoint(const gsl_vector* v, void* params, double* f, gsl_vector* df);

    // shared implementations
    virtual gsl_vector* getParamsToEst() const;
    virtual gsl_vector* getFixedParams() const;
    virtual int getNumParamsToEst() const;
    virtual int getNumFixedParams() const;
    virtual gsl_vector* getBestLLInitGuess() const; // get the initGuess that had the best final LL
    virtual std::vector<gsl_vector*>* getBFGSParamResults() const;
    void bfgs(gsl_vector* initGuess, Optimizable* bestGuessOptim, bool verbose = true) const; // bestGuessOptim is the dest Optimizable to copy into; should be specific to the class
    Optimizable* callBFGSNTimes(int numRuns, bool verbose = true, int seed = 43);

    // need a dummy function for testing
    virtual void miscFunctions() = 0;
};

#endif

