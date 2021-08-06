#ifndef ONECELL0TRPARAM2DEGPOLYHMM_HPP
#define ONECELL0TRPARAM2DEGPOLYHMM_HPP

#include "OneCell3TrParam2DegPolyHMM.hpp"

/*
 * One cell, all transition params (alpha/beta/gamma) are fixed. Estimate lib and t
 *
 * this->paramsToEst = [lib, t]
 * //this->fixedParams = [alpha, beta, gamma]
 * this->fixedParams = [alpha, beta, lambda]
 */
class OneCell0TrParam2DegPolyHMM : public OneCell3TrParam2DegPolyHMM {
  protected:
    OneCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches);

  public:
    // constructors and destructor
    OneCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy);
    virtual ~OneCell0TrParam2DegPolyHMM();

    // accessors and mutators

    // functions that depend on numbering and ordering of transition params
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) override;
    using HMM::setTransition;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual OneCell0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true) override;
};

#endif

