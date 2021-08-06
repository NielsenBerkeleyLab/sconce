#ifndef ALLINDFIXLIB3TRPARAM2DEGPOLYHMM_HPP
#define ALLINDFIXLIB3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include <vector>
#include <unordered_map>

#include "AllInd3TrParam2DegPolyHMM.hpp"
#include "OneCellFixLib3TrParam2DegPolyHMM.hpp"

/*
 * This class is for making a collection of HMMs that estimates transition params (beta/lambda) and branches for a given number of tumor cells.
 * libs (and alpha) are fixed
 *
 * Because each cell is represented in exactly one HMM, all cellIdx and hmmIdx variables are the same.
 *
 * this->paramsToEst = [beta, lambda, t_cell0, t_cell1, ..., t_cellN]
 * this->fixedParams = [lib0, lib1, ..., libN, alpha]
 */
class AllIndFixLib3TrParam2DegPolyHMM : public AllInd3TrParam2DegPolyHMM {
  private:
  protected:
    // protected ctors
    AllIndFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    AllIndFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst);
    using AllInd3TrParam2DegPolyHMM::makeHMMs;
    virtual void makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams);

  public:
    // constants

    // constructors and destructor
    virtual ~AllIndFixLib3TrParam2DegPolyHMM();
    static AllIndFixLib3TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    // accessors and mutators
    virtual void setLibScalingFactor(int cellIdx, double lib);
    virtual void setAllLibScalingFactors(double libScalingFactor);
    virtual double getLibScalingFactor(int cellIdx) const;

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setSimParamsToEst(gsl_vector* params) override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;

};

#endif

