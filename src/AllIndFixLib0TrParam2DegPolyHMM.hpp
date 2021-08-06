#ifndef ALLINDFIXLIB0TRPARAM2DEGPOLYHMM_HPP
#define ALLINDFIXLIB0TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include <vector>
#include <unordered_map>

#include "AllInd0TrParam2DegPolyHMM.hpp"
#include "OneCellFixLib0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making a collection of HMMs that estimates branches for a given number of tumor cells.
 * //transition params (beta/gamma) and libs are fixed.
 * transition params (beta/lambda) and libs are fixed.
 *
 * Because each cell is represented in exactly one HMM, all cellIdx and hmmIdx variables are the same.
 *
 * this->paramsToEst = [t_cell0, t_cell1, ..., t_cellN]
 * //this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, gamma]
 * this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, lambda]
 */
class AllIndFixLib0TrParam2DegPolyHMM : public AllInd0TrParam2DegPolyHMM {
  private:
  protected:
    // protected ctors
    AllIndFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    AllIndFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst);
    //virtual void makeHMMs();
    using AllInd0TrParam2DegPolyHMM::makeHMMs; // unhide parent method of same name https://stackoverflow.com/a/18100999
    virtual void makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams);

  public:
    // constants

    // constructors and destructor
    virtual ~AllIndFixLib0TrParam2DegPolyHMM();
    static AllIndFixLib0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    // accessors and mutators
    virtual void setLibScalingFactor(int cellIdx, double lib);
    virtual void setAllLibScalingFactors(double libScalingFactor);
    virtual double getLibScalingFactor(int cellIdx) const;
    //virtual double setAllTransition(gsl_vector* transitionParams) override;
    //virtual double setAllBranches(gsl_vector* branches);

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;

};

#endif



