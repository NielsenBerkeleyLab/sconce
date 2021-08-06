#ifndef ONECELLFIXLIB3TRPARAM2DEGPOLYHMM_HPP
#define ONECELLFIXLIB3TRPARAM2DEGPOLYHMM_HPP

#include "OneCell3TrParam2DegPolyHMM.hpp"

/*
 * This class is uesd in the second stage of AllPairs2Stages3TrParam2DegPolyHMM
 * by the AllPairs0TrParam2DegPolyHMM class.
 *
 * It's the same as OneCell3TrParam2DegPolyHMM but the library size is fixed
 *
 * this->paramsToEst = [beta, lambda, t]
 * this->fixedParams = [lib, alpha]
 */
class OneCellFixLib3TrParam2DegPolyHMM : public OneCell3TrParam2DegPolyHMM {
  protected:
    OneCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches);
    std::vector<gsl_matrix*>* totalLogEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]
    void updateTotalLogEmissionLookup();

  public:
    // constructors and destructor
    OneCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy);
    virtual ~OneCellFixLib3TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual void setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) override;

    // functions that depend on numbering and ordering of transition params
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual OneCellFixLib3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true) override;
};

#endif

