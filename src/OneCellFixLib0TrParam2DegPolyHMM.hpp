#ifndef ONECELLFIXLIB0TRPARAM2DEGPOLYHMM_HPP
#define ONECELLFIXLIB0TRPARAM2DEGPOLYHMM_HPP

#include "OneCell0TrParam2DegPolyHMM.hpp"

/*
 * Library size, alpha/beta/lambda are all fixed. Est branch length
 *
 * this->paramsToEst = [t]
 * this->fixedParams = [lib, alpha, beta, lambda]
 */
class OneCellFixLib0TrParam2DegPolyHMM : public OneCell0TrParam2DegPolyHMM {
  protected:
    std::vector<gsl_matrix*>* totalLogEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]. Copied directly from TwoCellFixLib3TrParam2DegPolyHMM
    void updateTotalLogEmissionLookup();

    OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches);

  public:
    // constructors and destructor
    OneCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy);
    virtual ~OneCellFixLib0TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual void setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) override;

    // functions that depend on numbering and ordering of transition params
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual OneCellFixLib0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true) override;

    virtual void miscFunctions() override;
};

#endif

