#ifndef ONECELL3TRPARAM2DEGPOLYHMM_HPP
#define ONECELL3TRPARAM2DEGPOLYHMM_HPP

#include "HMM.hpp"
#include <gsl/gsl_statistics_double.h>

/*
 * This class analyzes 1 cell for CNVs, by allowing library size, beta/lambda, and t to vary.
 * alpha is fixed.
 *
 * this->paramsToEst = [lib, beta, lambda, t]
 * this->fixedParams = [alpha]
 */
class OneCell3TrParam2DegPolyHMM : public HMM {
  private:
    // member variables
    std::vector<double>* logFacKVec; // k:[log(k!)]
    double getLogFacK(int k);
    void setLogFacK();

  protected:
    OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches);
    double setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t);

  public:
    // constructors and destructor
    OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy);
    virtual ~OneCell3TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void print(FILE* stream) override;
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    virtual double getAlpha() const;
    virtual void setAlpha(double alpha);

    // functions that depend on numbering and ordering of transition params
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) override;
    using HMM::setTransition;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual double getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual double getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual OneCell3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true) override;
    virtual void simulate() override;
    virtual void simulate(int seed) override;
    virtual void simulate(int seed, bool simDiploid, double diploid_lambda_i = 272.5568, int numDiploidCells = 35);
    virtual void setUpBaumWelchLeastSquares();
    virtual double baumWelchLeastSquares_f(gsl_vector* probs);

    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;

    virtual double checkOptimProbValidity(gsl_vector* probs) const override;
    
};

#endif

