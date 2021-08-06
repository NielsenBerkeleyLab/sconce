#ifndef ALLIND3TRPARAM2DEGPOLYHMM_HPP
#define ALLIND3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include <vector>
#include <unordered_map>

#include "Optimizable.hpp"
#include "OneCell3TrParam2DegPolyHMM.hpp"
#include "AllCells3TrParam2DegPolyHMM.hpp"

/*
 * This class is for making a collection of HMMs that jointly estimates beta and gammma across a given number of
 * tumor cells, where each cell is independent of the other cells (except for the shared transition parameters).
 * Branch length (t) and library size scaling factors are estimated for each cells.
 *
 * Because there are so many cells, it makes sense to make the likelihood calc (at each iter
 * in BFGS) parallelizable. So, this class will hold a collection of OneCell* HMMs where at each iteration
 * of BFGS, the appropriate parameters are set
 *
 * Because each cell is represented in exactly one HMM, all cellIdx and hmmIdx variables are the same.
 *
 * //this->paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t_cell0, t_cell1, ..., t_cellN]
 * this->paramsToEst = [lib0, lib1, ..., libN, beta, lambda, t_cell0, t_cell1, ..., t_cellN]
 * this->fixedParams = [alpha]
 */
class AllInd3TrParam2DegPolyHMM : public AllCells3TrParam2DegPolyHMM {
  private:
  protected:
    // protected ctors
    AllInd3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy);
    AllInd3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst);
    virtual void makeHMMs();
    virtual void makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams);

    // methods for least squares for solving for Baum Welch estimates
    virtual void baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual double baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) override;
    virtual void saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) override;

  public:
    // constants

    // constructors and destructor
    virtual ~AllInd3TrParam2DegPolyHMM();
    static AllInd3TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy);

    // accessors and mutators
    virtual void setLibScalingFactor(int cellIdx, double lib);
    virtual void setAllLibScalingFactors(double libScalingFactor);
    virtual double getLibScalingFactor(int cellIdx) const;
    virtual double setAllTransition(gsl_vector* transitionParams) override;

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;
    virtual void miscFunctions() override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;

};

#endif

