#ifndef ALLIND0TRPARAM2DEGPOLYHMM_HPP
#define ALLIND0TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include <vector>
#include <unordered_map>

#include "AllInd3TrParam2DegPolyHMM.hpp"
#include "OneCell0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making a collection of HMMs that estimates libs and branches for a given number of tumor cells.
 * //transition params (beta/gamma) are fixed.
 * transition params (beta/lambda) are fixed.
 *
 * Because each cell is represented in exactly one HMM, all cellIdx and hmmIdx variables are the same.
 *
 * this->paramsToEst = [lib0, lib1, ..., libN, t_cell0, t_cell1, ..., t_cellN]
 * //this->fixedParams = [alpha, beta, gamma]
 * this->fixedParams = [alpha, beta, lambda]
 */
class AllInd0TrParam2DegPolyHMM : public AllInd3TrParam2DegPolyHMM {
  private:
  protected:
    // protected ctors
    AllInd0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    AllInd0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst);
    virtual void makeHMMs();
    virtual void makeHMMs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams);

    //// methods for least squares for solving for Baum Welch estimates. // TODO not impl here since don't think it makes sense to run baum welch on classes that don't allow all variables to vary?
    //virtual void baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    //virtual void baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    //virtual double baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) override;
    //virtual void saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) override;

  public:
    // constants

    // constructors and destructor
    virtual ~AllInd0TrParam2DegPolyHMM();
    static AllInd0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy);
    // accessors and mutators
    //virtual void setLibScalingFactor(int cellIdx, double lib);
    //virtual void setAllLibScalingFactors(double libScalingFactor);
    //virtual double getLibScalingFactor(int cellIdx) const;
    //virtual double setAllTransition(gsl_vector* transitionParams) override;
    virtual double setAllBranches(gsl_vector* branches);

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;
    //virtual void miscFunctions() override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    AllInd0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true);

};

#endif


