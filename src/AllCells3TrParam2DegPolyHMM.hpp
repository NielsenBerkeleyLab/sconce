#ifndef ALLCELLS3TRPARAM2DEGPOLYHMM_HPP
#define ALLCELLS3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <vector>
#include <unordered_map>

#include "Optimizable.hpp"
#include "DepthPair.hpp"
#include "HMM.hpp"

/*
 * This class is for making a collection of HMMs that jointly estimates beta and gammma across a given number of
 * tumor cells, where each HMM is independent of the other pairs (except for the shared parameters
 * and cell specific library size scaling factor (this is used every time a cell is used). Branch
 * lengths (t1/t2/t3/t/splitTime) are estimated for each HMM.
 *
 * Parameters are:
 *   - //3 transition parameters (alpha=P(adjacent CNA); beta=P(any CNA); gamma=P(back to diploid))
 *   - 3 rate parameters (alpha=P(adjacent CNA); beta=P(any CNA); lambda=P(event happens to both cells))
 *   - a library size scaling factor for each cell
 *   - branch lengths (3 branches between 2 cells, 1 scaling "branch" for 1 cell, or split time between 2 cells (given total branch lengths))
 *
 * This class is abstract, and is the super class to AllPairs* and AllInd*. Parameter layouts are determined in subclasses
 */
class AllCells3TrParam2DegPolyHMM : public Optimizable {
  private:
  protected:
    // protected ctors
    //AllCells3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy);
    //AllCells3TrParam2DegPolyHMM(const AllCells3TrParam2DegPolyHMM& otherAllCells3TrParam2DegPolyHMM);
    AllCells3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst);

    // member variables
    std::vector<DepthPair*>* depthsVec;
    std::vector<std::string>* sampleList; // list of cell/sample names as determined by filename (corresponds to depthsVec)
    std::vector<HMM*>* hmmVec;
    std::vector<std::string>* hmmNames; // list of cell/sample names as determined by filename (corresponds to hmmVec, comma sep, must be set in subclass)
    std::vector<gsl_vector*>* baumWelchParamResults; // vector of paramsToEst from baum welch + least squares, starting from different library scaling factors, stored in priority order for bfgs to explore

    // methods for least squares for solving for Baum Welch estimates
    virtual double doLeastSquares(gsl_vector* initGuess, int maxNumAttempts = 3);
    static double baumWelchLeastSquares_f(const gsl_vector* v, void* params);
    static void baumWelchLeastSquares_df(const gsl_vector* v, void* params, gsl_vector* df);
    static void baumWelchLeastSquares_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df);
    virtual void baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const = 0;
    virtual void baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const = 0;
    virtual double baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) = 0; // internal method to actually calculate the sum of squared residuals for least squares
    virtual void saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) = 0;

  public:
    // constants
    const int MAX_PLOIDY;
    const int NUM_CELLS;
    const int NUM_HMMS;
    const int NUM_LIBS_TO_EST;
    const int NUM_SHARED_TRANSITION_PARAMS_TO_EST; // ie. beta, lambda. Does not count branch lengths or lib size scaling factors
    const int NUM_BRANCH_LENGTHS_TO_EST;
    const int NUM_FIXED_LIBS;
    const int NUM_FIXED_SHARED_TRANSITION_PARAMS;
    const int LIB_SIZE_SCALING_FACTOR_START_IDX;
    const int SHARED_TRANSITION_PROB_START_IDX; // estimated params only
    const int BRANCH_LENGTH_START_IDX;
    const int FIXED_TRANSITION_PROB_START_IDX;

    // constructors and destructor
    virtual ~AllCells3TrParam2DegPolyHMM();

    // accessors and mutators
    int getKploidy() const;
    std::vector<HMM*>* getHMMs();
    virtual void print(FILE* stream) override;
    std::vector<std::string>* getSampleList();
    void setSampleList(std::vector<std::string>* sampleList);
    virtual std::vector<gsl_vector*>* getBaumWelchParamResults() const;
    virtual void setBaumWelchParamResults(std::vector<gsl_vector*>* paramResults);
    /*virtual void setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1);
    virtual void setAllLibScalingFactors(int cellNumInPair, double libScalingFactor);*/
    virtual double getLibScalingFactor(int cellNum) const = 0;
    virtual double getAlpha() const;
    virtual void setAlpha(double alpha);
    virtual double setAllTransition(gsl_vector* transitionParams) = 0;
    virtual void setAllMeanVarianceFn(gsl_vector* meanVarianceCoefVec);

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override = 0;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override = 0;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override = 0;
    virtual double getLogLikelihood() override;
    virtual double getViterbiLogLikelihood();
    virtual double checkOptimProbValidity(gsl_vector* probs) const override;
    virtual double checkStateValidity(double epislon = 1e-8) const override; // should check if transition matrices are ok
    virtual double checkForTransientStates(); // check if states only show up transiently
    virtual void setSimParamsToEst(gsl_vector* params) override = 0;
    virtual void setSimFixedParams(gsl_vector* params) override = 0;
    virtual void miscFunctions() override = 0;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override = 0;
    virtual void setBaumWelchInitGuess(gsl_vector* initGuess, int numBWIters, int numLibStarts = 3, double libStartValue = 1.0);
    AllCells3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true);
    double runForwardAlg();

    // methods to save results
    void viterbiDecodeAll();
    void saveAllViterbiDecodedCNA(std::string filename);
};

#endif

