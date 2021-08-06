#ifndef ALLIND2STAGES3TRPARAM2DEGPOLYHMM_HPP
#define ALLIND2STAGES3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include "Optimizable.hpp"
#include "DepthPair.hpp"

#include "AllCells3TrParam2DegPolyHMM.hpp"

#include "AllInd3TrParam2DegPolyHMM.hpp"
#include "AllInd0TrParam2DegPolyHMM.hpp"
#include "AllIndFixLib3TrParam2DegPolyHMM.hpp"
#include "AllIndFixLib0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that independently process cells.
 * alpha is fixed, and beta and lambda are shared between pairs.
 * Library sizes are independently estimated, as are branch lengths.
 *
 * Beta and lambda are optionally estimated via baum welch + least squares,
 * as are library sizes can be estimated with baum welch
 *
 * Then branch lengths (and optionally shared transition params and library sizes) are
 * estimated through BFGS.
 * As of Tue 16 Feb 2021 02:50:35 PM PST, unclear if 2 stages is necessary for time saving reasons
 *
 * All paramter storage is delegated to its member variables.
 */
class AllInd2Stages3TrParam2DegPolyHMM {
  private:
    // member variables
    AllInd3TrParam2DegPolyHMM* bwAllInd;
    AllInd3TrParam2DegPolyHMM* bfgsAllInd;
    std::vector<DepthPair*>* depthsVec;
    std::vector<std::string>* sampleList; // list of cell/sample names as determined by filename (corresponds to depthsVec)
    bool gradientDebug;

  public:
    // constants
    const int NUM_CELLS;
    const int NUM_SHARED_TRANSITION_PARAMS_TO_EST; // ie. beta, lambda. Does not count branch lengths or lib size scaling factors
    const int MAX_PLOIDY;

    // constructors and destructor
    AllInd2Stages3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, bool gradientDebug);
    virtual ~AllInd2Stages3TrParam2DegPolyHMM();

    // accessors and mutators
    int getKploidy() const;
    void setUpBaumWelch();
    void setUpAllIndBFGS(bool fixLib, bool estTrParamsBFGS, gsl_vector* fixedParams);
    AllInd3TrParam2DegPolyHMM* getBWAllInd();
    AllInd3TrParam2DegPolyHMM* getBFGSAllInd();
    void print(FILE* stream);
    std::vector<std::string>* getSampleList();
    void setSampleList(std::vector<std::string>* sampleList);
    int getInitGuessSize() const;
    int getBWInitGuessSize() const;

    double getBWAllIndLogLikelihood();
    double getBFGSAllIndLogLikelihood();

    AllInd3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, bool verbose = true);
    AllInd3TrParam2DegPolyHMM* callBFGSNTimes(gsl_vector* initGuess, int numRuns, bool verbose = true);
    void setInitGuessNthTime(gsl_vector* initGuess, int iter) const;
    void setSharedInitGuess(gsl_vector* initGuess, double lib, double beta, double lambda, double t) const;
    void setBaumWelchInitGuess(gsl_vector* initGuess, int numBWIters, int numLibStarts = 3, double libStartVal = 1.0); // runs baum welch to find a good starting place for bfgs, saves into passed initGuess
    void copyBaumWelchIntoBFGSInitGuess(gsl_vector* bfgsInitGuess, gsl_vector* bwInitGuess);

    // methods to save results
    void viterbiDecode();
    void saveViterbiDecodedCNA(std::string filename);

};

#endif

