#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <boost/program_options.hpp>
#include "HMM.hpp"
#include "DepthPair.hpp"
#include "util.hpp"
#include "Optimizable.hpp"

#include "AllInd2Stages3TrParam2DegPolyHMM.hpp"

namespace po = boost::program_options;

/*
 * Tue 16 Feb 2021 01:35:17 PM PST
 * program to independently analyze all cells passed in (as opposed to pairs)
 * if you pass in a list of tumor cells, all transition parameters (beta/lambda) will be shared among them
 * if you pass in one tumor cell filename, it'll have its own transition parameters (this is how to make completely indepdendent runs)
 */

int main(int argc, char** argv) {
  // DEBUGGING: TURN OFF BUFFERING. from https://stackoverflow.com/a/1716621
  setbuf(stdout, NULL);

  // parse command line options
  int maxKploid = 0;
  std::string diploidFile;
  std::string tumorFilename;
  std::string outputBase;
  bool verbose = false;
  bool fixLib = false; // should bfgs estimate libs, or should those be fixed from BW?
  bool estTrParamsBFGS = true; // should transition params be estimated with BFGS?
  bool disableEstTrParamsBFGS = false; // command line flag for disabling estimating transition params using BFGS

  // optimization start points
  bool simStart = false;
  bool bwStart = true; // run baum welch first
  bool disableBWstart = false; // command line flag for disabling running baum welch first
  int numBWIters = 0;
  int numLibStarts = 0;
  double libStartVal = 0;

  // config file options if want to specify simulation params to compare likelihood
  // see https://github.com/mousebird/boost/blob/master/libs/program_options/example/multiple_sources.cpp
  std::string simParamFile;
  double alpha = 0;
  double beta = 0;
  double lambda = 0;
  double t = 0;
  double lib = 0;
  bool haveSimParams = false;

  // config file options if want to specify mean/var coefs for negative binomial. Otherwise, uses default from HMM::createMeanVarianceCoefVec()
  std::string meanVarCoefFile;
  double intercept = 0;
  double slope = 0;
  double poly2 = 0;
  gsl_vector* meanVarianceCoefVec = nullptr;
  try {
    po::options_description cmdLineOps("Command line options");
    cmdLineOps.add_options()
      ("help,h", "help message")
      ("diploid,d", po::value<std::string>(&diploidFile)->required(), "path to diploid depth file")
      ("tumorFilename,s", po::value<std::string>(&tumorFilename), "path to single tumor depth file (ie to do fully independent runs)")
      ("maxKploid,k", po::value<int>(&maxKploid)->default_value(10), "maximum allowed ploidy")
      ("fixLib,f", po::bool_switch(&fixLib)->default_value(false), "fix library sizes after baum welch estimation")
      ("verbose,v", po::bool_switch(&verbose)->default_value(false), "enable debugging statements")
      ("simParamFile,m", po::value<std::string>(&simParamFile), "path to simulation parameter file")
      ("meanVarCoefFile", po::value<std::string>(&meanVarCoefFile), "path to negative binomial mean/variance coefficients file")
      ("sim", po::bool_switch(&simStart), "start at simulation params")
      ("disableBWstart", po::bool_switch(&disableBWstart)->default_value(false), "disable running baum welch to get initial starting points for BFGS")
      ("disableEstTrParamsBFGS,e", po::bool_switch(&disableEstTrParamsBFGS)->default_value(false), "disable estimating transition params using BFGS")
      ("bwIters", po::value<int>(&numBWIters)->default_value(20), "number of baum welch iterations")
      ("numLibStarts", po::value<int>(&numLibStarts)->default_value(3), "number of library starting points for baum welch (should be 1 (for [1]) or 3 (for [0.5, 1, 2]))")
      ("libStartVal", po::value<double>(&libStartVal)->default_value(1.0), "if numLibStarts == 1, the value to start the library size scaling factor at")
      ("outputBase,o", po::value<std::string>(&outputBase)->required(), "path to output files");

    po::options_description simParamFileOps("Simulation Parameter configuration file options");
    simParamFileOps.add_options()
      ("alpha,a", po::value<double>(&alpha), "alpha (P(adj CNA))")
      ("beta,b", po::value<double>(&beta)->required(), "beta (P(any CNA))")
      ("lambda,L", po::value<double>(&lambda)->required(), "lambda (rate that events affect both bins)")
      ("t", po::value<double>(&t)->required(), "t branch length")
      ("lib", po::value<double>(&lib)->required(), "lib (cell library size scaling factor)");

    po::options_description meanVarCoefFileOps("Negative Binomial Mean and Variance Coefficient configuration file options");
    meanVarCoefFileOps.add_options()
      ("intercept", po::value<double>(&intercept)->required(), "var = {intercept} + slope * mean + poly2 * mean^2")
      ("slope", po::value<double>(&slope)->required(), "var = intercept + {slope} * mean + poly2 * mean^2")
      ("poly2", po::value<double>(&poly2)->required(), "var = intercept + slope * mean + {poly2} * mean^2");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdLineOps), vm);

    if(vm.count("help") || argc == 1) {
      std::cout << cmdLineOps << std::endl;;
      std::cout << simParamFileOps << std::endl;
      std::cout << meanVarCoefFileOps << std::endl;;
      return 0;
    }
    po::notify(vm); // deal with any command line errors

    if(vm.count("simParamFile")) {
      std::ifstream cfg(simParamFile.c_str());
      if(!cfg) {
        std::cerr << "Error: cannot open " << simParamFile << std::endl;
        exit(EXIT_FAILURE);
      }
      po::store(po::parse_config_file(cfg, simParamFileOps, true), vm); // true is to allow unknown options (ie for alpha to be in the masterParams file without messing up param parsing; see https://stackoverflow.com/a/31922646
      haveSimParams = true;
      cfg.close();
    }
    po::notify(vm);

    if(vm.count("meanVarCoefFile")) {
      std::ifstream cfg(meanVarCoefFile.c_str());
      if(!cfg) {
        std::cerr << "Error: cannot open " << meanVarCoefFile << std::endl;
        exit(EXIT_FAILURE);
      }
      po::store(po::parse_config_file(cfg, meanVarCoefFileOps, true), vm); // true is to allow unknown options (ie for alpha to be in the masterParams file without messing up param parsing; see https://stackoverflow.com/a/31922646
      cfg.close();
      meanVarianceCoefVec = gsl_vector_alloc(3);
    }
    po::notify(vm); // store params from file into local vars, then store into meanVarianceCoefVec
    if(meanVarianceCoefVec != nullptr) {
      gsl_vector_set(meanVarianceCoefVec, 0, intercept);
      gsl_vector_set(meanVarianceCoefVec, 1, slope);
      gsl_vector_set(meanVarianceCoefVec, 2, poly2);
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
  if(disableEstTrParamsBFGS) {
    estTrParamsBFGS = false;
  }
  if(disableBWstart) {
    bwStart = false;
  }

  // ##### set up AllInd2Stages obj #####
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // make a DepthPair for the passed tumor file
  std::vector<DepthPair*>* depthsVec = new std::vector<DepthPair*>();
  DepthPair* firstDepthPair = nullptr;
  std::vector<std::string>* sampleList = new std::vector<std::string>();
  if(!boost::filesystem::exists(diploidFile)) {
    std::cerr << "Error: " << diploidFile << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!boost::filesystem::exists(tumorFilename)) {
    std::cerr << "Error: " << tumorFilename << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }
  firstDepthPair = new DepthPair(diploidFile, tumorFilename);
  depthsVec->push_back(firstDepthPair);
  sampleList->push_back(parseSampleName(tumorFilename));
  int numCells = depthsVec->size();
  // ##### end file reading #####

  // create all HMMs of tumor cells
  AllInd2Stages3TrParam2DegPolyHMM* allIndHMM = new AllInd2Stages3TrParam2DegPolyHMM(depthsVec, sampleList, maxKploid, verbose);
  double bfgsInitLl = 0;
  double optimLL = 0;

  // ##### run baum welch #####
  gsl_vector* bwInitGuess = nullptr;
  if(bwStart) {
    std::cout << "######## STARTING BAUM WELCH SET UP ########" << std::endl;
    allIndHMM->setUpBaumWelch();
    bwInitGuess = gsl_vector_alloc(allIndHMM->getBWInitGuessSize());
    AllInd3TrParam2DegPolyHMM* bwHMM = allIndHMM->getBWAllInd();
    if(meanVarianceCoefVec != nullptr) {
      bwHMM->setAllMeanVarianceFn(meanVarianceCoefVec);
    }
    if(haveSimParams) {
      gsl_vector* simParamsForBw = gsl_vector_alloc(bwHMM->getNumParamsToEst());
      gsl_vector* fixedParamsForBw = gsl_vector_alloc(bwHMM->getNumFixedParams());
      int simParamsIdx = 0;
      // libs
      for(; simParamsIdx < numCells; simParamsIdx++) {
        gsl_vector_set(simParamsForBw, simParamsIdx, lib);
      }

      // beta/lambda
      gsl_vector_set(simParamsForBw, simParamsIdx, beta);
      simParamsIdx++;
      gsl_vector_set(simParamsForBw, simParamsIdx, lambda);
      simParamsIdx++;

      // branch lengths
      for(int i = 0; i < numCells; i++) {
        gsl_vector_set(simParamsForBw, simParamsIdx + i, t);
      }
      gsl_vector_set(fixedParamsForBw, fixedParamsForBw->size - 1, alpha);
      bwHMM->setSimParamsToEst(simParamsForBw);
      bwHMM->setSimFixedParams(fixedParamsForBw);
    }
    allIndHMM->setBaumWelchInitGuess(bwInitGuess, numBWIters, numLibStarts, libStartVal);

    std::cout << "######## DONE WITH BAUM WELCH INITIALIZATION. BWINITGUESS IS ########" << std::endl;
    
    // INSERT bwInitGuess DEBUGGING SHORTCUTS HERE

    printColVector(bwInitGuess);
    allIndHMM->print(stdout);
  }
  // ##### end baum welch #####

  gsl_vector* fixedParams = nullptr;
  gsl_vector* paramsToEstStartFromBW = nullptr;
  int fixedIdx = 0;
  unsigned int bwInitGuessIdx = 0;
  // if estimate transition params (beta/lambda) with bfgs
  if(estTrParamsBFGS) {
    // assumes bw libs are available 
    if(fixLib) {
      fixedParams = gsl_vector_alloc(numCells + 1); // all libs + alpha are fixed
      gsl_vector_set_zero(fixedParams);

      // lib scaling factors
      for(; bwInitGuessIdx < depthsVec->size(); fixedIdx++, bwInitGuessIdx++) {
        gsl_vector_set(fixedParams, fixedIdx, gsl_vector_get(bwInitGuess, bwInitGuessIdx));
      }
      gsl_vector_set(fixedParams, fixedIdx + 0, allIndHMM->getBWAllInd()->getAlpha()); // alpha

      // save transition parameters + tree branch lengths from bw for bfgs
      paramsToEstStartFromBW = gsl_vector_alloc(2 + numCells); // beta/lambda + one branch per cell
      gsl_vector_set(paramsToEstStartFromBW, 0, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 0)); // beta
      gsl_vector_set(paramsToEstStartFromBW, 1, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 1)); // lambda
      for(int cellIdx = 0; cellIdx < numCells; cellIdx++) {
        gsl_vector_set(paramsToEstStartFromBW, 2 + cellIdx, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 2 + cellIdx)); // t
      }
    }
    // estimating everything with bfgs, no fixedParams (alpha is set automatically in ctor)
    else {
      fixedParams = nullptr;

      // save libs + transition parameters + tree branch lengths from bw for bfgs (ie everything)
      paramsToEstStartFromBW = gsl_vector_alloc(numCells + 2 + numCells); // libs + beta/lambda + one branch per cell
      gsl_vector_memcpy(paramsToEstStartFromBW, bwInitGuess);
    }

  }
  // else using transition param estimates from baum welch
  else {
    if(fixLib) {
      fixedParams = gsl_vector_alloc(numCells + 3); // all libs + alpha/beta/lambda are fixed
      gsl_vector_set_zero(fixedParams);

      // lib scaling factors
      for(; bwInitGuessIdx < (unsigned int) numCells; fixedIdx++, bwInitGuessIdx++) {
        gsl_vector_set(fixedParams, fixedIdx, gsl_vector_get(bwInitGuess, bwInitGuessIdx));
      }

      // save branch length ests into paramsToEstStartFromBW
      paramsToEstStartFromBW = gsl_vector_alloc(numCells); // one branch per cell
    }
    else {
      fixedParams = gsl_vector_alloc(3); // alpha/beta/lambda are fixed
      gsl_vector_set_zero(fixedParams);
      bwInitGuessIdx = numCells; // skip over libs

      // save lib size ests and branch length ests into paramsToEstStartFromBW
      paramsToEstStartFromBW = gsl_vector_alloc(numCells + numCells); // one lib + one branch per cell
      for(; bwInitGuessIdx < (unsigned int) numCells; fixedIdx++, bwInitGuessIdx++) {
        gsl_vector_set(paramsToEstStartFromBW, fixedIdx, gsl_vector_get(bwInitGuess, bwInitGuessIdx));
      }
    }

    // shared transition parameters into fixedParams
    gsl_vector_set(fixedParams, fixedIdx + 0, allIndHMM->getBWAllInd()->getAlpha()); // alpha
    gsl_vector_set(fixedParams, fixedIdx + 1, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 0)); // beta
    gsl_vector_set(fixedParams, fixedIdx + 2, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 1)); // lambda

    // save branch length ests into paramsToEstStartFromBW
    for(int cellIdx = 0; cellIdx < numCells; cellIdx++) {
      gsl_vector_set(paramsToEstStartFromBW, cellIdx, gsl_vector_get(bwInitGuess, bwInitGuessIdx + 2 + cellIdx)); // t, skip over beta/lambda in indexing
    }

  }
  if(fixedParams != nullptr) {
    std::cout << "fixedParams:" << std::endl;
    printColVector(fixedParams);
  }
  allIndHMM->setUpAllIndBFGS(fixLib, estTrParamsBFGS, fixedParams);
  AllInd3TrParam2DegPolyHMM* bfgsAllInd = allIndHMM->getBFGSAllInd();
  if(paramsToEstStartFromBW != nullptr) {
    std::cout << "paramsToEstStartFromBW:" << std::endl;
    printColVector(paramsToEstStartFromBW);
    bool bwParamAdjusted = false;
    for(unsigned int paramIdx = 0; paramIdx < paramsToEstStartFromBW->size; paramIdx++) {
      double bwParamValue = gsl_vector_get(paramsToEstStartFromBW, paramIdx);
      if(bwParamValue < 1e-4) {
        bwParamAdjusted = true;
        bwParamValue = 1e-3;
        gsl_vector_set(paramsToEstStartFromBW, paramIdx, bwParamValue);
      }
    }
    if(bwParamAdjusted) {
      std::cout << "adjusted paramsToEstStartFromBW values < 1e-4 to be 1e-3" << std::endl;
      printColVector(paramsToEstStartFromBW);
    }
    bfgsAllInd->setParamsToEst(paramsToEstStartFromBW);
  }
  if(meanVarianceCoefVec != nullptr) {
    bfgsAllInd->setAllMeanVarianceFn(meanVarianceCoefVec);
  }

  // if have bw lib ests, and reestimating them with bfgs, start bfgs with the bw ests
  if(bwStart && !fixLib) {
    double lib = 0;
    for(int cellIdx = 0; cellIdx < bfgsAllInd->NUM_LIBS_TO_EST; cellIdx++) {
      lib = gsl_vector_get(bwInitGuess, cellIdx);
      bfgsAllInd->setLibScalingFactor(cellIdx, lib);
    }
  }


  //bfgsAllInd->setLibScalingFactor(0, 1.3912459941748089864432813556049950420856); // debugging

  std::cout << "######## STARTING INDCELLS STAGE 2 BFGS ########" << std::endl;
  bfgsAllInd->print(stdout);
  bfgsInitLl = bfgsAllInd->getLogLikelihood();
  gsl_vector* tmpInitGuess = gsl_vector_alloc(bfgsAllInd->getNumParamsToEst());
  gsl_vector_memcpy(tmpInitGuess, bfgsAllInd->getParamsToEst());

  // ##### run bfgs on indCells stage 2 #####
  allIndHMM->callBFGSNTimes(tmpInitGuess, 1, verbose);
  std::cout << "######## DONE WITH INDCELLS STAGE 2 BFGS ########" << std::endl;
  bfgsAllInd->print(stdout);
  optimLL = allIndHMM->getBFGSAllIndLogLikelihood();


  std::chrono::steady_clock::time_point decodeStart = std::chrono::steady_clock::now();
  allIndHMM->viterbiDecode();
  allIndHMM->saveViterbiDecodedCNA(outputBase);
  std::chrono::steady_clock::time_point decodeEnd = std::chrono::steady_clock::now();
  double decodeElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(decodeEnd - decodeStart).count() / 1000000.0;
  printf("TOTAL DECODING TIME ELAPSED (sec): %.5f\n", decodeElapsedSec);
  printf("AVERAGE DECODING TIME ELAPSED (sec): %.5f\n", decodeElapsedSec / (double) numCells);

  // #### save the optimized HMM ####
  std::cout << "######## SAVING HMM ########" << std::endl;
  FILE* oFile = fopen((outputBase + ".hmm").c_str(), "w");
  allIndHMM->print(oFile);
  fclose(oFile);

  // print optimized params for this param set
  printf("OPTIMIZED FINAL LOGLIKELIHOOD: %.40f\n", optimLL);
  gsl_vector* bestOptimLLInitGuess = bfgsAllInd->getBestLLInitGuess();
  if(bestOptimLLInitGuess != nullptr) {
    gsl_vector* finalParams = gsl_vector_alloc(bfgsAllInd->getNumParamsToEst());
    gsl_vector_memcpy(finalParams, bfgsAllInd->getParamsToEst());
    bfgsAllInd->setParamsToEst(bestOptimLLInitGuess);
    bfgsInitLl = bfgsAllInd->getLogLikelihood();
    bfgsAllInd->setParamsToEst(finalParams);
    gsl_vector_free(finalParams);
  }
  printf("DIFFERENCE IN LOGLIKELIHOOD: %.40f\n", optimLL - bfgsInitLl);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
  printf("TOTAL TIME ELAPSED (sec): %.5f\n", elapsedSec);

  // clean up
  // TODO
  
  return 0;
}

