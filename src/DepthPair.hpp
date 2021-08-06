#ifndef DEPTHFILE_HPP
#define DEPTHFILE_HPP

#include <string>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

#include <gsl/gsl_matrix.h>

#include <boost/filesystem.hpp> // for checking if files exist

#include "util.hpp" // for printMatrix

class DepthPair {
  public:
    // constructors and destructor
    DepthPair(std::string diploidFilename, std::string tumorFilename);
    DepthPair(DepthPair* otherDepths, std::string tumorFilename);
    DepthPair(int numWindows, int numChr = 1, int windowSize = 250000); // ctor for simulation
    DepthPair(bool shareDiploid, DepthPair* otherDepths); // ctor for simulation for subsequent cells
    DepthPair(std::string windowsFilename); // ctor for simulation with a windows file (ie if want to match a reference genome)
    DepthPair(const DepthPair& other);
    ~DepthPair();

    // member variables
    double diploidLibrarySize;
    double tumorLibrarySize;
    unsigned int numWindows;
    int maxWindowSize;
    double maxDiploidDepth;
    double maxTumorDepth;
    std::vector<std::string>* chrVec; // list of all chromosomes
    std::vector<std::string>* allRegions; // list of all regions
    std::unordered_map<std::string, std::vector<std::string>*>* regions; // chr:vector of region
    std::unordered_map<std::string, std::vector<double>*>* chrToDiploidDepthMap; // chr:vector of depths
    std::unordered_map<std::string, std::vector<double>*>* chrToDiploidVarMap; // chr:vector of depths
    std::unordered_map<std::string, std::vector<double>*>* chrToTumorDepthMap; // chr:vector of depths
    std::unordered_map<std::string, std::vector<std::string>*>* chrToDiploidSimStateMap; // chr:vector of simulated states
    std::unordered_map<std::string, std::vector<std::string>*>* chrToTumorSimStateMap; // chr:vector of simulated states
    std::unordered_map<std::string, std::vector<int>*>* chrToViterbiPathMap; // chr:vector of viterbi decoded paths; assigned after viterbiDecode() is called in the HMM
    std::unordered_map<std::string, gsl_matrix*>* forBackMargMatMap; // chr:marginalized likelihood matrix from forBackAlg
    std::vector<double>* avgDiploidDepth;
    std::vector<double>* avgTumorDepth;

    // functions
    void print(FILE* stream);
    double getDiploidDepthAt(std::string chr, int i);
    double getTumorDepthAt(std::string chr, int i);
    std::vector<double>* getAverageDepth(bool getDiploid);
    double getTotalDiploidDepth();
    double getTotalTumorDepth();
    void saveSimDiploidAsDepthFile(std::string filename);
    void saveSimTumorAsDepthFile(std::string filename);
};

#endif

