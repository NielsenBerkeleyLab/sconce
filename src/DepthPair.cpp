#include "DepthPair.hpp"

DepthPair::DepthPair(std::string diploidFilename, std::string tumorFilename) {
  this->diploidLibrarySize = 0;
  this->tumorLibrarySize = 0;
  this->maxWindowSize = -1;
  this->maxDiploidDepth = 0;
  this->maxTumorDepth = 0;
  this->chrVec = new std::vector<std::string>(); // to maintain insertion order of chromosomes
  this->allRegions = new std::vector<std::string>(); // to maintain insertion order of keys (regions) across genome
  this->regions = new std::unordered_map<std::string, std::vector<std::string>*>(); // chr:region; to maintain insertion order of keys (regions) per chr
  this->chrToDiploidDepthMap = new std::unordered_map<std::string, std::vector<double>*>(); // chr:vector of depths
  this->chrToDiploidVarMap = new std::unordered_map<std::string, std::vector<double>*>(); // chr:vector of variances
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>(); // chr:vector of depths
  this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>(); // chr:vector of simulated diploid states
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>(); // chr:vector simulated tumor states

  if(!boost::filesystem::exists(diploidFilename)) {
    std::cerr << "Error: " << diploidFilename << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!boost::filesystem::exists(tumorFilename)) {
    std::cerr << "Error: " << tumorFilename << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream diploid(diploidFilename);
  std::ifstream tumor(tumorFilename);
  std::string diploidLine;
  std::string tumorLine;
  std::string key;
  std::string testKey;
  std::string chr;
  long long int start;
  long long int end;
  double depth; // number of reads that overlap this region
  double variance;
  while(getline(diploid, diploidLine)) {
    // process diploid line
    std::istringstream iss(diploidLine); // TODO make more memory efficient?
    iss >> chr >> start >> end >> depth >> variance;
    key = chr + ":" + std::to_string(start) + "-" + std::to_string(end);

    // if chr doesn't have a depth vector yet, make one
    std::unordered_map<std::string, std::vector<double>*>::iterator chrRegionItr = this->chrToDiploidDepthMap->find(chr);
    if(chrRegionItr == this->chrToDiploidDepthMap->end()) {
      this->chrVec->push_back(chr);
      this->chrToDiploidDepthMap->insert(make_pair(chr, new std::vector<double>()));
      this->chrToDiploidVarMap->insert(make_pair(chr, new std::vector<double>()));
      this->regions->insert(make_pair(chr, new std::vector<std::string>()));
    }

    // insert into appropriate depth vector
    (*this->chrToDiploidDepthMap)[chr]->push_back(depth);
    (*this->chrToDiploidVarMap)[chr]->push_back(variance);

    (*this->regions)[chr]->push_back(key);
    this->diploidLibrarySize += depth;
    if(this->maxDiploidDepth < depth) {
      this->maxDiploidDepth = depth;
    }
    this->allRegions->push_back(key);

    // check if should update maxWindowSize
    if(this->maxWindowSize < (end - start)) {
      this->maxWindowSize = end - start;
    }
  }

  //std::string tmpSimState; // Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
  std::vector<std::string>::iterator it = this->allRegions->begin();
  while(getline(tumor, tumorLine) && it != this->allRegions->end()) {
    // process tumor line
    std::istringstream iss(tumorLine);
    iss >> chr >> start >> end >> depth;
    //iss >> chr >> start >> end >> depth >> tmpSimState; // Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
    key = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    if(key.compare(*it) != 0) {
      std::cerr << "Error: Regions in depth files are not equal. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }
    this->tumorLibrarySize += depth;

    // if chr doesn't have a depth map yet, make one
    std::unordered_map<std::string, std::vector<double>*>::iterator chrRegionItr = this->chrToTumorDepthMap->find(chr);
    if(chrRegionItr == this->chrToTumorDepthMap->end()) {
      this->chrToTumorDepthMap->insert(make_pair(chr, new std::vector<double>()));
    }

    // insert into appropriate depth map
    (*this->chrToTumorDepthMap)[chr]->push_back(depth);

    /*// Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging. save into chrToTumorSimStateMap
    std::unordered_map<std::string, std::vector<std::string>*>::iterator chrSimStateItr = this->chrToTumorSimStateMap->find(chr);
    if(chrSimStateItr == this->chrToTumorSimStateMap->end()) {
      this->chrToTumorSimStateMap->insert(make_pair(chr, new std::vector<std::string>()));
    }
    //(*this->chrToTumorSimStateMap)[chr]->push_back(tmpSimState);
    std::string* tmpSimCopy = new std::string(tmpSimState);
    (*this->chrToTumorSimStateMap)[chr]->push_back(*tmpSimCopy);
    //std::cout << (*(*this->chrToTumorSimStateMap)[chr])[(*(*this->chrToTumorSimStateMap)[chr]).size() - 1] << std::endl;
    */



    if(this->maxTumorDepth < depth) {
      this->maxTumorDepth = depth;
    }
    ++it;
  }
  this->numWindows = this->allRegions->size();

  // get avg depths
  this->avgDiploidDepth = nullptr; // get rid of warning about conditional jump on uninit values
  this->avgTumorDepth = nullptr;
  this->avgDiploidDepth = this->getAverageDepth(true);
  this->avgTumorDepth = this->getAverageDepth(false);

  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>();
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>();

  // close filestreams
  diploid.close();
  tumor.close();
}

/*
 * constructor for the second tumor cell. Shares as much data from first DepthPair as possible
 */
DepthPair::DepthPair(DepthPair* otherDepths, std::string tumorFilename) {
  this->diploidLibrarySize = otherDepths->diploidLibrarySize;
  this->tumorLibrarySize = 0;
  this->numWindows = otherDepths->numWindows;
  this->maxWindowSize = otherDepths->maxWindowSize;
  this->maxDiploidDepth = otherDepths->maxDiploidDepth;
  this->maxTumorDepth = 0;
  this->chrVec = otherDepths->chrVec; // to maintain insertion order of chromosomes
  this->allRegions = otherDepths->allRegions; // to maintain insertion order of keys (regions) across genome
  this->regions = otherDepths->regions; // chr:region; to maintain insertion order of keys (regions) per chr
  this->chrToDiploidDepthMap = otherDepths->chrToDiploidDepthMap; // chr:vector of depths
  this->chrToDiploidVarMap = otherDepths->chrToDiploidVarMap; // chr:vector of variances
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>(); // chr:vector of depths
  this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->avgDiploidDepth = otherDepths->avgDiploidDepth;

  if(!boost::filesystem::exists(tumorFilename)) {
    std::cerr << "Error: " << tumorFilename << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream tumor(tumorFilename);
  std::string tumorLine;
  std::string key;
  std::string testKey;
  std::string chr;
  long long int start;
  long long int end;
  double depth; // number of reads that overlap this region

  //std::string tmpSimState; // Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
  std::vector<std::string>::iterator it = this->allRegions->begin();
  while(getline(tumor, tumorLine) && it != this->allRegions->end()) {
    // process tumor line
    std::istringstream iss(tumorLine);
    iss >> chr >> start >> end >> depth;
    //iss >> chr >> start >> end >> depth >> tmpSimState; // Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging
    key = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    if(key.compare(*it) != 0) {
      std::cerr << "Error: Regions in depth files are not equal. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }
    this->tumorLibrarySize += depth;

    // if chr doesn't have a depth map yet, make one
    std::unordered_map<std::string, std::vector<double>*>::iterator chrRegionItr = this->chrToTumorDepthMap->find(chr);
    if(chrRegionItr == this->chrToTumorDepthMap->end()) {
      this->chrToTumorDepthMap->insert(make_pair(chr, new std::vector<double>()));
    }

    // insert into appropriate depth map
    (*this->chrToTumorDepthMap)[chr]->push_back(depth);

    /*// Tue 24 Mar 2020 05:53:24 PM PDT simStateEmi debugging. save into chrToTumorSimStateMap
    std::unordered_map<std::string, std::vector<std::string>*>::iterator chrSimStateItr = this->chrToTumorSimStateMap->find(chr);
    if(chrSimStateItr == this->chrToTumorSimStateMap->end()) {
      this->chrToTumorSimStateMap->insert(make_pair(chr, new std::vector<std::string>()));
    }
    //(*this->chrToTumorSimStateMap)[chr]->push_back(tmpSimState);
    std::string* tmpSimCopy = new std::string(tmpSimState);
    (*this->chrToTumorSimStateMap)[chr]->push_back(*tmpSimCopy);
    //std::cout << (*this->chrToTumorSimStateMap)[chr] << std::endl;*/


    if(this->maxTumorDepth < depth) {
      this->maxTumorDepth = depth;
    }
    ++it;
  }

  // get avg depths
  this->avgDiploidDepth = otherDepths->avgDiploidDepth;
  this->avgTumorDepth = nullptr; // must set to null first
  this->avgTumorDepth = this->getAverageDepth(false);

  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>();
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>();

  // close filestreams
  tumor.close();
}

// ctor for simulation
DepthPair::DepthPair(int numWindows, int numChr, int windowSize) {
  this->diploidLibrarySize = 0;
  this->tumorLibrarySize = 0;
  this->numWindows = numWindows;
  this->maxDiploidDepth = 0;
  this->maxTumorDepth = 0;
  this->chrVec = new std::vector<std::string>();
  this->allRegions = new std::vector<std::string>();
  int winPerChr = numWindows / numChr;
  this->regions = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->maxWindowSize = windowSize;

  this->chrToDiploidDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToDiploidVarMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>();
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>();
  this->avgDiploidDepth = new std::vector<double>();
  this->avgTumorDepth = new std::vector<double>();

  for(int i = 1; i <= numChr; i++) {
    std::string chr = "chr" + std::to_string(i);
    this->chrVec->push_back(chr);
    this->chrToDiploidDepthMap->insert(make_pair(chr, new std::vector<double>()));
    this->chrToDiploidVarMap->insert(make_pair(chr, new std::vector<double>()));
    this->chrToTumorDepthMap->insert(make_pair(chr, new std::vector<double>()));
    this->regions->insert(make_pair(chr, new std::vector<std::string>()));
    for(long long int j = 0; j < winPerChr; j++) {
      long long int start = j * windowSize + 1;
      long long int end = (j+1) * windowSize;
      // key is chr#:start-end
      std::string key = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
      //std::cerr << key << std::endl;
      this->allRegions->push_back(key);
      (*this->regions)[chr]->push_back(key);
      (*this->chrToDiploidDepthMap)[chr]->push_back(0);
      (*this->chrToDiploidVarMap)[chr]->push_back(0);
      (*this->chrToTumorDepthMap)[chr]->push_back(0);
    }
  }
}

// ctor for simulation for subsequent cells. shareDiploid should always be true, here as a safety measure to not accidentally call copy ctor
//DepthPair::DepthPair(DepthPair* otherDepths, bool shareDiploid) {
DepthPair::DepthPair(bool shareDiploid, DepthPair* otherDepths) {
  this->diploidLibrarySize = otherDepths->diploidLibrarySize;
  this->tumorLibrarySize = otherDepths->tumorLibrarySize;
  this->numWindows = otherDepths->numWindows;
  this->maxDiploidDepth = otherDepths->maxDiploidDepth;
  this->maxTumorDepth = otherDepths->maxTumorDepth;
  this->chrVec = otherDepths->chrVec;
  this->allRegions = otherDepths->allRegions;
  this->regions = otherDepths->regions;
  this->maxWindowSize = otherDepths->maxWindowSize;

  if(shareDiploid) {
    this->chrToDiploidDepthMap = otherDepths->chrToDiploidDepthMap;
    this->chrToDiploidVarMap = otherDepths->chrToDiploidVarMap;
    this->chrToDiploidSimStateMap = otherDepths->chrToDiploidSimStateMap;
    this->avgDiploidDepth = otherDepths->avgDiploidDepth;
  }
  else {
    this->chrToDiploidDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
    this->chrToDiploidVarMap = new std::unordered_map<std::string, std::vector<double>*>();
    this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
    this->avgDiploidDepth = new std::vector<double>();
  }
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->avgTumorDepth = new std::vector<double>();
  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>();
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>();

  std::string currChr;
  for(unsigned int i = 0; i < this->chrVec->size(); i++) {
    currChr = (*chrVec)[i];
    if(!shareDiploid) {
      this->chrToDiploidDepthMap->insert(make_pair(currChr, new std::vector<double>()));
      this->chrToDiploidVarMap->insert(make_pair(currChr, new std::vector<double>()));
    }
    this->chrToTumorDepthMap->insert(make_pair(currChr, new std::vector<double>()));
    for(unsigned int j = 0; j < (*this->regions)[currChr]->size(); j++) {
      if(!shareDiploid) {
        (*this->chrToDiploidDepthMap)[currChr]->push_back(0);
        (*this->chrToDiploidVarMap)[currChr]->push_back(0);
      }
      (*this->chrToTumorDepthMap)[currChr]->push_back(0);
    }
  }
}

// ctor for simulation with a windows file (ie if want to match a reference genome)
DepthPair::DepthPair(std::string windowsFilename) {
  this->diploidLibrarySize = 0;
  this->tumorLibrarySize = 0;
  this->maxDiploidDepth = 0;
  this->maxTumorDepth = 0;
  this->chrVec = new std::vector<std::string>();
  this->allRegions = new std::vector<std::string>();
  this->regions = new std::unordered_map<std::string, std::vector<std::string>*>();

  this->chrToDiploidDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToDiploidVarMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>();
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>();
  this->avgDiploidDepth = new std::vector<double>();
  this->avgTumorDepth = new std::vector<double>();

  std::ifstream windows(windowsFilename);
  std::string windowLine;
  std::string key;
  std::string chr;
  long long int start;
  long long int end;
  while(getline(windows, windowLine)) {
    std::istringstream iss(windowLine);
    iss >> chr >> start >> end;
    key = chr + ":" + std::to_string(start) + "-" + std::to_string(end);

    // if chr doesn't have a depth vector yet, make one
    std::unordered_map<std::string, std::vector<double>*>::iterator chrRegionItr = this->chrToDiploidDepthMap->find(chr);
    if(chrRegionItr == this->chrToDiploidDepthMap->end()) {
      this->chrVec->push_back(chr);
      this->chrToDiploidDepthMap->insert(make_pair(chr, new std::vector<double>()));
      this->chrToDiploidVarMap->insert(make_pair(chr, new std::vector<double>()));
      this->chrToTumorDepthMap->insert(make_pair(chr, new std::vector<double>()));
      this->regions->insert(make_pair(chr, new std::vector<std::string>()));
    }

    // insert filler 0's into appropriate depth vector
    (*this->chrToDiploidDepthMap)[chr]->push_back(0);
    (*this->chrToDiploidVarMap)[chr]->push_back(0);
    (*this->chrToTumorDepthMap)[chr]->push_back(0);

    (*this->regions)[chr]->push_back(key);
    this->allRegions->push_back(key);

    // check if should update maxWindowSize
    if(this->maxWindowSize < (end - start)) {
      this->maxWindowSize = end - start;
    }
  }
  this->numWindows = this->allRegions->size();

  // close filestreams
  windows.close();
}


// deep copy ctor
DepthPair::DepthPair(const DepthPair& other) {
  this->diploidLibrarySize = other.diploidLibrarySize;
  this->tumorLibrarySize = other.tumorLibrarySize;
  this->numWindows = other.numWindows;
  this->maxWindowSize = other.maxWindowSize;
  this->maxDiploidDepth = other.maxDiploidDepth;
  this->maxTumorDepth = other.maxTumorDepth;
  this->chrVec = new std::vector<std::string>(*other.chrVec);
  this->allRegions = new std::vector<std::string>(*other.allRegions);
  this->regions = new std::unordered_map<std::string, std::vector<std::string>*>();
  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = other.regions->begin(); it != other.regions->end(); ++it) {
    this->regions->insert(make_pair((*it).first, new std::vector<std::string>(*(*it).second)));
  }
  this->chrToDiploidDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = other.chrToDiploidDepthMap->begin(); it != other.chrToDiploidDepthMap->end(); ++it) {
    this->chrToDiploidDepthMap->insert(make_pair((*it).first, new std::vector<double>(*(*it).second)));
  }
  this->chrToDiploidVarMap = new std::unordered_map<std::string, std::vector<double>*>();
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = other.chrToDiploidVarMap->begin(); it != other.chrToDiploidVarMap->end(); ++it) {
    this->chrToDiploidVarMap->insert(make_pair((*it).first, new std::vector<double>(*(*it).second)));
  }
  this->chrToTumorDepthMap = new std::unordered_map<std::string, std::vector<double>*>();
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = other.chrToTumorDepthMap->begin(); it != other.chrToTumorDepthMap->end(); ++it) {
    this->chrToTumorDepthMap->insert(make_pair((*it).first, new std::vector<double>(*(*it).second)));
  }
  this->chrToDiploidSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = other.chrToDiploidSimStateMap->begin(); it != other.chrToDiploidSimStateMap->end(); ++it) {
    this->chrToDiploidSimStateMap->insert(make_pair((*it).first, new std::vector<std::string>(*(*it).second)));
  }
  this->chrToTumorSimStateMap = new std::unordered_map<std::string, std::vector<std::string>*>();
  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = other.chrToTumorSimStateMap->begin(); it != other.chrToTumorSimStateMap->end(); ++it) {
    this->chrToTumorSimStateMap->insert(make_pair((*it).first, new std::vector<std::string>(*(*it).second)));
  }
  this->chrToViterbiPathMap = new std::unordered_map<std::string, std::vector<int>*>(other.chrToViterbiPathMap->size());
  for(std::unordered_map<std::string, std::vector<int>*>::iterator it = other.chrToViterbiPathMap->begin(); it != other.chrToViterbiPathMap->end(); ++it) {
    this->chrToViterbiPathMap->insert(make_pair((*it).first, new std::vector<int>(*(*it).second)));
  }
  this->forBackMargMatMap = new std::unordered_map<std::string, gsl_matrix*>(other.forBackMargMatMap->size());
  for(std::unordered_map<std::string, gsl_matrix*>::iterator it = other.forBackMargMatMap->begin(); it != other.forBackMargMatMap->end(); ++it) {
    gsl_matrix* otherMat = (*it).second;
    gsl_matrix* copyMat = gsl_matrix_alloc(otherMat->size1, otherMat->size2);
    gsl_matrix_memcpy(copyMat, otherMat);
    this->forBackMargMatMap->insert(make_pair((*it).first, copyMat));
  }
  this->avgDiploidDepth = new std::vector<double>(*other.avgDiploidDepth);
  this->avgTumorDepth = new std::vector<double>(*other.avgTumorDepth);
}

// destructor
DepthPair::~DepthPair() {
  delete this->chrVec;
  delete this->allRegions;
  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = this->regions->begin(); it != this->regions->end(); ++it) {
    //delete (*it).second;
    //(*it).second->clear();
    delete (*it).second;
  }
  delete this->regions;
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = this->chrToDiploidDepthMap->begin(); it != this->chrToDiploidDepthMap->end(); ++it) {
    //std::cerr << "deleting in dip depth mat" << std::endl;
    //(*it).second->clear();
    delete (*it).second;
  }

  if(this->chrToDiploidDepthMap != nullptr) {
    delete this->chrToDiploidDepthMap;
  }
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = this->chrToDiploidVarMap->begin(); it != this->chrToDiploidVarMap->end(); ++it) {
    //(*it).second->clear();
    delete (*it).second;
  }
  if(this->chrToDiploidVarMap != nullptr) {
    delete this->chrToDiploidVarMap;
  }
  for(std::unordered_map<std::string, std::vector<double>*>::iterator it = this->chrToTumorDepthMap->begin(); it != this->chrToTumorDepthMap->end(); ++it) {
    //(*it).second->clear();
    delete (*it).second;
  }
  delete this->chrToTumorDepthMap;

  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = this->chrToDiploidSimStateMap->begin(); it != this->chrToDiploidSimStateMap->end(); ++it) {
    //(*it).second->clear();
    delete (*it).second;
  }
  delete this->chrToDiploidSimStateMap;
  for(std::unordered_map<std::string, std::vector<std::string>*>::iterator it = this->chrToTumorSimStateMap->begin(); it != this->chrToTumorSimStateMap->end(); ++it) {
    //(*it).second->clear();
    delete (*it).second;
  }
  delete this->chrToTumorSimStateMap;

  if(this->chrToViterbiPathMap != nullptr) {
    for(std::unordered_map<std::string, std::vector<int>*>::iterator it = this->chrToViterbiPathMap->begin(); it != this->chrToViterbiPathMap->end(); ++it) {
      //(*it).second->clear();
      delete (*it).second;
    }
    delete this->chrToViterbiPathMap;
  }

  if(this->forBackMargMatMap != nullptr) {
    for(std::unordered_map<std::string, gsl_matrix*>::iterator it = this->forBackMargMatMap->begin(); it != this->forBackMargMatMap->end(); ++it) {
      gsl_matrix_free((*it).second);
    }
    delete this->forBackMargMatMap;
  }

  if(this->avgDiploidDepth != nullptr) {
    delete this->avgDiploidDepth;
  }
  if(this->avgTumorDepth != nullptr) {
    delete this->avgTumorDepth;
  }
}

// function to print important fields from this object, should be called from HMM::print
void DepthPair::print(FILE* stream) {
  std::vector<std::string>::iterator chrItr = this->chrVec->begin();
  /*fprintf(stream, "observed diploid sequence:\n");
  std::vector<double>* diploidSeq = nullptr;
  for(; chrItr != this->chrVec->end(); ++chrItr) {
    diploidSeq = (*this->chrToDiploidDepthMap)[*chrItr];
    fprintf(stream, "%s:\n", (*chrItr).c_str());
    for(std::vector<double>::iterator it = diploidSeq->begin(); it != diploidSeq->end(); ++it) {
      fprintf(stream, "%f ", *it);
    }
    fprintf(stream, "\n\n");
  }
  fprintf(stream, "\n");*/

  /*fprintf(stream, "observed diploid variance:\n");
  chrItr = this->chrVec->begin();
  std::vector<double>* diploidVar = nullptr;
  for(; chrItr != this->chrVec->end(); ++chrItr) {
    diploidVar = (*this->chrToDiploidVarMap)[*chrItr];
    fprintf(stream, "%s:\n", (*chrItr).c_str());
    for(std::vector<double>::iterator it = diploidVar->begin(); it != diploidVar->end(); ++it) {
      fprintf(stream, "%f ", *it);
    }
    fprintf(stream, "\n\n");
  }
  fprintf(stream, "\n");*/

  /*fprintf(stream, "observed tumor sequence:\n");
  chrItr = this->chrVec->begin();
  std::vector<double>* tumorSeq = nullptr;
  for(; chrItr != this->chrVec->end(); ++chrItr) {
    tumorSeq = (*this->chrToTumorDepthMap)[*chrItr];
    fprintf(stream, "%s:\n", (*chrItr).c_str());
    for(std::vector<double>::iterator it = tumorSeq->begin(); it != tumorSeq->end(); ++it) {
      fprintf(stream, "%g ", *it);
    }
    fprintf(stream, "\n\n");
  }
  fprintf(stream, "\n");*/

  /*if(this->chrToDiploidSimStateMap->size() > 0) {
    fprintf(stream, "simulated diploid copy number sequence:\n");
    chrItr = this->chrVec->begin();
    std::vector<std::string>* simStates = nullptr;
    for(; chrItr != this->chrVec->end(); ++chrItr) {
      simStates = (*this->chrToDiploidSimStateMap)[*chrItr];
      fprintf(stream, "%s:\n", (*chrItr).c_str());
      for(std::vector<std::string>::iterator it = simStates->begin(); it != simStates->end(); ++it) {
        fprintf(stream, "%s ", (*it).c_str());
      }
      fprintf(stream, "\n\n");
    }
    fprintf(stream, "\n");
  }*/

  /*if(this->chrToTumorSimStateMap->size() > 0) {
    fprintf(stream, "simulated tumor copy number sequence:\n");
    chrItr = this->chrVec->begin();
    std::vector<std::string>* simStates = nullptr;
    for(; chrItr != this->chrVec->end(); ++chrItr) {
      simStates = (*this->chrToTumorSimStateMap)[*chrItr];
      fprintf(stream, "%s:\n", (*chrItr).c_str());
      for(std::vector<std::string>::iterator it = simStates->begin(); it != simStates->end(); ++it) {
        fprintf(stream, "%s ", (*it).c_str());
      }
      fprintf(stream, "\n\n");
    }
    fprintf(stream, "\n");
  }*/
  if(this->chrToViterbiPathMap->size() > 0) {
    fprintf(stream, "viterbi decoded copy number sequence:\n");
    chrItr = this->chrVec->begin();
    std::vector<int>* path = nullptr;
    for(; chrItr != this->chrVec->end(); ++chrItr) {
      path = (*this->chrToViterbiPathMap)[*chrItr];
      fprintf(stream, "%s:\n", (*chrItr).c_str());
      for(std::vector<int>::iterator it = path->begin(); it != path->end(); ++it) {
        fprintf(stream, "%i ", *it);
      }
      fprintf(stream, "\n\n");
    }
    fprintf(stream, "\n");
  }
  /*if(this->forBackMargMatMap->size() > 0) {
    fprintf(stream, "marginalized forwardBackward matrix:\n");
    chrItr = this->chrVec->begin();
    for(; chrItr != this->chrVec->end(); ++chrItr) {
      printMatrix(stream, (*this->forBackMargMatMap)[*chrItr]);
    }
  }*/
}

// function to return diploid depth (symbol) from ith window in the passed chr
double DepthPair::getDiploidDepthAt(std::string chr, int i) {
  return (*(*this->chrToDiploidDepthMap)[chr])[i];
}

// function to return tumor depth (symbol) from ith window in the passed chr
double DepthPair::getTumorDepthAt(std::string chr, int i) {
  return (*(*this->chrToTumorDepthMap)[chr])[i];
}

/*
 * helper function to get the average depth per chromosome
 */
std::vector<double>* DepthPair::getAverageDepth(bool getDiploid) {
  if(getDiploid && this->avgDiploidDepth != nullptr) {
    return this->avgDiploidDepth;
  }
  if(!getDiploid && this->avgTumorDepth != nullptr) {
    return this->avgTumorDepth;
  }
  std::vector<double>* avgDepth = new std::vector<double>(this->chrVec->size(), 0.0);
  std::string currChr;
  std::vector<double>* currDepthVec = nullptr;
  for(unsigned int chrIdx = 0; chrIdx < this->chrVec->size(); chrIdx++) {
    currChr = (*this->chrVec)[chrIdx];
    if(getDiploid) {
      currDepthVec = (*this->chrToDiploidDepthMap)[currChr];
    } else {
      currDepthVec = (*this->chrToTumorDepthMap)[currChr];
    }
    (*avgDepth)[chrIdx] = std::accumulate(currDepthVec->begin(), currDepthVec->end(), 0.0);
    (*avgDepth)[chrIdx] = (*avgDepth)[chrIdx] / currDepthVec->size();
  }
  return avgDepth;
}

/*
 * helper function to get total diploid depth
 */
double DepthPair::getTotalDiploidDepth() {
  double sum = 0;
  std::string currChr;
  std::vector<double>* currDepthVec = nullptr;
  for(unsigned int chrIdx = 0; chrIdx < this->chrVec->size(); chrIdx++) {
    currChr = (*this->chrVec)[chrIdx];
    currDepthVec = (*this->chrToDiploidDepthMap)[currChr];
    sum += std::accumulate(currDepthVec->begin(), currDepthVec->end(), 0.0);
  }
  return sum;
}
/*
 * helper function to get total tumor depth
 */
double DepthPair::getTotalTumorDepth() {
  double sum = 0;
  std::string currChr;
  std::vector<double>* currDepthVec = nullptr;
  for(unsigned int chrIdx = 0; chrIdx < this->chrVec->size(); chrIdx++) {
    currChr = (*this->chrVec)[chrIdx];
    currDepthVec = (*this->chrToTumorDepthMap)[currChr];
    sum += std::accumulate(currDepthVec->begin(), currDepthVec->end(), 0.0);
  }
  return sum;
}

/*
 * function to save simulated depths like a depth file so that simulations can be read
 * back in in subsequent tests
 * format is tab separated:
 *   chr | start | end | depth [| diploid variance]
 */
void DepthPair::saveSimDiploidAsDepthFile(std::string filename) {
  std::ofstream outFile(filename);
  std::string sep = "\t";
  std::vector<std::string>* currRegions = nullptr;
  std::vector<double>* currDiploidDepth = nullptr;
  std::vector<double>* currDiploidVar = nullptr;
  for(std::vector<std::string>::iterator chrItr = this->chrVec->begin(); chrItr != this->chrVec->end(); ++chrItr) {
    currRegions = (*this->regions)[*chrItr];
    currDiploidDepth = (*this->chrToDiploidDepthMap)[*chrItr];
    currDiploidVar = (*this->chrToDiploidVarMap)[*chrItr];

    for(unsigned int i = 0; i < currRegions->size(); i++) {
      // deparse key stored earlier (chr:start-end) into tab separated values
      std::string s = (*currRegions)[i];
      std::replace(s.begin(), s.end(), ':', *(sep.c_str()));
      std::replace(s.begin(), s.end(), '-', *(sep.c_str()));
      outFile << s << sep << (*currDiploidDepth)[i] << sep << (*currDiploidVar)[i] << std::endl;
    }
  }
  outFile.close();
}
void DepthPair::saveSimTumorAsDepthFile(std::string filename) {
  std::ofstream outFile(filename);
  std::string sep = "\t";
  std::vector<std::string>* currRegions = nullptr;
  std::vector<double>* currTumorDepth = nullptr;
  for(std::vector<std::string>::iterator chrItr = this->chrVec->begin(); chrItr != this->chrVec->end(); ++chrItr) {
    currRegions = (*this->regions)[*chrItr];
    currTumorDepth = (*this->chrToTumorDepthMap)[*chrItr];

    for(unsigned int i = 0; i < currRegions->size(); i++) {
      // deparse key stored earlier (chr:start-end) into tab separated values
      std::string s = (*currRegions)[i];
      std::replace(s.begin(), s.end(), ':', *(sep.c_str()));
      std::replace(s.begin(), s.end(), '-', *(sep.c_str()));
      outFile << s << sep << (*currTumorDepth)[i] << std::endl;
    }
  }
  outFile.close();
}
