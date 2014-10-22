#ifndef INC_MULTI_DYNAMICLM_H
#define INC_MULTI_DYNAMICLM_H
#include "onlineRLM.h"
#include "dirList.h"

template<typename T>
class MultiOnlineRLM {
  // constructor builds all ORLMs that are in directory
public: 
  MultiOnlineRLM(string sdir, int order, string weightsFile) {
    DirList dir(sdir, "*.orlm.gz");
    assert(dir.cFileList_.size() > 0);
    numOfLMs_ = dir.cFileList_.size();
    std::cerr << "numOfLms_ = " << numOfLMs_ << endl;
    readWeights(weightsFile);
    lms_ = new OnlineRLM<T>*[numOfLMs_];
    for(int i=0; i < numOfLMs_; ++i) lms_[i] = 0;
    for(int i=0; i < numOfLMs_; ++i) {
      if(weights_[i] == 0) continue;
      std::cerr << "Loading ORLM at " << dir.cFileList_[i] << std::endl;
      FileHandler fin(dir.cFileList_[i], std::ios::in|std::ios::binary, true);
      lms_[i] = new OnlineRLM<T>(&fin, order);
      fin.close();
    }
  }
  ~MultiOnlineRLM() {
    for(int i=0; i < numOfLMs_; ++i) {
      if(lms_[i]) delete lms_[i];
    }
    delete[] lms_;
    delete[] weights_;
  }
  float getProb(const std::vector<string>& sngram, int len, const void** state) {
    //iterate(sngram, itr) std::cerr << *itr << "|";
    // gets each lms score
    float scores[numOfLMs_];
    for(int i=0; i < numOfLMs_; ++i) scores[i] = 0;
    for(int i=0; i < numOfLMs_; ++i) {
      if(weights_[i] == 0) continue;
      wordID_t ngramIDs[len]; 
      for(int j=0; j < len; ++j) { // translate to vocab
        ngramIDs[j] = lms_[i]->vocab_->getWordID(sngram[j]);
      }
      scores[i] = lms_[i]->getProb(&ngramIDs[0], len, state);
    }
    float score(0);
    for(int i=0; i < numOfLMs_; ++i) 
      score += (weights_[i] * scores[i]);
    return score;
  }
  void clearCache() {
    for(int i=0; i < numOfLMs_; ++i) {
      if(weights_[i] == 0) continue;
      else lms_[i]->clearCache();
    }
  }
private:
  OnlineRLM<T>** lms_;
  float* weights_;
  int numOfLMs_;
  void readWeights(string weightsFilePath) {
    //const string weightsFilePath = "/home/abby/workspace/experiments/current/weights.txt";
    cerr << "Reading weights from " << weightsFilePath << endl;
    weights_ = new float[numOfLMs_];
    FileHandler fweights(weightsFilePath, std::ios::in, true);
    string line;
    for(int i=0; i < numOfLMs_; ++i) {
      if(!getline(fweights, line)) {
        std::cerr << "Under weights\n" << std::endl;
        exit(1);
      }
      else {
        weights_[i] = atof(line.c_str());
        std::cerr << "LM[" << i << "] has weight " << weights_[i] << endl;
      }
    }
    fweights.close();
  }
};
#endif
