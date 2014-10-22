#ifndef online_hmm_h
#define online_hmm_h
#include <algorithm>
#include "model1.h"
using std::vector;
// emission probs are translation probs (seeded from model1)
// observations/states (source, target sentences).

class HMM: public Model1 {
public:
  HMM(string source, string target): Model1(source, target), 
  alpha_(0), beta_(0), maxSntLen_(100) {
    alignPrb_ = new vector<ttable_t>(maxSntLen_ + 1);
    alignCnts_ = new vector<ttable_t>(maxSntLen_ + 1);
    wordID_t s = trgVcb->getWordID("s"), t = trgVcb->getWordID("t");
    (*wrdProbs_)[wrdPair_t(srcVcb->getWordID("A"), s)] = .4;
    (*wrdProbs_)[wrdPair_t(srcVcb->getWordID("B"), s)] = .6;
    (*wrdProbs_)[wrdPair_t(srcVcb->getWordID("A"), t)] = .5;
    (*wrdProbs_)[wrdPair_t(srcVcb->getWordID("B"), t)] = .5;
    (*alignPrb_)[0][wrdPair_t(0,s)] = .85; 
    (*alignPrb_)[0][wrdPair_t(0,t)] = .15; 
    (*alignPrb_)[0][wrdPair_t(s,s)] = .3; 
    (*alignPrb_)[0][wrdPair_t(s,t)] = .7; 
    (*alignPrb_)[0][wrdPair_t(t,t)] = .9; 
    (*alignPrb_)[0][wrdPair_t(t,s)] = .1;
  }
  ~HMM() {
    if(alignPrb_) delete alignPrb_;
    if(alignCnts_) delete alignCnts_;
  }
  int EM(int iterations); 
private: 
  float **alpha_, **beta_;
  const int maxSntLen_;
  float pr0_;
  vector<ttable_t>* alignPrb_; //jump probabilites (transition probs, alignment probs)
  vector<ttable_t>* alignCnts_;  // to collect counts for jump probs 
  std::map<wordID_t, float> totTrgPrb_;
  float forward(vector<wordID_t>&, vector<wordID_t>&);
  float backward(vector<wordID_t>&, vector<wordID_t>&);
  float count(vector<wordID_t>&, vector<wordID_t>&, const float);
  void newTrellis(const int rows, const int cols);
  void freeTrellis(const int rows);
  void clearCounts();
  float normalize();
};
#endif
