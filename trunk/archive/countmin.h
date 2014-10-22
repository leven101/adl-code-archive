#ifndef INC_COUNTMINSKETCH_S0674000_H
#define INC_COUNTMINSKETCH_S0674000_H

#include "hash.h"
#include "types.h"

class CMSketch {
public:
  CMSketch(int depth, int width): depth_(depth), width_(width) {
    assert(depth > 0 && width > 0);
    //get memory for counts with depth and width
    counts_ = new int*[depth_];
    for(int i = 0; i < depth_; ++i) { 
      counts_[i] = new int[width_];
      for(int j = 0; j < width_; ++j) 
        counts_[i][j] = 0; // initialize to 0 
    }
    //initialize UHF with number of function and range
    h_ = new UnivHash_linear<unsigned>(width_, depth_, PRIME);
  }
  ~CMSketch() {
    // free counts
    for(int i = 0; i < depth_; ++i) {  
      delete[] counts_[i];
    }
    delete[] counts_;
    // free hash functions
    delete h_;
  }
  void update(wordID_t* item, int len, int diff) {
    for(int i = 0; i < depth_; ++i) 
      counts_[i][h_->hash(item, len, i)] += diff;
  }
  int pointEst(wordID_t* query, int len) {
    // return min value in counts
    int est = counts_[0][h_->hash(query, len, 0)];
    for(int i = 1; i < depth_; ++i) 
      est = std::min(est, counts_[i][h_->hash(query, len, i)]);
    return est;
  }
private:
  const int depth_;
  const int width_;
  int** counts_;
  UnivHash_linear<unsigned> *h_;
};



#endif
