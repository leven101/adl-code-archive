#ifndef inc_orlm_stats_h
#define inc_orlm_stats_h
#include "types.h"


class Stats {
public:
  Stats(int streamNum, int order) {
    streamNum_ = streamNum;
    order_ = order;
    order_counts_ = new count_t*[streamNum];
    for(int i=0; i < streamNum; ++i) {
      order_counts_[i] = new count_t[order];
      for(int j=0; j < order; ++j) 
        order_counts_[i][j] = 0;
    }
  }
  ~Stats() {
    for(int i=0; i < streamNum_; ++i) {
      delete order_counts_[i];
    }
    delete order_counts_;
  }
  count_t get(int stream, int order) {
    assert((stream <= streamNum_) && (order <= order_));
    assert((stream > 0) && (order > 0));
    return order_counts_[stream-1][order-1];
  }
  void add(int stream, int order) {
    assert((stream <= streamNum_) && (order <= order_));
    assert((stream > 0) && (order > 0));
    ++order_counts_[stream-1][order-1];
  }
  void print() {
    for(int i=0; i < order_; ++i)
      std::cout << "\t" << (i+1);
    std::cout << endl;
    for(int i=0; i < streamNum_; ++i) {
      std::cout << (i+1);
      for(int j=0; j < order_; ++j)
        std::cout << "\t" << order_counts_[i][j];
      std::cout << endl;
    }
  }
private:
  count_t** order_counts_;
  int streamNum_, order_;
};
#endif
