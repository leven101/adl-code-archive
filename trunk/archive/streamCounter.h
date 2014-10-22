#ifndef STREAM_COUNTER_H
#define STREAM_COUNTER_H

#include <iostream>
#include <sstream>
#include <map>
#include <pthread.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <unicode/brkiter.h>
#include "types.h"
#include "file.h"

class StreamCounter {
  public:
    StreamCounter();
    ~StreamCounter();
    void countStreamItems(std::istringstream& stream);
    void extractRawData(std::istringstream& stream);
    void printRaw();
    void print(string prefix = "");
    void printSentences();
  private:
    size_t total_tokens_;
    date_t period_;
    const short uni_idx_, bi_idx_, tri_idx_;
    date_cnt_t* ngrm_cntrs_[NGRAM_ORDER];
    void printCounts(const string& what, const string& filename, 
                        date_cnt_t& cntr);
    void clearCntrs();
    void sort(const string& filename);
    word_t longstring_;
};

#endif //STREAM_COUNTER_H
