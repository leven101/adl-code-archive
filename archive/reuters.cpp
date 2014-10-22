#include <iostream>
#include <iomanip>
#include <ctime>
#include <pthread.h>
#include "reuters.h"
#include "streamCounter.h"

const std::string corpus = "/group/corpora/restricted/reuters/english/",
                  local = "../data/";
// function pointer definition
typedef void (StreamCounter::*FN_PTR)(std::istringstream&);

int main( void ) {
  std::string reutersDir = corpus, tmpDir = "../data/tmp/";
  StreamCounter strmCnt;
  ReutersCorpusStream reutStream(reutersDir, tmpDir, false);
  FN_PTR fptr = &StreamCounter::extractRawData; 
  time_t begin = time(NULL);
  reutStream.stream<StreamCounter, FN_PTR> (strmCnt, fptr);
  strmCnt.printRaw();
  time_t end = time(NULL);
  //std::cout << "Time taken = " << std::setprecision(5) << difftime(end,begin) / 60 << " minutes\n";
  return 0;
}
