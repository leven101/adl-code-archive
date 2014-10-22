#include "multiOnlineRLM.h"

int main(int argc, char** argv) {
  MultiOnlineRLM<unsigned> mlm("/home/abby/workspace/experiments/current/orlms/", 5);
 
  FileHandler ftest("./test", std::ios::in);
  string line;
  std::vector<string> ngram, vtmp;
  while(getline(ftest, line)) {
    Utils::trim(line);
    if(Utils::splitToStr(line, ngram, " ") > 0) { // split tokens in ngram
      int len = ngram.size();
      cout <<  mlm.getProb(ngram, len, NULL) << endl;
    }
  }
  ftest.close();
  return 0;
}
