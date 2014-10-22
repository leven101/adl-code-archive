#ifndef online_model1_h
#define online_model1_h
#include <functional>
#include "file.h"
#include "utils.h"
#include "vocab.h"
#include "types.h"
#include "google/sparse_hash_map"

using google::sparse_hash_map;
#define MAX_W 457979
/* I must become the model */

/*----------- Defnition of Hash Function (from GIZA) ------------------*/
typedef std::pair<wordID_t, wordID_t> wrdPair_t;
class hashpair : public std::unary_function<std::pair<wordID_t, wordID_t>, count_t >
{
public:
  count_t operator() (const std::pair<wordID_t, wordID_t>& key) const {
    return (count_t) MAX_W*key.first + key.second; //MAX_W is prime 
  }
};

class Model1 {
public:
  Model1(string source, string target): srcPath(source), trgPath(target) {
    srcVcb = new Vocab(false);
    trgVcb = new Vocab(false);
    wrdProbs_ = new ttable_t;
    counts_ = new ttable_t;
    loadCorpora(srcPath, trgPath);
  }
  virtual ~Model1() {
    delete srcVcb;
    delete trgVcb;
    if(counts_) delete counts_;
    if(wrdProbs_) delete wrdProbs_;
  }
  virtual int EM(const int noIters);
  void printTTable(string fname = "");
protected:
  typedef std::pair<wrdPair_t, float> pairStat_t;
  typedef sparse_hash_map<wrdPair_t, float, hashpair, std::equal_to<wrdPair_t> > ttable_t;
  Vocab *srcVcb, *trgVcb;
  int totSrcToks, totTrgToks;
  ttable_t *wrdProbs_;
  ttable_t *counts_;
  string srcPath, trgPath;
  std::vector<std::vector<wordID_t> > src_sents, trg_sents; // holds all aligned sentences in memory
private:
  void loadCorpora(string srcCorpus, string trgCorpus);
  void initializeUniformly();
  void initZeroCounts();
};
#endif
