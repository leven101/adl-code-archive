#ifndef bloomier_h_
#define bloomier_h_

#include <set>
#include <algorithm>
#include "hash.h"
#include "vocab.h"
#include "quantizer.h"
#include "RandLMFilter.h"
#include "RandLMCache.h"
#include "dirList.h"

using randlm::Filter;
using randlm::Cache;

static bool strict_checks_ = true;
class BloomierXX {
public:
  BloomierXX(string dataDir, int order) {
    cerr << "loading bloomier filter lm...\n";
    // load general info
    FileHandler finfo(dataDir + "/filter.info.gz", std::ios::in);
    vocab_ = new Vocab(&finfo);
    fgh_ = new UnivHash_linear<count_t>(&finfo);
    finfo.read((char*)&cells_, sizeof(cells_));
    finfo.read((char*)&cellWidth_, sizeof(cellWidth_));
    cerr << "fingerprint bits = " << cellWidth_ << endl;
    finfo.read((char*)&maxCode_, sizeof(maxCode_));
    cerr << "maxcode = " << maxCode_ << endl;
    finfo.read((char*)&corpusSize_, sizeof(corpusSize_));
    cerr << "corpus size = " << corpusSize_ << endl;
    finfo.read((char*)&tot_subsets_, sizeof(tot_subsets_));
    finfo.read((char*)&order_, sizeof(order_));
    qtizer_ = new LogQtizer(&finfo);
    lastGramIDs = new wordID_t*[tot_subsets_];
    for(int i = 0; i < tot_subsets_; ++i) {
      lastGramIDs[i] = new wordID_t[order_];
      for(int j = 0; j < order_; ++j)
        finfo.read((char*)&lastGramIDs[i][j], sizeof(lastGramIDs[i][j]));
    }
    finfo.close();
    // load all subset filters and associated hash functions
    DirList dir(dataDir, "*.bfilter.gz");
    assert((int)dir.cFileList_.size() == tot_subsets_);
    FileHandler* fdata(0);
    filter2_ = new Filter<count_t>*[tot_subsets_];
    h2_ = new UnivHash_linear<count_t>*[tot_subsets_];
    for(int subset = 0; subset < tot_subsets_; ++subset) {
      string fname = dir.cFileList_[subset];
      fdata = new FileHandler(fname, std::ios::in|std::ios::binary);
      filter2_[subset] = new Filter<count_t>(fdata);
      h2_[subset] = new UnivHash_linear<count_t>(fdata);
      fdata->close();
      delete fdata;
    }
    cerr << "filter mem usage = " << filter2_[0]->size() * tot_subsets_ << "MB\n";
    // lm stuff
    alpha_ = new float[order + 1];
    for(int i = 0; i <= order_; ++i) 
      alpha_[i] = i * log10(0.4);
    oovprob_ = log10(1.0 / (static_cast<float>(vocab_->size()) - 1));
    cache_ = new Cache<float>(8888.8888, 9999.9999); // unknown_value, null_value
    // get string representation of final ngram of each subset 
    lastGramStrs = new string[tot_subsets_];
    for(int i = 0; i < tot_subsets_; ++i) {
      lastGramStrs[i] = "";
      for(int j = 0; j < order_; ++j)
        lastGramStrs[i] += lastGramIDs[i][j] == 0 ? "" 
          : vocab_->getWord(lastGramIDs[i][j]) + " ";
      Utils::trim(lastGramStrs[i]);
      //cerr << "lastGramStr[" << i << "] = " << lastGramStrs[i] << endl;
    }        
  }
  ~BloomierXX() {
    delete vocab_;
    delete fgh_;
    delete cache_;
    for(int i = 0; i < tot_subsets_; ++i) {
      delete[] lastGramIDs[i];
      delete filter2_[i];
      delete h2_[i];
    }
    if(lastGramIDs) delete[] lastGramIDs;
    if(h2_) delete h2_;
    if(filter2_) delete filter2_;
    if(lastGramStrs) delete[] lastGramStrs;
  }
  static count_t nonZeroSignature(wordID_t* IDs, int len, 
    UnivHash_linear<count_t>* fh) {
    count_t fingerprint(0);
    int h(0);
    do {
      fingerprint = fh->hash(IDs, len, h);
    } while((fingerprint == 0) && (++h < 10));
    if(fingerprint == 0) 
      cerr << "WARNING: Unable to find non-zero signature for ngram\n" << endl;
    return fingerprint;
  }
  int findSubset(wordID_t* ngram, int len) {
    string sngram("");
    for(int i = 0; i < len; ++i)
      sngram += vocab_->getWord(ngram[i]) + " ";
    Utils::trim(sngram);
    for(int i = 0; i < tot_subsets_; ++i)
      if(sngram.compare(lastGramStrs[i]) <= 0)
        return i;
    return -1;
  }
  int query(wordID_t* ngram, int len) {
    // find which filter to query
    int subset = findSubset(ngram, len);
    if(subset < 0) return 0;
    // hash ngram for indexes and fingerprint
    count_t fgp = BloomierXX::nonZeroSignature(ngram, len, fgh_);
    std::set<count_t> indexes;
    for(int k = 0; k < 3; ++k) // for all k hash functions
      indexes.insert(h2_[subset]->hash(ngram, len, k));  // get each unique index associated with this ngram
    count_t g = cells_ + 1;
    iterate(indexes, itr) {
      count_t arrVal = filter2_[subset]->read(*itr);   // get value at index
      g = (g == (count_t)cells_ + 1 ? arrVal : (g ^ arrVal));  // xor everything together
    }
    g ^= fgp;
    // retreive code
    if(g <= (count_t)qtizer_->maxcode()) return (int)qtizer_->value(g);  // if code <= maxCode then return value
    //if(g <= (count_t)maxCode_) return (int)qtizer_->value(g);  // if code <= maxCode then return value
    else return 0;
  }
  void clearCache() {if(cache_) cache_->clear();}
  Vocab* vocab_;
  Filter<count_t>** filter2_;
  UnivHash_linear<count_t>** h2_;
  UnivHash_linear<count_t>* fgh_;
  int maxCode_;
  int cells_;
  int cellWidth_;
  int corpusSize_;
  int order_;
  int tot_subsets_;
  wordID_t** lastGramIDs;
  string* lastGramStrs;
  float* alpha_;
  float oovprob_;
  Cache<float>* cache_;
  LogQtizer* qtizer_;
  float getProb(wordID_t* ngram, int len, const void** state) {
    float logprob(0);
    const void* context = (state) ? *state : 0;
    // if full ngram and prob not in cache
    if(!cache_->checkCacheNgram(ngram, len, &logprob, &context)) {
      // get full prob and put in cache
      int num_fnd(0), den_val(0);
      int in[len]; // in[] keeps counts of increasing order numerator 
      for(int i = len - 1; i >= 0; --i) {
        if(ngram[i] == Vocab::kOOVWordID) break;  // no need to query if OOV
        in[i] = query(&ngram[i], len - i);
        if(in[i] > 0) {
          num_fnd = len - i;
        }
        else if(strict_checks_) break;
      }
      while(num_fnd > 1) {
        //get sub-context of size one less than length found (exluding target) 
        if(((den_val = query(&ngram[len - num_fnd], num_fnd - 1)) > 0) &&
            den_val >= in[len - num_fnd]) {
          break;
        }
        else --num_fnd; // else backoff to lower ngram order 
      }
      // find prob
      switch(num_fnd) {
        case 0: // OOV
          logprob = alpha_[len] + oovprob_;
          break;
        case 1: // unigram found only
          assert(in[len - 1] > 0);
          logprob = alpha_[len - 1] + 
            log10(static_cast<float>(in[len - 1]) / static_cast<float>(corpusSize_));
          break;
        default:
          assert(den_val > 0);
          //if(subgram == in[len - found]) ++subgram; // avoid returning zero probs????
          logprob = alpha_[len - num_fnd] + 
            log10(static_cast<float>(in[len - num_fnd]) / static_cast<float>(den_val));
          break;
      }
      // need unique context
      context = getContext(&ngram[len - num_fnd], num_fnd);
      // put whatever was found in cache
      cache_->setCacheNgram(ngram, len, logprob, context);
    } // end checkCache
    return logprob; 
  }
  const void* getContext(wordID_t* ngram, int len) {
    int dummy(0);
    float* addresses[len];  // only interested in addresses of cache
    assert(cache_->getCache2(ngram, len, &addresses[0], &dummy) == len);
    // return address of cache node
    return (const void*)addresses[0]; 
  }
};
#endif
