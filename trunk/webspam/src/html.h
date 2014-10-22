#include <map>
#include <set>
#include "htmlParser.h"
#include "dirList.h"
#include "onlineRLM.h"
#include "vocab.h"

static int order_ = 3;
map<int, bool> train_labels_;
vector<vector<float> > features_;
Vocab vcb_(false);
static OnlineRLM<unsigned> *spLM=0, *nspLM=0, *seedLM=0;
static const string punc = " \t.\",+/*:-()'><;#?][%&";
map<string, int> countStreamNgrams(string&, const int, bool = false);
static void extractLMFeatures(string&, const int, vector<float>&);
static void loadTrainData() {
  vector<string> tmpV;
  string line;
  string infile = "../data/2007/training/webspam-uk2007-set1-labels.txt"; 
  FileHandler f2(infile, ios::in);
  while(getline(f2, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int id = atoi(tmpV[0].c_str());
    Utils::trim(tmpV[1]);
    train_labels_[id] = tmpV[1] == "spam" ? 1 : 0;
  }
  f2.close();
}
static void processWarcFile(int id, vector<float>& feats) {
  string infile = "../data/2007/raw/sites/" + Utils::IntToStr(id) + ".warc.gz";
  if(!Utils::fileExists(infile)) {
    FileHandler f(infile, ios::in);
    bool add(false);
    string shtml(""), line, text;
    while(getline(f, line)) {
      Utils::trim(line);
      if(line.empty()) continue;
      if(COMPARE("warc/0.9", line.substr(0,8).c_str())) { // for each webpage
        if(!shtml.empty()) {
          text += parseHtml(shtml);
          shtml.clear();
        }
        add = false;
      }
      else if(COMPARE("<html", line.substr(0,5).c_str())) {
        add = true;
      }
      if(add) {
        shtml +=  " " + line;
      }
    } // end while
    extractLMFeatures(text, id, feats);
  }
}
static void extractLMFeatures(string& text, const int nodeID, vector<float>& feats) {
  // for each website get entropy, kl-divergence from base, kl-divergence from spamLM, nsLM
  if(!(seedLM && spLM && nspLM)) {
    spLM = new OnlineRLM<unsigned>(50, 20, 50, 3, &vcb_, 8); // LM for nonspam   
    nspLM = new OnlineRLM<unsigned>(300, 20, 50, 3, &vcb_, 8); // LM for nonspam   
    FileHandler fseed("../data/2007/lms/seedLM.orlm.gz", ios::in|ios::binary);
    seedLM = new OnlineRLM<unsigned>(&fseed, order_);
    fseed.close();
  }
  assert(seedLM != 0);
  assert(spLM != 0);
  assert(nspLM != 0);
  if(train_labels_.size() == 0)
    loadTrainData();
  std::vector<string> ngram;
  static bool bSpam(false), bNonspam(false);
  float tot1grams(0), tot2grams(0), tot3grams(0); // for counts
  float unigramH(0), bigramH(0), trigramH(0); // for entropy
  float uniKLSeedToSp(0), uniKLSeedToNsp(0), uniKLNspToSp(0);
  float biKLSeedToSp(0), biKLSeedToNsp(0), biKLNspToSp(0);
  float triKLSeedToSp(0), triKLSeedToNsp(0), triKLNspToSp(0);
  map<int, bool>::iterator lblIt = train_labels_.find(nodeID);
  if(lblIt != train_labels_.end()) {
    if(lblIt->second) { // spam
      bSpam = true;
      bNonspam = false;
    }
    else  { // nonspam
      bNonspam = true;
      bSpam = false;
    }
  }
  map<string, int> allgrams;
  string line;
  for(size_t i=0; i < text.size(); ++i) {
    if(text[i] != '\n')
      line += text[i];
    else {
      map<string, int> grams = countStreamNgrams(line, order_, true);
      iterate(grams, git)
        allgrams[git->first] += git->second;
      line.clear();
    }
  }
  iterate(allgrams, git) { // for each ngram 
    Utils::splitToStr(git->first, ngram, " ");
    //cerr << "'" << git->first << "'" << endl;
    // add to stream LMs
    if(bSpam) {
      spLM->update(ngram, git->second);
    }
    else if(bNonspam) {
      nspLM->update(ngram, git->second);
    }
    int len = ngram.size();
    wordID_t wrdIDs[len];
    for(int i=0; i < len; i++) {
      Utils::trim(ngram[i]);
      wrdIDs[i] = seedLM->vocab_->getWordID(ngram[i]);
    }
    double seedLgPrb = seedLM->getProb(wrdIDs, len, NULL);
    double seedPrb = pow(10, seedLgPrb); // LM return log10 prob
    for(int i=0; i < len; i++) {
      Utils::trim(ngram[i]);
      wrdIDs[i] = vcb_.getWordID(ngram[i]);
    }
    double spLgPrb = spLM->getProb(wrdIDs, len, NULL);
    double spPrb = pow(10, spLgPrb); // LM return log10 prob
    double nspLgPrb = nspLM->getProb(wrdIDs, len, NULL);
    double nspPrb = pow(10, nspLgPrb); // LM return log10 prob
    switch(len) {
      case 1:
        tot1grams += git->second;
        unigramH += seedPrb * seedLgPrb;
        uniKLSeedToSp += seedPrb * log10(seedPrb / spPrb);
        uniKLSeedToNsp += seedPrb * log10(seedPrb /nspPrb);
        uniKLNspToSp += nspPrb * log10(nspPrb / spPrb);
        break;
      case 2:
        tot2grams += git->second;
        bigramH += seedPrb * seedLgPrb;
        biKLSeedToSp += seedPrb * log10(seedPrb / spPrb);
        biKLSeedToNsp += seedPrb * log10(seedPrb /nspPrb);
        biKLNspToSp += nspPrb * log10(nspPrb / spPrb);
        break;
      case 3:
        tot3grams += git->second;
        trigramH += seedPrb * seedLgPrb;
        triKLSeedToSp += seedPrb * log10(seedPrb / spPrb);
        triKLSeedToNsp += seedPrb * log10(seedPrb /nspPrb);
        triKLNspToSp += nspPrb * log10(nspPrb / spPrb);
        break;
      default: // do nada
        break;
    }
  }
  feats.push_back(tot1grams);
  feats.push_back(tot2grams);
  feats.push_back(tot3grams);
  feats.push_back(unigramH * -1);
  feats.push_back(bigramH * -1);
  feats.push_back(trigramH* -1);
  feats.push_back(uniKLSeedToSp);
  feats.push_back(uniKLSeedToNsp);
  feats.push_back(uniKLNspToSp);
  feats.push_back(biKLSeedToSp);
  feats.push_back(biKLSeedToNsp);
  feats.push_back(biKLNspToSp);
  feats.push_back(triKLSeedToSp);
  feats.push_back(triKLSeedToNsp);
  feats.push_back(triKLNspToSp);
}
map<string, int> countStreamNgrams(string& text, const int order, bool add2Vcb) {
  map<string, int> grams;
  istringstream stream(text);
  string tokens[order], token;
  for(int i=0; i < order; ++i)
    tokens[i] = "";
  int tokenCnt(0), idx(0); 
  while(tokenCnt < order-1) { // get first ngram 
    if(stream >> token) {
      Utils::trim(token, punc); //remove punctuation
      if(token.empty()) continue;
      Utils::strToLowercase(token);
      tokens[tokenCnt] = token; 
      ++tokenCnt;
      if(add2Vcb) vcb_.getWordID(token);
    }
    else {
      break;
    }
  }
  if(stream.eof()) { // handle case where stream size < order
    string ngram("");
    for(int n=0; n < tokenCnt; ++n) {
      ngram += " " + tokens[n];
      Utils::trim(ngram);
      ++grams[ngram];
      //cerr << n << ": \"" << ngram << "\"" << endl;
    }
    return grams;
  }
  // while stream exists build all ngrams from order 1 to order-1.
  while(stream >> token) {
    Utils::trim(token, punc); //remove punctuation
    if(token.empty()) continue;
    Utils::strToLowercase(token);
    tokens[tokenCnt % order] = token;
    if(add2Vcb) vcb_.getWordID(token);
    string ngram("");
    //build and record all order ngrams
    for(int n = order-1; n >= 0; --n) {
      idx = (tokenCnt - n) % order;
      ngram += " " + tokens[idx];
      Utils::trim(ngram);
      ++grams[ngram];
      //cerr << n << ": \"" << ngram << "\"" << endl;
    }
    ++tokenCnt; 
  }
  int lastidx = (idx + 1) % order;
  tokens[lastidx] = "";
  // Here we have all order n grams. get remainder of ngrams in cache
  for(int i=order-1; i >=0; --i) {
    string ngram("");
    for(int n = i; n >= 1; --n) {
      idx = (tokenCnt - n) % order;
      ngram += " " + tokens[idx];
      Utils::trim(ngram);
      ++grams[ngram];
      //cerr << n << ": \"" << ngram << "\"" << endl;
    }
  }
  return grams;
}
