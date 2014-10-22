#include <iostream>
#include <cstdlib>
#include <algorithm>
#include "utils.h"
#include "params.h"
#include "dirList.h"
#include "types.h"
#include "vocab.h"
#include "file.h"

#include <Trie.cc>

struct nodeValue {
  size_t count;
  bool hit;
};

#ifdef INSTANTIATE_TEMPLATES
INSTANTIATE_TRIE(wordID_t, nodeValue);
#endif

using namespace std;

const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"top-dir", "/disk/scratch/s0674000/data/reuters/tokenized/", "dir", Parameters::kStringValue, "Directory with training data."},
  {"data-ext", "week*.grams.gz", "ext", Parameters::kStringValue, "Data file's extensions."},
  {"order", "3", "ord", Parameters::kIntValue, "ngram order."},
  {"mfiOrder", "1", "mfi", Parameters::kIntValue, "mfi bound."},
  {"baseline", "0", "bl", Parameters::kBoolValue, "if this is a baseline run"},
  {"time", "1", "tf", Parameters::kIntValue, "the time window to stream over"},
  {"memory", "15", "m", Parameters::kIntValue, "Upper bound size in MBs."}
};

// variables and types
typedef Trie<wordID_t, nodeValue> node2_t;
static node2_t trie2_;
typedef Trie<wordID_t, size_t> node_t;
static node_t MFIs_;
unsigned maxBytes_, 
         totalEntries_ = 0,
         cur_period_ = 0, 
         last_period_ = 0,
         time_frame_ = 2,
         corpusSize_ = 0,
         order_,
         mfiOrder_;
Parameters* params;
Vocab* vocab;
vector<double> blScores;

// function definitions
int processData();
void insert(const size_t* IDs, const size_t count);
void insertLM(const size_t* IDs, const size_t count);
void batchUpdateLM();
void dynamicLMStuff(); 
void memStats(); 
void dumpTrie(node2_t&);
void clearHits(node2_t&);
size_t countChildren(node2_t& head);
size_t countChildren(node_t& head);
size_t heuristicDelete(size_t numToDel, size_t depth);
void trackMFIs(const size_t* IDs, const size_t count);
size_t corpusSize();

void printGramSet(node2_t& head, FileHandler* fout) {
  static int idx = 0;
  static vector<int> gram(order_);
  TrieIter<wordID_t, nodeValue> trieIter_(head);
  wordID_t key;
  node2_t* node2;
  while(node2 = trieIter_.next(key)) {
    string sout = "";
    gram[idx] = key;
    for(int i = 0; i <= idx; i++) 
       sout += vocab->getWord(gram[i]) + " ";
    Utils::trim(sout);
    sout += "\t" + Utils::IntToStr(node2->value().count);
    *fout << sout << endl;
    idx++;
    printGramSet(*node2, fout);
    idx--;
  }
}

float trieMem() { 
  return totalEntries_ * (sizeof(node_t)) * 1.7;
}  

size_t getPeriod(std::string str) {
  const date_t period = "week";
  return atoi(str.substr(str.rfind(period) 
    + period.size()).c_str());
}

int processData() {
  unsigned no_key;
  Map_noKey(no_key);
  DirList dir_(params->getParamValue("top-dir"),    // which directory
               params->getParamValue("data-ext"));  // which file extensions
  iterate(dir_.cFileList_, itr) {                   // open each file
    //extract time period
    cur_period_ = getPeriod(*itr);
    cerr << "Processing week " << cur_period_ << endl;
    FileHandler gramsFile(*itr, ios::in);
    string line; 
    while(getline(gramsFile, line)) {
      // split up ngram and associated count
      vector<string> linedata;
      assert(Utils::splitToStr(line, linedata, "\t") == 2);
      line = linedata[0];
      size_t count = atoi(linedata[1].c_str());
      // split up ngram
      Utils::splitToStr(line, linedata, " ");
      size_t wrdIDs[linedata.size() + 1];
      for(size_t i=0; i < linedata.size(); i++) {
        // map to vocab
        wrdIDs[i] = vocab->getWordID(linedata[i]);
      }
      wrdIDs[linedata.size()] = no_key;
      insert(wrdIDs,count);
    }
  }
  batchUpdateLM();
  printGramSet(trie2_, new FileHandler("52.subset.350MB.gz", ios::out));
  return 0;
}

size_t corpusSize() {
  // training corpus size is sum of unigram counts
  size_t count = 0;
  node2_t* node2;
  size_t key;
  TrieIter<wordID_t, nodeValue> itr(trie2_);
  while(node2 = (itr.next(key))) 
    count += node2->value().count;
  return count;
}

/*NOTE: Could also update LM when a certain number of new items have been found */

void insert(const size_t* IDs, const size_t count) {
// Fill trie untill space param is reached. (Time period 0)
  static bool full = false;
  static unsigned first_period = 0;
  if(!full) {
    if(trieMem() < maxBytes_) {
      insertLM(IDs, count);
    }
    else {
      cerr << "Memory bound reached. (period=" << cur_period_ << ")\n";
      memStats();
      last_period_ = cur_period_ - time_frame_;
      first_period = cur_period_;
      full = true;
      insertLM(IDs, count);
    }
  } // finish inserting current data to avoid smoothing problems
  else if(cur_period_ == first_period) {
    insertLM(IDs, count);
  }
  else {
    dynamicLMStuff();
    trackMFIs(IDs,count); //track MFI stream items
  }
}

void insertLM(const size_t* IDs, const size_t count) {
  nodeValue* nodeVal;
  bool fnd = true;
  nodeVal = trie2_.insert(IDs, fnd);
  nodeVal->count += count;
  nodeVal->hit = false;
  if(!fnd) totalEntries_++;
}

void dynamicLMStuff() {
  // if memory bound reached update after designated time period
  if(time_frame_ <= abs(int(cur_period_ - last_period_))) {
    batchUpdateLM();
    //corpusSize_ = corpusSize();
    if(cur_period_ == 20 || cur_period_ == 36) {
      string fname(Utils::IntToStr(cur_period_ - 1) + ".subset");
      fname += "." + params->getParamValue("memory") + "MB.gz";
      FileHandler fout(fname, ios::out, false);
      cout << "printing " << fname << endl;
      printGramSet(trie2_, &fout);
      cout << "done printing " << fname << endl;
      clearHits(trie2_);
    }
    last_period_ = cur_period_;
  }
}

void batchUpdateLM() {
/* This funcion
 *   - iterates through the current MFIs
 *   - adds each to LM 
 *   - deletes nodes from the main LM
 */
  //cerr << "start batch update...\n";
  if(countChildren(MFIs_) < 100) return;
  size_t inserted = 0, deleted = 0,
    order = order_; 
  node_t* newItem;
  node2_t* node2;
  size_t* count;
  bool fnd;
  // insert each item in the MFIs data structure
  for(size_t i=1; i <= order; i++) {
    deleted = 0;
    wordID_t keys[i];  
    TrieIter2<wordID_t, size_t> iter2(MFIs_, keys, i);
    while(newItem = iter2.next()) {
      node2 = trie2_.insertTrie(keys, fnd);
      if(!fnd) deleted++;
      node2->value().hit = true;
      if(node2->value().count < newItem->value())
        node2->value().count = newItem->value();
    }
    deleted = heuristicDelete(deleted, i);
  }
  // clear old MFIs
  MFIs_.clear();
  //cerr << "finish batch update...\n";
}

void testDelete() {
  FileHandler fout("___stdout___", std::ios::out);
  printGramSet(trie2_, &fout);
  dumpTrie(trie2_);
  trie2_.remove(3);
  printGramSet(trie2_, &fout);
  dumpTrie(trie2_);
  node2_t* node2 = trie2_.findTrie(9)->findTrie(10)->findTrie(11);
  cerr << "NumEntries for ``pricings'' = "<< node2->numEntries() << endl;
}

size_t heuristicDelete(size_t num2Del, size_t depth) {
  /*
   * This deletes only those nodes which
   * haven't been "hit" in test time.
   */
  size_t deleted = 0, order = order_; 
  node2_t* node2;
  //for(int i = order; i > 0; i--) {
  for(int i = depth; i <= order; i++) {
    wordID_t keys[i];
    TrieIter2<wordID_t, nodeValue> itr2(trie2_, keys, i);
    while((node2 = itr2.next()) && (deleted < num2Del)) {
      if((!node2->value().hit) && (node2->numEntries() < 1)) {
        trie2_.remove(keys);
        deleted++;
      }
    }
    if(deleted >= num2Del) break;
  }
  return deleted;
}

void clearHits(node2_t& head) {
  TrieIter<wordID_t, nodeValue> trieIter_(head);
  wordID_t key;
  node2_t* node2;
  while(node2 = trieIter_.next(key)) {
    node2->value().hit = false; 
    clearHits(*node2);
  }
}

void dumpTrie(node2_t& head) {
  TrieIter<wordID_t, nodeValue> trieIter_(head);
  wordID_t key;
  node2_t* node2;
  static int indent = 0;
  indent += 2;
  while(node2 = trieIter_.next(key)) {
    for(int i = 0; i < indent; i++) cerr << " ";
    cout << key << " / " << node2->value().count  << endl;
    dumpTrie(*node2);
  }
  indent -= 2;
}

void memStats() {
  MemStats memuse;
  trie2_.memStats(memuse);
  memuse.print();
}

size_t countChildren(node2_t& head) {
  TrieIter<wordID_t, nodeValue> trieIter_(head);
  wordID_t key;
  node2_t *node2;
  size_t children = head.numEntries();
  while(node2 = trieIter_.next(key)) {
    children += countChildren(*node2);
  }
  return children;
}
size_t countChildren(node_t& head) {
  TrieIter<wordID_t, size_t> trieIter_(head);
  wordID_t key;
  node_t *node;
  size_t children = head.numEntries();
  while(node = trieIter_.next(key)) {
    children += countChildren(*node);
  }
  return children;
}

void trackMFIs(const size_t* IDs, const size_t count) {
  // add to current MFIs DS if not in prior MFIs
  // with a minimum count.
  if(count >= mfiOrder_)
     *MFIs_.insert(IDs) += count;
}

/*int main(int argc, char** argv) {
  params = new Parameters(argc, argv, paramdefs, NumOfParams(paramdefs));
  maxBytes_ = atoi(params->getParamValue("memory").c_str()) * (1UL << 20);
  order_ = atoi(params->getParamValue("order").c_str());
  mfiOrder_ = atoi(params->getParamValue("mfiOrder").c_str());
  time_frame_ = atoi(params->getParamValue("time").c_str());
  vocab = new Vocab();
  processData();
  delete vocab;
  delete params;
  return 0;
}*/
int main(int argc, char**argv) {
  unsigned no_key;
  Map_noKey(no_key);
  vocab = new Vocab();
  FileHandler gramsFile("/disk/scratch/s0674000/data/reuters/aclRandLMExperiments/270milTrainingGrams.gz", ios::in);
  string line;
  int r(0), oneoffs(0);
  vector<string> linedata;
  while(getline(gramsFile, line) && (++r < 7000000)) {
    // split up ngram and associated count
    assert(Utils::splitToStr(line, linedata, "\t") == 2);
    line = linedata[0];
    size_t count = atoi(linedata[1].c_str());
    if(count == 1) ++oneoffs;
    // split up ngram
    Utils::splitToStr(line, linedata, " ");
    size_t wrdIDs[linedata.size() + 1];
    for(size_t i=0; i < linedata.size(); i++) {
      // map to vocab
      wrdIDs[i] = vocab->getWordID(linedata[i]);
    }
    wrdIDs[linedata.size()] = no_key;
    insertLM(wrdIDs,count);
  }
  gramsFile.close();
  memStats();
  return 1;
  //dumpTrie(trie2_);
  node2_t* node2;
  cerr << "trie2_.numEntries() = " << trie2_.numEntries() << endl;
  cerr << "Total number of entries with frequency 1 = " << oneoffs 
    << " of " << r << " events " << endl;
  for(int i = 0; i < 5; i++) {
    wordID_t keys[i];
    int lvlcnt(0), gtone(0), gtthree(0), gtten(0), max(0);
    TrieIter2<wordID_t, nodeValue> itr2(trie2_, keys, i);
    while(node2 = itr2.next()) { 
      lvlcnt += node2->numEntries();
      if(node2->numEntries() > max) max = node2->numEntries();
      if(node2->numEntries() > 1) ++gtone;
      if(node2->numEntries() > 3) ++gtthree;
      if(node2->numEntries() > 10) ++gtten;
    }
    cerr << "Level " << i << " :  " << lvlcnt << " nodes";
    cerr << " (" << gtone << " > 1, " << gtthree << " > 3, " << gtten << " > 10)\n";
    cerr << "max = " << max << endl;
  }
  delete vocab;
  return 1;
}


/***********************************
 **** CODE USED FOR MSE TESTS *****
***********************************
double getScore(size_t* keys, size_t order) {
  float alpha = 1.0;
  size_t gramCnt = 0, subgramCnt = 0, i; 
  node2_t* node2 = 0;
  for(i = 0; i < order; i++) {
    int idx = i;  // start from history keys[i]
    gramCnt = subgramCnt = 0;
    node2 = &trie2_;
    // descend to the longest available trie path
    while(idx < order) {
      node2 = node2->findTrie(keys[idx]);
      if(node2) {
        node2->value().hit = true;  // record hit
        subgramCnt = gramCnt;       // remember history count
        gramCnt = node2->value().count; // get current count
        idx++;
      }
      else {
        break;
      }
    }
    if(idx == order) 
      break;  // stop if found full ngram
    else 
      alpha *= 0.4;
  } 
  // if nothing was found
  if(gramCnt == 0) {
    return 0;
  } // handle unigram case
  else if(subgramCnt == 0 && gramCnt > 0 
      && i == order - 1) {
    subgramCnt = corpusSize_;
    alpha = 1;
  }
  double result = alpha * (float(gramCnt) / float(subgramCnt));
  if(gramCnt > subgramCnt) {
    cerr << "Result=" << result;
    cerr << ", i=" << i << endl;
    for(int j=0; j < order; j++)
      cerr << vocab->getWord(keys[j]) << " "; 
    cerr << endl << "GramCnt=" << gramCnt;
    cerr << ", SubGramCnt=" << subgramCnt << endl << endl;
    return 0;
  }
  return result; 
}

void testLM() {
  unsigned no_key;
  Map_noKey(no_key);
  vector<double> scores;
  DirList testDir(params->getParamValue("top-dir") + "test/",
                  params->getParamValue("data-ext"));
  iterate(testDir.cFileList_, itr) {
    // for each testFile
    FileHandler testFile(*itr, ios::in);
    string line;
    while(getline(testFile, line)) {
      // split up ngram and associated count
      vector<string> linedata;
      Utils::splitToStr(line, linedata, "\t");
      line = linedata[0];
      // split up ngram
      Utils::splitToStr(line, linedata, " ");
      size_t wrdIDs[linedata.size() + 1];
      for(int i=0; i < linedata.size(); i++) {
        // map to vocab
        wrdIDs[i] = vocab->getWordID(linedata[i]);
      }
      wrdIDs[linedata.size()] = no_key;
      //cerr << getScore(wrdIDs, linedata.size()) << endl;
      scores.push_back(getScore(wrdIDs, linedata.size())); 
    }
  }
  
  if(params->getParamValue("baseline") == "1") {
    // save scores to file
    FileHandler scoreFile("baseline.scores", ios::out, false);
    iterate(scores, itr)
      scoreFile << *itr << endl;
  }
  else { 
    // if scores haven't been read in yet
    if(blScores.empty()) {
      //cerr << "Reading in baseline scores...\n";
      // read in baseline scores from file
      FileHandler scoreFile("baseline.scores", ios::in);
      string line;
      while(getline(scoreFile, line)){
        blScores.push_back(atof(line.c_str()));
      }
    }
    // output MSE
    double mse = 0, diff = 0;
    for(int i=0; i < scores.size(); i++) {
      diff = blScores[i] - scores[i];
      mse += (diff * diff);
    }
    //print stuff
    cerr << "Current period=" << cur_period_ << "\nMSE=" << mse/float(scores.size()) << endl;
  }
}
 */ 
