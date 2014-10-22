#ifndef hiero_h
#define hiero_h

#include <boost/functional/hash.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <pthread.h>
#include "file.h"
#include "params.h"
#include "vocab.h"
#include "types.h"
#include "pyp.hh"
#include "tree.hh"

using std::vector;
using std::map;
using std::pair;
using std::set;

#define MAX_SENT_LEN 1000;

struct hieroRule;
// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"source", "/Users/ablev/work/data/fbis-mini/chinese", "s", Parameters::kStringValue, ""},
  {"target", "/Users/ablev/work/data/fbis-mini/english", "t", Parameters::kStringValue, ""},
  {"model1", "/Users/ablev/work/data/fbis-mini/zh-en.t1.15", "m1", Parameters::kStringValue, ""},
  {"model1-inv", "/Users/ablev/work/data/fbis-mini/en-zh.t1.15", "m1i", Parameters::kStringValue, ""},
  {"outputdir", "", "o", Parameters::kStringValue, ""},
  {"init-aligns", "/Users/ablev/work/data/fbis-mini/model1.gdfa", "m1a", Parameters::kStringValue, ""},
  {"anneal", Parameters::kFalseValue, "", Parameters::kBoolValue, ""},
  {"no-init", Parameters::kFalseValue, "", Parameters::kBoolValue, ""},
  {"phr-init", Parameters::kFalseValue, "", Parameters::kBoolValue, ""},
  {"restart", Parameters::kFalseValue, "", Parameters::kBoolValue, ""},
  {"discount", "0.8", "a", Parameters::kFloatValue, ""},
  {"strength", "1", "b", Parameters::kFloatValue, ""},
  {"poisson-mean", "1", "p", Parameters::kFloatValue, ""},
  {"terminal-penalty", ".9", "tp", Parameters::kFloatValue, ""},
  {"iterations", "10", "i", Parameters::kIntValue, ""},
  {"maxSnts", "0", "tot", Parameters::kIntValue, ""},
  {"threads", "1", "", Parameters::kIntValue, ""},
  {"snt2Debug", "0", "s2d", Parameters::kIntValue, ""},
};
std::ostream& operator<<(std::ostream& out, const hieroRule& X);
inline bool isTerm(int i) {
  return i > 0;
}
struct hieroRule {
  vector<int> s,t;
  size_t operator()(const hieroRule& hr) const {
    size_t h(0);
    boost::hash_combine(h, boost::hash_value(hr.s));
    boost::hash_combine(h, boost::hash_value(hr.t));
    return h;
  }
  bool operator()(const hieroRule& hr1, const hieroRule& hr2) const {
    // only equal if all symbols in rule match in order
    if(hr1.s == hr2.s && (hr1.t == hr2.t))
      return true;
    return false;
  }
  bool operator==(const hieroRule& hr) const {
    if(this->s == hr.s && (this->t == hr.t))
      return true;
    return false;
  }
  int arity() const {
    int a(0);
    iterate(s, sit) {
      if(!isTerm(*sit))
        ++a;
    }
    return a;
  }
};

// global variables
Vocab *srcVcb_, *trgVcb_;
Parameters *params_;
vector<vector<wordID_t> > src_sents, trg_sents;
unordered_map<pair<wordID_t, wordID_t>, pair<float, float>, boost::hash<pair<wordID_t, wordID_t> > > m1_;
PYP<hieroRule, hieroRule>* restr_;
unordered_map<hieroRule, int, hieroRule> timesSampled;
int snt2debug_, noSamplesChecked_, noSamplesBuilt_, currIter_, totIter_;
bool initPhrs_, noInit_;
double tp_;
std::stringstream outputDir_;
struct node_t {
  node_t(): NTIndex(0) {}
  hieroRule rule;
  int label;
  int NTIndex;
  int sstr, send;
  // int tstr, tend; currently not used
  bool operator==(const int& lbl) const {
    return (this->label == lbl);
  }
};
vector<tree<node_t>* > trees_;
pthread_mutex_t mutexsum;
pthread_attr_t attr;
struct threadData_t {
  int sntBegin;
  int sntEnd;
  float anneal;
  int maxspan;
  bool skipDecr; 
};
//functions
void freeStuff() {
  delete srcVcb_;
  delete trgVcb_;
  delete params_;
  delete restr_;
}
wordID_t getWordID(string word, Vocab* vcb) {
  if(word == "NULL")
    return Vocab::kOOVWordID;
  else
    return vcb->getWordID(word);
}
double log_factorial(size_t n) { 
  // return log of factorial
  return lgamma(n + 1);
  //return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
double geomDist(float k, float p) {
  //k = num getting prob for, p = prob of success
  return pow(1.0f-p, k - 1.0f) * (p);
}
double log_geomDist(float k, float p) {
  return log(1.0f-p) * (k-1) + log(p); 
}
double poissDist(size_t k, float lambda) {
  // k = number getting probability for, lambda = mean of distribution
  double num = pow(lambda, k) * exp(-lambda);
  double denom = exp(log_factorial(k)); 
  return num / denom;
}
double log_poissDist(size_t k, float lambda) {
  // k = number getting probability for, lambda = mean of distribution
  double num = (log(lambda) * k)  - lambda;
  double denom = log_factorial(k); 
  return num - denom;
}
double get_p0(const hieroRule* hr, const bool);
void loadData();
void loadM1Params();
void loadSentences();
void initPYP();
vector<int> getMostProbAlgs(const int);
void fullSentenceInit();
tree<node_t>* buildXRuleGrid(const vector<int>&, const int);
bool sampleChildOn(const hieroRule&, const int, const int, const double, 
  const double, tree<node_t>::iterator, tree<node_t>*, const float, const bool, const int,
  const int);
bool sampleChildOff(const int, const int, const int, const double, 
  const double, tree<node_t>::iterator, tree<node_t>*, const float, const bool);
void* sample(void*);
void printTreeNode(tree<node_t>::iterator node);
void printAlignments(int);
void getTermIdxs(const tree<node_t>::iterator, map<int, set<int> >&, int*, 
  const tree<node_t>&, bool); 
void runSampling();
void readTreeState();
void saveTreeState();
double dRand() {
  return drand48();
  //return rand() * (1/double(RAND_MAX));
}
int iRand(int modval) {
  return lrand48() % modval;
  //return rand() % modval;
}
float termPenalty() {
  if(((float)currIter_/(float)totIter_) < .9)
    return tp_;
  else 
    return .55;
}
int maxSpan() {
  int span = MAX_SENT_LEN;
  if(currIter_ <= 3) {
    span = 1;
  }
  else if(currIter_ <= 6) {
    span= 2;
  }
  else if(currIter_ <= 10) {
    span = 3;
  }
  return span;
}
float annealT() {
  static const bool withAnneal = params_->getParam("anneal") 
    == Parameters::kTrueValue ? true: false;
  // need an annealing schedule dependento on currIter_ation
  if(withAnneal) {
    if(currIter_ <= 10) {
      return .2;
    }
    else if(currIter_ <= 50) {
      return .5;
    }
    else if(currIter_ <= 100) {
      return .75;
    }
    else if(currIter_ <= 200) {
      return 1;
    }
  }
  return 1;
}
bool alreadyOn(const hieroRule& rule) {
  bool on(false);
  if((rule.s == rule.t) && (rule.s.size() == 1) && (!isTerm(rule.s[0]))) // same nonterm
    on = true;
  return on;
}
bool validRule(const hieroRule& rule, const hieroRule& parent, bool* newSrcSide,
  bool* newTrgSide, set<int>& srcNTs, set<int>& trgNTs, int* oldLeftmostTrgSym) {
  if(rule.s == parent.s && (rule.t == parent.t)) { 
    return false;
  }
  // all nonterminals must match
  bool valid(false);
  if(*newSrcSide) {
    srcNTs.clear();
    trgNTs.clear();
    iterate(rule.s, sit) {
      if(!isTerm(*sit))
        srcNTs.insert(*sit);
    }
    iterate(rule.t, tit) {
      if(!isTerm(*tit))
        trgNTs.insert(*tit);
    }
    *newSrcSide = false;
    assert(*newTrgSide);
    *oldLeftmostTrgSym = rule.t[0];
    *newTrgSide = false;
  }
  else if(*newTrgSide) { // same source, new target
    trgNTs.clear();
    iterate(rule.t, tit) {
      if(!isTerm(*tit))
        trgNTs.insert(*tit);
    }
    *oldLeftmostTrgSym = rule.t[0];
    *newTrgSide = false;
  }
  else { 
    // rule has 'slid' one symbol to the left 
    if(!isTerm(*oldLeftmostTrgSym))
      trgNTs.erase(*oldLeftmostTrgSym);
    if(!isTerm(rule.t.back()))
      trgNTs.insert(rule.t.back());
    *oldLeftmostTrgSym = rule.t[0];
  }
  // all NTs match 
  valid = (srcNTs == trgNTs); 
  return valid;
}
hieroRule combineChildParent(const hieroRule& prtSpan, const hieroRule& ocr, 
  const int sIdx, const int tIdx, const bool debug) {
  // have child span so merge with parent span 
  hieroRule npr = prtSpan;
  npr.s.erase(npr.s.begin() + sIdx);
  npr.s.insert(npr.s.begin() + sIdx, ocr.s.begin(),
    ocr.s.end());
  npr.t.erase(npr.t.begin() + tIdx);
  npr.t.insert(npr.t.begin() + tIdx, ocr.t.begin(),
    ocr.t.end());
  return npr;
}
hieroRule newParentRule1(const hieroRule& oldRule, const hieroRule& newChild,
  const int sIdx, const int tIdx, const int lbl) {
  hieroRule npr = oldRule;
  npr.s.erase(npr.s.begin() + sIdx, npr.s.begin() + sIdx + newChild.s.size());
  npr.s.insert(npr.s.begin() + sIdx, lbl);
  npr.t.erase(npr.t.begin() + tIdx, npr.t.begin() + tIdx + newChild.t.size());
  npr.t.insert(npr.t.begin() + tIdx, lbl);
  return npr;
}
void updateTree(tree<node_t>::iterator node, tree<node_t>* fullTree, 
  const hieroRule& npr, const hieroRule& ncr, const int lbl, const int ss,
  const int se) {
  node->rule = npr;
  node_t newNode; // create new node
  newNode.rule = ncr;
  newNode.label = lbl;
  newNode.sstr = ss;
  newNode.send = se;
  tree<node_t>::iterator newNodeItr, ntNodeItr;
  newNodeItr = fullTree->append_child(node, newNode); // add new node to tree
  // move all nodes with labels in the new child span so new child becomes their parent 
  iterate(ncr.s, sit) { // assume good rule here since checked previously
    if(!isTerm(*sit)) {
      ntNodeItr = std::find(node.begin(), node.end(), *sit);
      assert(ntNodeItr != node.end());
      fullTree->append_child(newNodeItr, ntNodeItr);
      fullTree->erase(ntNodeItr);
    }
  }
}
hieroRule transform(const hieroRule& hr) {
  // transform rule so NTs match in order but reordering is consistent 
  // between source and target sides for all individual rules
  hieroRule copy = hr;
  int NTidx(-1);
  for(size_t i=0; i < hr.s.size(); ++i) {
    if(!isTerm(hr.s[i])) {
      for(size_t j=0; j < hr.t.size(); ++j) {
        if(hr.t[j] == hr.s[i]) {
          copy.s[i] = NTidx;
          copy.t[j] = NTidx;
          --NTidx;
          break;
        }
      }
    }
  }
  return copy;
}
void printTree(int snt) {
  cout << endl;
  tree<node_t>::iterator trit2 = trees_[snt]->begin();
  while(trit2 != trees_[snt]->end()) {
    printTreeNode(trit2);
    ++trit2;
  }
  cout << endl;
}
int increment(const hieroRule& hr, const double p0) {
  pthread_mutex_lock (&mutexsum);
  int val = restr_->increment(hr, p0); 
  pthread_mutex_unlock (&mutexsum);
  return val;
}
int decrement(const hieroRule& hr) {
  pthread_mutex_lock (&mutexsum);
  int val = restr_->decrement(hr); 
  pthread_mutex_unlock (&mutexsum);
  return val;
}
double dataLLAndOtherStuff(const int iter) {
  std::stringstream fname;
  fname << outputDir_.str() << "/rules." << iter;
  FileHandler fout(fname.str(), std::ios::out, false); 
  const double restrProb = restr_->log_restaurant_prob();
  double allPriorProbs(0);
  float avgSrcSideLen(0), avgTrgSideLen(0);
  map<pair<int,int>, int> ruleLens;
  PYP<hieroRule, hieroRule>::iterator rit = restr_->begin();
  for(; rit != restr_->end(); ++rit) { // for each rule type
    double numTables = restr_->num_tables(rit->first); // get the number of tables it's at
    // calculate its prior times the number of tables
    double prb = numTables * get_p0(&(rit->first), false); 
    if(prb != -std::numeric_limits<double>::infinity())
      allPriorProbs += prb;
    int srclen = rit->first.s.size();
    int trglen = rit->first.t.size();
    avgSrcSideLen += srclen; 
    avgTrgSideLen += trglen;
    ruleLens[std::make_pair(srclen, trglen)] += rit->second;
    fout << rit->first << "\t" << rit->second << "\t" 
      << restr_->log_prob(rit->first, get_p0(&(rit->first), false)) << endl;
  }
  fout.close();
  cout << "Average src/trg lengths: " << (avgSrcSideLen / restr_->num_types()) 
    << " / " << (avgTrgSideLen / restr_->num_types()) << endl; 
  /*int currSrcLen(1);
  iterate(ruleLens, lit) {
    if(lit->first.first != currSrcLen) {
      currSrcLen = lit->first.first;
      cout << endl;
    }
    cout << "<" << lit->first.first << "," << lit->first.second << "> = "
      << lit->second << "  ";
  }
  cout << endl;*/
  return restrProb + allPriorProbs;
}
vector<int> precomputeSpans(tree<node_t>::iterator trit, bool debug) {
  const hieroRule& hr = trit->rule;
  vector<int> symSpans(hr.s.size() + 1, 0);
  tree<node_t>::iterator it; 
  /*if(debug) {
    cout << "PrecomputeSpans" << endl;
    printTreeNode(trit);
  }*/
  for(size_t j=0; j < symSpans.size(); ++j) { // precompute span of each index in rule 
    if(j==0) {
      symSpans[0] = trit->sstr;
    }
    else {
      int i = j-1;
      if(!isTerm(hr.s[i])) {
      it = std::find(trit.begin(), trit.end(), hr.s[i]); 
      assert(it != trit.end());
      symSpans[j] = it->send;
      }
      else {
        symSpans[j] = 1 + symSpans[i];
      }
    }
    //if(debug) cout << j << "--" << symSpans[j] << endl;
  }
  return symSpans;
}
map<int, pair<double, double> > cacheChildPrbs(tree<node_t>::iterator node,
  double* allChildPrbs, bool debug) {
  map<int, pair<double, double> > mchildPrbs;
  *allChildPrbs = 0;
  // get probability of dependent part of sentence derivation
  iterate(node->rule.s, s) {
    if(!isTerm(*s)) {
      tree<node_t>::iterator child = std::find(node.begin(), node.end(), *s);
      hieroRule crTr = transform(child->rule);
      decrement(crTr);
      double p0 = get_p0(&child->rule, debug);
      double p = restr_->log_prob(crTr, p0);
      mchildPrbs[*s] = std::make_pair(p, p0);
      *allChildPrbs += p;
      increment(crTr, exp(p0));
    }
  }
  return mchildPrbs;
}
struct resample_term_penalty {
  typedef map<int, pair<int, int> > len_t;
  len_t lengths_; //<total terminals, total NTs> in restr for rules of this length
  resample_term_penalty() {
    histogram();
  }
  void histogram() {
    // build histogram of src lengths
    PYP<hieroRule, hieroRule>::iterator rit = restr_->begin();
    for(; rit != restr_->end(); ++rit) { // for each rule type
      const hieroRule& r = rit->first;
      int arity = r.arity();
      int noTables = restr_->num_tables(r); 
      len_t::iterator lit = lengths_.find(r.s.size());
      if(lit == lengths_.end()) {
        lengths_[r.s.size()] = pair<int, int>(0, 0);
        lit = lengths_.find(r.s.size());
      }
      lit->second.first += ((r.s.size() - arity) * noTables); 
      lit->second.second += (arity * noTables); 
    }
  }
  double operator() (double prop_tp) const {
    // for each src length get the probability with the proposed tp 
    double log_prob=0;
    iterate(lengths_, lit) {
      double tp = pow(prop_tp, lit->first);
      double noTrms = lit->second.first;
      double noNTs = lit->second.second;
      log_prob += (noTrms * log(tp)) + (noNTs * log(1-tp)); 
    }
    return log_prob;
  }
};
void resampleTermPenalty() {
  int niterations = 5;
  cerr << "Before resampling tp: " << tp_ << endl;
  resample_term_penalty rtp;
  tp_ = slice_sampler1d(rtp, tp_, mt_genrand_res53, std::numeric_limits<double>::min(), 
    (double) 1.0, (double) 0.0, niterations, 100*niterations);
  cerr << "After resampling tp: " << tp_ << endl;
}
void printArityHist() {
  map<int, float> hist;
  float tot(0);
  PYP<hieroRule, hieroRule>::iterator rit = restr_->begin();
  for(; rit != restr_->end(); ++rit) { // for each rule type
    const hieroRule& r = rit->first;
    int count = rit->second;
    int arity = r.arity();
    hist[arity] += count;
    tot += count;
  }
  iterate(hist, hit) 
    cout << hit->first << "\t" << (hit->second/tot) << "\t" << hit->second << endl;
}
#endif
