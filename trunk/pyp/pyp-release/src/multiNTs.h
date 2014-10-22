#ifndef hiero_h
#define hiero_h

#include <boost/functional/hash.hpp>
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
using std::make_pair;


// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, commandline abbreviation, data type, help message
  {"source", "/Users/ablev/work/data/fbis-mini/chinese", "s", Parameters::kStringValue, "source langauge file (one sentence per line)"},
  {"target", "/Users/ablev/work/data/fbis-mini/english", "t", Parameters::kStringValue, "target language file (one sentence per line)"},
  {"model1", "/Users/ablev/work/data/fbis-mini/zh-en.t1.15", "m1", Parameters::kStringValue, "IBM Model1 source to target word translation probabilities"},
  {"model1-inv", "/Users/ablev/work/data/fbis-mini/en-zh.t1.15", "m1i", Parameters::kStringValue, "IBM Model1 target to source word translation probabilities"},
  {"init-aligns", "/Users/ablev/work/data/fbis-mini/model1.gdfa", "m1a", Parameters::kStringValue, "word alignments for initialization"},
  {"outputdir", "", "o", Parameters::kStringValue, "directory to save output files to"},
  {"anneal", Parameters::kFalseValue, "", Parameters::kBoolValue, "turn on annealing"},
  {"no-init", Parameters::kFalseValue, "", Parameters::kBoolValue, "use full sentence initialization"},
  {"phr-init", Parameters::kFalseValue, "", Parameters::kBoolValue, "use multiword initialization"},
  {"restart", Parameters::kFalseValue, "", Parameters::kBoolValue, "restart sampling from prior sampling run"},
  {"block-sample", Parameters::kFalseValue, "bs", Parameters::kBoolValue, "run the block sampler"},
  {"bs-every", "8", "", Parameters::kIntValue, "how often to run the block sampler between binary sample iterations"},
  {"discount", "0.8", "a", Parameters::kFloatValue, "discount parameter for the PYP"},
  {"strength", "1", "b", Parameters::kFloatValue, "strength parameter for the PYP"},
  {"poisson-mean", "1", "p", Parameters::kFloatValue, "poisson mean parameter (used in the base distribution)"},
  {"terminal-penalty", ".55", "tp", Parameters::kFloatValue, "terminal penalty parameter (used in the base distribution)"},
  {"iterations", "10", "i", Parameters::kIntValue, "total number of sampling iterations"},
  {"maxSnts", "0", "tot", Parameters::kIntValue, "maximum number of sentences to load from the data files"},
  {"threads", "1", "", Parameters::kIntValue, "number of threads to use"},
  {"snt2Debug", "0", "s2d", Parameters::kIntValue, "sentence number to debug (set to -1 to turn off debugging)"},
  {"hierarchical", Parameters::kFalseValue, "hier", Parameters::kBoolValue, "use multiple nonterminals for the SCFG"},
  {"restrict-scfg", Parameters::kTrueValue, "X", Parameters::kBoolValue, "use a single nonterminal X for the SCFG"},
  {"src-nts", "/Users/ablev/work/data/fbis-mini/mkcls/zh.out", "swc", Parameters::kStringValue, "source nonterminal class file"},
  {"trg-nts", "/Users/ablev/work/data/fbis-mini/mkcls/en.out", "twc", Parameters::kStringValue, "target nonterminal class file"},
};
//#define MAX_SENT_LEN 1000;
inline int MAX_SENT_LEN() {return 1000;}
inline bool isTerm(int i) {return i>0;}

struct hieroRule {
  int NT;
  vector<int> s,t;
  size_t operator()(const hieroRule& hr) const {
    size_t h(0);
    boost::hash_combine(h, boost::hash_value(NT));
    boost::hash_combine(h, boost::hash_value(hr.s));
    boost::hash_combine(h, boost::hash_value(hr.t));
    return h;
  }
  bool operator()(const hieroRule& hr1, const hieroRule& hr2) const {
    // only equal if all symbols in rule match in order
    if((hr1.s == hr2.s && (hr1.t == hr2.t)) && (hr1.NT == hr2.NT))
      return true;
    return false;
  }
  bool operator==(const hieroRule& hr) const {
    if((this->s == hr.s && (this->t == hr.t)) && (this->NT == hr.NT))
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
struct node_t {
  hieroRule rule;
  int label; // node's label in tree
  int nodeIndex; // usually empty. highest label for full tree
  int sstr, send;
  // int tstr, tend; currently not used
  bool operator==(const int& lbl) const {
    return (label == lbl);
  }
};
class MultiNT;
class GrammarCache;
struct threadData_t {
  int sntBegin;
  int sntEnd;
  float anneal;
  int maxspan;
  bool skipDecr; 
  MultiNT* mnt;
  //GrammarCache* cache;
};
std::ostream& operator<<(std::ostream& out, const hieroRule& X);
// global variables
Vocab *srcVcb_, *trgVcb_, *nts_;
pthread_mutex_t mutexsum;
pthread_mutex_t mutexsum_NT; // mutex for NT creation
pthread_attr_t attr;
class MultiNT {
public: 
  void loadLongSnts();
  MultiNT() {};
  ~MultiNT() { 
    delete params_;
    delete srcVcb_;
    delete trgVcb_;
    delete restr_;
    if(hierch_) delete base_restr_;
    iterate(trees_, trit) delete *trit;
  }
  void loadData();
  void initPYP();
  void runSampling();
  void initialize(int, char**);
// rest private? 
  void loadM1Params();
  void loadSentences();
  void loadNonterminals();
  void fullSentenceInit();
  void printAlignments(int);
  int readTreeState();
  void saveTreeState();
  //derivation and probs
  typedef pair<double, double> p0_t;
  p0_t get_p0(const hieroRule*, const bool);
  double base_p0(const hieroRule*, const bool, double&);
  tree<node_t>* buildXRuleGrid(const vector<int>&, const int);
  //variables
  std::stringstream outputDir_;
  PYP<hieroRule, hieroRule>* base_restr_;
  vector<PYP<hieroRule, hieroRule> >* restr_;
  vector<tree<node_t>* > trees_;
  vector<vector<wordID_t> > src_sents, trg_sents;
  unordered_map<pair<wordID_t, wordID_t>, pair<float, float>, boost::hash<pair<wordID_t, wordID_t> > > m1_;
  map<wordID_t, int> srcWrdCls_, trgWrdCls_;
  int snt2debug_, noSamplesChecked_, noSamplesBuilt_, currIter_, numSrcWrdCls_, totIter_;
  bool initPhrs_, noInit_, hierch_, singleNT_, blockSample_;
  double a_, b_, tp_; // hyper-parameters
  void getTermIdxs(const tree<node_t>::iterator, map<int, set<int> >&, int*, 
    const tree<node_t>&, bool); 
  Parameters *params_;
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
  void printSntPair(const int snt) {
    for(size_t i=0; i < src_sents[snt].size(); ++i) {
      cerr << srcVcb_->getWord(src_sents[snt][i]) << " ";
    }
    cerr << "||| ";
    for(size_t i=0; i < trg_sents[snt].size(); ++i) {
      cerr << trgVcb_->getWord(trg_sents[snt][i]) << " ";
    }
    cerr << endl;
  }
  float termPenalty() {
    return tp_;
  }
  int maxSpan() {
    int span = MAX_SENT_LEN();
    if(currIter_ <= 3) {
      span = 3;
    }
    else if(currIter_ <= 6) {
      span = 4;
    }
    else if(currIter_ <= 10) {
      span= 5;
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
  void printTreeNode(tree<node_t>::iterator node) const {
    cout << node->label << ": " << node->rule;
    cout << "[" << node->sstr << "-" << node->send << " | ?-? ]";
    if(node.node->parent)
      cout << "  (parent = " << node.node->parent->data.label << ")";
    cout << endl;
  }
  void printTree(int snt) const {
    cout << endl;
    tree<node_t>::iterator trit2 = trees_[snt]->begin();
    while(trit2 != trees_[snt]->end()) {
      printTreeNode(trit2);
      ++trit2;
    }
    cout << endl;
  }
  void printTree(const tree<node_t>* trp) const {
    cout << endl;
    tree<node_t>::iterator trit2 = trp->begin();
    while(trit2 != trp->end()) {
      printTreeNode(trit2);
      ++trit2;
    }
    cout << endl;
  }
  hieroRule transform(const hieroRule& hr, const tree<node_t>::iterator trit,
     const int missingNT=-1) const {
    if(singleNT_) 
      return transform(hr, false);
    // LHS -> negative number is NT symbol
    // RHS -> negative number is index relation on LHS 
    hieroRule copy = hr;
    tree<node_t>::iterator it; 
    for(int i=0; i < (int)hr.s.size(); ++i) { 
      if(!isTerm(hr.s[i])) {  // each tree index on src side get's replaced by it's NT symbol
        //it = std::find(trit.begin(), trit.end(), hr.s[i]);  // doesn't work!?!?
        it = trit.begin(); 
        while(it != trit.end()) {
          if(it->label == hr.s[i]) break;
          ++it;
        }
        if(it != trit.end()) {
          copy.s[i] = (it->rule.NT) * -1; // store -NT in the PYP
        }
        else {
          assert(missingNT != -1);
          copy.s[i] = missingNT * -1;
        }
        // each NT on trg side get's index of src side
        for(int j=0; j < (int)hr.t.size(); ++j) {
          if(hr.t[j] == hr.s[i]) {
            copy.t[j] = (-1 * i); // is this right?  
          }
        }
      }
    }
    return copy;
  }
  hieroRule transform(const hieroRule& hr, bool restr_format) const {
    // transform hiero-style rule so Xs match in order but reordering is consistent 
    // between source and target sides for all individual rules
    hieroRule copy = hr;
    int NTidx(-1);
    for(size_t i=0; i < hr.s.size(); ++i) {
      if(!isTerm(hr.s[i])) {
        for(size_t j=0; j < hr.t.size(); ++j) {
          if(restr_format && (hr.t[j] == ((int)i * -1))) {
            copy.s[i] = NTidx;
            copy.t[j] = NTidx;
            --NTidx;
            break;
          }
          else if(hr.t[j] == hr.s[i]) {
            copy.s[i] = NTidx;
            copy.t[j] = NTidx;
            --NTidx;
            break;
          }
        }
      }
    }
    copy.NT = 0; //reset to single NT
    return copy;
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
  int NTSymbol(const wordID_t w1, const wordID_t w2) {
    if(singleNT_) // if just a single NT symbol 
      return 0;
    map<wordID_t, int>::iterator itr = srcWrdCls_.find(w1);
    assert(itr != srcWrdCls_.end());
    int id1 = itr->second;
    itr = srcWrdCls_.find(w2);
    assert(itr != srcWrdCls_.end());
    int id2 = itr->second;
    string nt = Utils::IntToStr(id1);
    nt += "-" + Utils::IntToStr(id2);
    pthread_mutex_lock(&mutexsum_NT);
    int ret = nts_->getWordID(nt);
    while(ret > (int)restr_->size()) { // ensure PYP for every NT
      restr_->push_back(PYP<hieroRule, hieroRule>(a_,b_));
    }
    pthread_mutex_unlock(&mutexsum_NT);
    return (ret-1); // keep in sync with zero-index vector restr_
  }
  int increment(const hieroRule& hr, const p0_t p0) {
    int nt = hr.NT;
    assert(nt < (int)restr_->size());
    pthread_mutex_lock(&mutexsum);
    int val = restr_->at(nt).increment(hr, exp(p0.second));
    if(val && hierch_) { // new table
      base_restr_->increment(transform(hr, true), exp(p0.first));
      // assert(base_restr_->num_customers() == total_tables_);
    }
    pthread_mutex_unlock(&mutexsum);
    return val;
  }
  int decrement(const hieroRule& hr) {
    int nt = hr.NT;
    assert((int)restr_->size() > nt);
    pthread_mutex_lock(&mutexsum);
    int val = restr_->at(nt).decrement(hr);
    if((val == -1) && hierch_) {
      base_restr_->decrement(transform(hr, true));
    }
    pthread_mutex_unlock(&mutexsum);
    return val;
  }
  void resampleHyperParams() {
    cerr << "Resampling hyperparameters..\n";
    if(hierch_) base_restr_->resample_prior(mt_genrand_res53);
    for(size_t i=0; i < restr_->size(); ++i) {
      restr_->at(i).resample_prior(mt_genrand_res53);
      //cerr << "a=" << restr_->a() << "\tb=" << restr_->b() << endl;
    }
    resampleTermPenalty();
  }
  double log_prob(const hieroRule& hr, double logp0) {
    int nt = hr.NT;
    pthread_mutex_lock(&mutexsum);
    assert(nt < (int)restr_->size());
    double prb = restr_->at(nt).log_prob(hr, logp0);
    pthread_mutex_unlock(&mutexsum);
    return prb;
  }
  double log_prob_minus1(const hieroRule& hr, double logp0) {
    int nt = hr.NT;
    pthread_mutex_lock(&mutexsum);
    assert(nt < (int)restr_->size());
    double prb = restr_->at(nt).log_prob_minus1(hr, logp0);
    pthread_mutex_unlock(&mutexsum);
    return prb;
  }
  map<int, pair<double, double> > cacheChildPrbs(tree<node_t>::iterator node,
    double* allChildPrbs, bool debug) {
    map<int, pair<double, double> > mchildPrbs;
    *allChildPrbs = 0;
    // get probability of dependent part of sentence derivation
    iterate(node->rule.s, s) {
      if(!isTerm(*s)) {
        tree<node_t>::iterator child = std::find(node.begin(), node.end(), *s);
        assert(child != node.end());
        hieroRule crTr = transform(child->rule, child);
        decrement(crTr);
        p0_t p0 = get_p0(&child->rule, debug);
        double p = log_prob(crTr, p0.second);
        mchildPrbs[*s] = std::make_pair(p, p0.second);
        *allChildPrbs += p;
        increment(crTr, p0);
      }
    }
    return mchildPrbs;
  }
  void printRules(const int iter) {
    std::stringstream fname;
    fname << outputDir_.str() << "/rules." << iter;
    FileHandler fout(fname.str(), std::ios::out, false); 
    for(size_t i=0; i < restr_->size(); ++i) {
      PYP<hieroRule, hieroRule>::iterator rit = restr_->at(i).begin();
      for(; rit != restr_->at(i).end(); ++rit) { // for each rule type
        fout << rit->first << "\t" << rit->second << "\t" 
          << log_prob(rit->first, get_p0(&(rit->first), false).second) << endl;
      }
    }
    fout.close();
  }
  double dataLLAndOtherStuff(const int iter) {
    float avgSrcLen(0), avgTrgLen(0), totTypes(0), badRules(0);
    double dataLL(0);
    for(size_t i=0; i < restr_->size(); ++i) {
      double lgruleProbs(0);
      PYP<hieroRule, hieroRule>::iterator rit = restr_->at(i).begin();
      for(; rit != restr_->at(i).end(); ++rit) { // for each rule type
        double numTables = restr_->at(i).num_tables(rit->first); // get the number of tables it's at
        // calculate its prior times the number of tables
        double lgprb = numTables * get_p0(&(rit->first), false).second; 
        if(lgprb == -std::numeric_limits<double>::infinity()) {
          cout << "WARNING: " << rit->first << " is a BAD rule" << endl;
          ++badRules;
        }
        else lgruleProbs += lgprb;
        int srclen = rit->first.s.size();
        int trglen = rit->first.t.size();
        avgSrcLen += srclen; 
        avgTrgLen += trglen;
      } 
      totTypes += restr_->at(i).num_types();
      dataLL += (lgruleProbs +  restr_->at(i).log_restaurant_prob()); 
      cout << "Rest. " << i << " stats: ";
      restr_->at(i).debug_info(cout);
    }
    cout << "Data LL (iteration " << iter << "): " << dataLL;
    cout << "  (with " << badRules << " -inf rules not counted!)" << endl;
    cout << "Average src/trg lengths: " << (avgSrcLen / totTypes) 
      << " / " << (avgTrgLen/totTypes) << endl; 
    return 1; 
  }
  struct resample_term_penalty {
    vector<PYP<hieroRule, hieroRule> >* restr_;
    bool hierch_;
    PYP<hieroRule, hieroRule>* base_restr_;
    typedef map<int, pair<int, int> > len_t;
    len_t lengths_; //<total terminals, total NTs> in restr for rules of this length
    resample_term_penalty(vector<PYP<hieroRule, hieroRule> >* restr, bool hier, 
      PYP<hieroRule, hieroRule>* base_restr): restr_(restr), hierch_(hier), 
      base_restr_(base_restr) {
      histogram();
    }
    void histogram() {
      // build histogram of src lengths
      for(size_t i=0; i < restr_->size(); ++i) {
        PYP<hieroRule, hieroRule>::iterator rit = restr_->at(i).begin();
        for(; rit != restr_->at(i).end(); ++rit) { // for each rule type
          const hieroRule& r = rit->first;
          int arity = r.arity();
          int noTables = restr_->at(i).num_tables(r); 
          len_t::iterator lit = lengths_.find(r.s.size());
          if(lit == lengths_.end()) {
            lengths_[r.s.size()] = pair<int, int>(0, 0);
            lit = lengths_.find(r.s.size());
          }
          lit->second.first += ((r.s.size() - arity) * noTables); 
          lit->second.second += (arity * noTables); 
        }
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
    resample_term_penalty rtp(restr_, hierch_, base_restr_);
    tp_ = slice_sampler1d(rtp, tp_, mt_genrand_res53, std::numeric_limits<double>::min(), 
      (double) 1.0, (double) 0.0, niterations, 100*niterations);
    cerr << "After resampling tp: " << tp_ << endl;
  }
};
#endif
