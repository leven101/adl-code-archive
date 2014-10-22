#ifndef wa_block_sampler_h 
#define wa_block_sampler_h

#include <queue>
#include <tr1/unordered_set>
#include "verbose.h" 
#include "grammar.h"
#include "multiNTs.h"
#include "utils.h"
#include "lattice.h"
#include "bottom_up_parser.h"
#include "hg_intersect.h"
#include "hg.h"
#include "sampler.h"
#include "wordid.h"
#include "inside_outside.h"
#include "nt_span.h"
#include "sampler_utils.h"
#include "trule.h"

using std::stringstream;
using std::queue;
using std::tr1::unordered_set;

struct PYPGrammar : public TextGrammar {
  PYPGrammar(const int max_span) : TextGrammar() {
    SetMaxSpan(max_span);
  }
  PYPGrammar(const string& in) : TextGrammar(in) {
    SetMaxSpan(MAX_SENT_LEN());
  }
};
struct TRuleHash {
  size_t operator()(const TRulePtr r) const {
    return hash_value(*r);
  }
  bool operator()(const TRulePtr& a, const TRulePtr& b) const {
    return (*a == *b); 
  }
};
inline bool operator==(const TRulePtr& a, const TRulePtr& b) {
  return (*a == *b); 
}
class GrammarCache {
// cdec uses shared dictionary between source and target 
public: 
  GrammarCache(MultiNT* mnt) {
    mnt_ = mnt;
    ruleDelim_= "||| ";
    cachePYPGrammar();
  }
  ~GrammarCache() {
  }
  void cachePYPGrammar() {
    map<WordID, int> wrd2cnt;
    PYP<hieroRule, hieroRule>::iterator rit = mnt_->restr_->at(0).begin();
    for(; rit != mnt_->restr_->at(0).end(); ++rit) { // for each rule type
      const hieroRule& r = rit->first;
      if((isInsertionRule(r)) || (r.arity() > 2)) continue;
      stringstream rss;
      convertPYPRule(r, rss); // convert to cdec format
      pthread_mutex_lock(&mutexsum);
      TRulePtr rp(new TRule(rss.str())); 
      pthread_mutex_unlock(&mutexsum);
      if(allNTs(r)) {
        ntRuleCache_.insert(rp);
      }
      else {
        iterate(rp->f_, sit) {
          if(isTerm(*sit)) {
            ruleCache_[*sit].insert(rp);
            ++wrd2cnt[*sit];
          }
        }
        iterate(rp->e_, tit) {
          if(isTerm(*tit)) {
            ruleCache_[*tit].insert(rp);
            ++wrd2cnt[*tit];
          }
        }
      }
    }
    // get most frequent terminals
    /*std::multimap<int, WordID> cnt2wrd;
    std::multimap<int, WordID>::reverse_iterator ritr;
    iterate(wrd2cnt, wit)
      cnt2wrd.insert(make_pair(wit->second, wit->first));
    int n(0);
    for(ritr=cnt2wrd.rbegin(); ritr != cnt2wrd.rend(); ++ritr) {
      if(n < 100) {
        freqWords_.insert(ritr->second);
        ++n;
      }
    }*/
  }
  void convertPYPRule(const hieroRule& r, stringstream& ruleStr) {
    // source to target
    ruleStr.clear();
    ruleStr << "[X] " << ruleDelim_;
    for(size_t i=0; i < r.s.size(); ++i) {
      if(!isTerm(r.s[i]))
        ruleStr << (convertNT(r.s[i])) << " ";
      else
        ruleStr << srcVcb_->getWord(r.s[i]) << " ";
    }
    ruleStr << ruleDelim_;
    for(size_t i=0; i < r.t.size(); ++i) {
      if(!isTerm(r.t[i]))
        ruleStr << (convertNT(r.t[i])) << " ";
      else
        ruleStr << trgVcb_->getWord(r.t[i]) << " ";
    }
    ruleStr << ruleDelim_;
    ruleStr << "Feature_1=" << getQProb(r) << endl; // add probability 
  }
  string convertNT(int nt) {
    string nts = "[X,";
    nts += Utils::IntToStr(nt * -1);
    nts += "]";
    return nts;
  }
  double getQProb(const hieroRule& r) {
    // returns proposal distribution. decrements current rule by one 
    MultiNT::p0_t p0 = mnt_->get_p0(&r, false);
    return mnt_->log_prob_minus1(r, p0.second);
  }
  bool allNTs(const hieroRule& r) const {
    return ((r.s.size() == r.t.size()) && (r.s.size() == r.arity()));
  }
  map<WordID, unordered_set<TRulePtr, TRuleHash> > ruleCache_;
  unordered_set<TRulePtr, TRuleHash> ntRuleCache_;
  //vector<unordered_set<hieroRule, hieroRule> > insRules_;
  string ruleDelim_;
  MultiNT* mnt_;
  set<WordID> freqWords_;
};
class BlockSampler {
public:
  BlockSampler(void* threadArgs) {
    threadData_t* td = (threadData_t*)threadArgs; 
    goalSymbol_="X";
    mnt_ = td->mnt;
    prng.reset(new MT19937(iRand()));
    SetSilent(true); // sets cdec to silent 
    grcache_ = td->cache;
    pthread_mutex_lock(&mutexsum);
    cerr << "Block sampling sentence batch " << 
      td->sntBegin << " to " << td->sntEnd << endl;
    pthread_mutex_unlock(&mutexsum);
    //blockSample();
    blockSample(td->sntBegin, td->sntEnd);
    pthread_mutex_lock(&mutexsum);
    cerr << "Done block sampling batch " << 
      td->sntBegin << " to " << td->sntEnd << endl;
    pthread_mutex_unlock(&mutexsum);
    mnt_=NULL;
    grcache_ = NULL;
  }
  ~BlockSampler() {
    //mnt_=0;
  }
  void blockSample() {
    MT19937& rng = *prng;
    grammar_ = new PYPGrammar("/Users/ablev/work/src/pyp/test-block/scfg.txt"); 
    vgrams_.clear();
    GrammarPtr gptr(grammar_);
    vgrams_.push_back(gptr);
    for(size_t s=0; s < 1; ++s) {
      Hypergraph hg;
      currSrcSnt_ = sntToCdecIDs(mnt_->src_sents[s], false); // english
      currTrgSnt_ = sntToCdecIDs(mnt_->trg_sents[s], true);  // chinese
      bool parsed = parseSntPair(s, hg); // run 2-parse algorithm
      if(!parsed) {
        cerr << "Error: sentence " << s << " failed to parse." << endl;
        continue;
      }
      //hg.PrintGraphviz(); 
      vector<unsigned> deriv;
      sampleDerivation(hg, &rng, &deriv); 
    }
  }
  void blockSample(const int str, const int end) {
    MT19937& rng = *prng;
    int failed(0);
    tot_accepted_=0;
    for(size_t s=str; s < end; ++s) {
      if(mnt_->src_sents[s].size() == 1) continue;
      if(mnt_->trg_sents[s].size() == 1) continue;
      if(mnt_->src_sents[s].size() >= 40) continue;
      if(mnt_->trg_sents[s].size() >= 40) continue;
      time_t start, finish;
      time(&start);
      debug_= (mnt_->snt2debug_ == (int)s);
      //debug_ = (s == 499); 
      if(debug_) {
        cerr << "\n\nParsing sentence " << (s+1) << ": ";  
        cerr << "Source length: " << mnt_->src_sents[s].size();
        cerr << "\tTarget length: " << mnt_->trg_sents[s].size() << endl;
        mnt_->printSntPair(s);
      }
      filterPYPGrammar(s);
      Hypergraph hg;
      bool parsed = parseSntPair(s, hg); // run 2-parse algorithm
      if(!parsed) {
        ++failed;
        continue;
      }
      vector<unsigned> deriv;
      sampleDerivation(hg, &rng, &deriv); 
      metropolisHastings(s, hg, deriv);
      time(&finish);
      double secs = difftime(finish, start);
      if(debug_)
        cerr << "Total seconds for succesful parse: " << secs << endl; 
    }
    if(end-str > 1) {
      cerr << "Accepted " << tot_accepted_ << " new samples. (";
      cerr << failed << " sentences failed to parse correctly.)" << endl;
    }
  }
private:
  PYPGrammar* grammar_;
  vector<GrammarPtr> vgrams_;
  MultiNT* mnt_;
  boost::shared_ptr<MT19937> prng;
  set<WordID> srcWrds_, trgWrds_;
  vector<WordID> currSrcSnt_, currTrgSnt_;
  bool debug_;
  int tot_accepted_, tot_zeroes_;
  string goalSymbol_;
  GrammarCache* grcache_;
  vector<WordID> sntToCdecIDs(const vector<wordID_t>& vID, bool target) {
    vector<WordID> cdecIDs(vID.size());
    for (int i=0; i < vID.size(); ++i) {
      if(target)
        cdecIDs[i] = TD::Convert(trgVcb_->getWord(vID[i])); 
      else
        cdecIDs[i] = TD::Convert(srcVcb_->getWord(vID[i])); 
    }
    return cdecIDs;
  }
  void printRawTRule(const TRulePtr trp) const {
    cerr << trp->lhs_ << " ||| ";
    iterate(trp->f_, fit) 
      cerr << *fit << " ";
    cerr << "||| ";
    iterate(trp->e_, eit) 
      cerr << *eit << " ";
    cerr << endl;
  }
  void filterPYPGrammar(const int s) {
    currSrcSnt_ = sntToCdecIDs(mnt_->src_sents[s], false); // foreign 
    currTrgSnt_ = sntToCdecIDs(mnt_->trg_sents[s], true);  // english 
    int grsize(0);
    grammar_ = new PYPGrammar(10000); // new grammar for each sentence
    srcWrds_.clear();
    srcWrds_.insert(currSrcSnt_.begin(), currSrcSnt_.end());
    trgWrds_.clear();
    trgWrds_.insert(currTrgSnt_.begin(), currTrgSnt_.end());
    unordered_set<TRulePtr, TRuleHash> rules2Chk;
    iterate(srcWrds_, sit) {
      rules2Chk.insert(grcache_->ruleCache_[*sit].begin(), grcache_->ruleCache_[*sit].end());
    }
    iterate(trgWrds_, tit) {
      rules2Chk.insert(grcache_->ruleCache_[*tit].begin(), grcache_->ruleCache_[*tit].end());
    }
    rules2Chk.insert(grcache_->ntRuleCache_.begin(), grcache_->ntRuleCache_.end());
    //int fw_added(0);
    iterate(rules2Chk, rit) {
      if(overlap(*rit)) {
        /*if((hasOnlyFreqWords(*rit)) && (!ruleInCurrDerv(*rit, s))) {
          continue;
          //if(++fw_added > 3) continue;
        }
        */
        grammar_->AddRule(*rit);
        ++grsize;
      }
    }
    //grsize += addInsertionRules2Grammar(s);
    vgrams_.clear();
    GrammarPtr gptr(grammar_);
    vgrams_.push_back(gptr);
    if(debug_) 
      cerr << "Added " << grsize << " rules for sentence " << (s+1) << "." << endl;
    if(s % 1000 == 0) { 
      pthread_mutex_lock(&mutexsum);
      cerr << "Added " << grsize << " rules for sentence " << (s+1) << "." << endl;
      pthread_mutex_unlock(&mutexsum);
    }
  }
  bool ruleInCurrDerv(const TRulePtr r, const int snt) const {
    const int sl = r->f_.size();
    const tree<node_t>* td = mnt_->trees_[snt];
    tree<node_t>::iterator it = td->begin(); 
    while(it != td->end()) {
      if(sl == it->rule.s.size()) {
        const hieroRule rc = mnt_->transform(it->rule, false);
        stringstream rss;
        grcache_->convertPYPRule(rc, rss); // convert to cdec format
        pthread_mutex_lock(&mutexsum);
        TRule tr(rss.str()); // does it need mutexes coz rule already seen?!?!
        pthread_mutex_unlock(&mutexsum);
        if(tr == *r) {
          return true;
        }
      }
      ++it;
    }
    return false;
  }
  bool hasOnlyFreqWords(const TRulePtr r) const {
    iterate(r->f_, sit) {
      if((isTerm(*sit)) && 
        (grcache_->freqWords_.find(*sit) == grcache_->freqWords_.end())) {
        return false;
      }
    }
    iterate(r->e_, tit) {
      if((isTerm(*tit)) && 
        (grcache_->freqWords_.find(*tit) == grcache_->freqWords_.end())) {
        return false;
      }
    }
    return true;
  }
  bool overlap(const TRulePtr& r) const {
    // return true if there are terminal matches on both source and target sides for current sentence
    iterate(r->f_, sit) {
      if((isTerm(*sit)) && (srcWrds_.find(*sit) == srcWrds_.end())) {
        return false;
      }
    }
    iterate(r->e_, tit) {
      if((isTerm(*tit)) && (trgWrds_.find(*tit) == trgWrds_.end())) {
        return false;
      }
    }
    return true; 
  }
  bool overlap(const hieroRule& r) const {
    // return true if there are terminal matches on both source and target sides for current sentence
    iterate(r.s, sit) {
      if((isTerm(*sit)) && (srcWrds_.find(*sit) == srcWrds_.end())) {
        return false;
      }
    }
    iterate(r.t, tit) {
      if((isTerm(*tit)) && (trgWrds_.find(*tit) == trgWrds_.end())) {
        return false;
      }
    }
    return true; 
  }
  Lattice sntToCdecLattice(const vector<WordID>& IDs) {
    Lattice lsentence;
    lsentence.resize(IDs.size());
    for (int i=0; i<IDs.size(); ++i) {
      lsentence[i].push_back(LatticeArc(IDs[i], 0.0, 1));  
    }
    return lsentence;
  }
  bool parseSntPair(const int snt, Hypergraph& hg) {
    Lattice lsrc = sntToCdecLattice(currSrcSnt_); 
    ExhaustiveBottomUpParser parser = ExhaustiveBottomUpParser(goalSymbol_, vgrams_);
    if(!parser.Parse(lsrc, &hg)) {
      if(debug_)
        cerr<<"WARNING: source sentence not fully parsed by the grammar!"<<endl;
      return false; 
    }
    //intersect the hg with the target sentence
    Lattice ltrg = sntToCdecLattice(currTrgSnt_);
    if(!HG::Intersect(ltrg, &hg)) {
      if(debug_)
        cerr<<"WARNING: target sentence not fully parsed by the grammar!"<<endl;
      return false; 
    }
    SparseVector<double> reweight;
    reweight.set_value(FD::Convert("Feature_1"), 1.0);
    //reweight.set_value(FD::Convert("PhraseModel_0"), 1.0);
    hg.Reweight(reweight);
    /*cerr << "Final bi-parse: " << endl;
    for(int i=0; i < hg.edges_.size(); ++i) {
      cerr << *hg.edges_[i].rule_ << endl;
    }*/
    return true;
  }
  void sampleDerivation(const Hypergraph& hg, MT19937* rng, vector<unsigned>* sampled_deriv) {
    vector<prob_t> node_probs;
    Inside<prob_t, EdgeProb>(hg, &node_probs);
    queue<unsigned> q;
    q.push(hg.nodes_.size()-1);
    while(!q.empty()) {
      unsigned cur_node_id = q.front();
      //cerr << "NODE=" << cur_node_id << endl;
      q.pop();
      const Hypergraph::Node& node = hg.nodes_[cur_node_id];
      const unsigned num_in_edges = node.in_edges_.size();
      //cerr << "num_in_edges: " << num_in_edges << endl;
      unsigned sampled_edge = 0;
      if (num_in_edges == 1) {
        sampled_edge = node.in_edges_[0];
      } 
      else {
        //prob_t z;
        assert(num_in_edges > 1);
        SampleSet<prob_t> ss;
        for (unsigned j = 0; j < num_in_edges; ++j) {
          const Hypergraph::Edge& edge = hg.edges_[node.in_edges_[j]];
          prob_t p = edge.edge_prob_;
          for (unsigned k = 0; k < edge.tail_nodes_.size(); ++k)
            p *= node_probs[edge.tail_nodes_[k]];
          ss.add(p);
          //z += p;
        }
        //for (unsigned j = 0; j < num_in_edges; ++j) {
        //const Hypergraph::Edge& edge = hg.edges_[node.in_edges_[j]];
        //cerr << exp(log(ss[j] / z)) << " ||| " << edge.rule_->AsString() << endl;
        //}
        //cerr << " --- \n";
        sampled_edge = node.in_edges_[rng->SelectSample(ss)];
      }
      //cerr << "sampled_edge: " << sampled_edge << endl;
      sampled_deriv->push_back(sampled_edge);
      const Hypergraph::Edge& edge = hg.edges_[sampled_edge];
      for (unsigned j = 0; j < edge.tail_nodes_.size(); ++j) {
        q.push(edge.tail_nodes_[j]);
      }
    }
    /*cerr << endl << "**** Sampled Derivation ****" << endl;
    for (unsigned i = 0; i < sampled_deriv->size(); ++i) {
      cerr << *hg.edges_[sampled_deriv->at(i)].rule_ << endl;
    }*/
  }
  double decrementTree(const int snt, unordered_set<hieroRule, hieroRule>& zrules) {
    double prb=0;
    const tree<node_t>* cur_tree = mnt_->trees_[snt];
    // get the current derivation's actual probability 
    tree<node_t>::iterator trit = cur_tree->end();
    while(trit != cur_tree->begin()) {
      --trit;
      const hieroRule rtr = mnt_->transform(trit->rule, false);
      MultiNT::p0_t p0 = mnt_->get_p0(&rtr, false);
      mnt_->decrement(rtr);
      if(mnt_->restr_->at(0).count(rtr)==0) zrules.insert(rtr);
      prb += mnt_->log_prob(rtr, p0.second);
    }
    return prb;
  }
  void incrementTree(const int snt) {
    const tree<node_t>* cur_tree = mnt_->trees_[snt];
    tree<node_t>::iterator trit = cur_tree->begin();
    while(trit != cur_tree->end()) { 
      const hieroRule& r = mnt_->transform(trit->rule, false);
      mnt_->increment(r, mnt_->get_p0(&r, false));
      ++trit;
    }
  }
  void metropolisHastings(const int snt, const Hypergraph& hg, 
    const vector<unsigned>& new_deriv) {
    unordered_set<hieroRule, hieroRule> zrules;
    double smplPropPrb(0), smplActPrb(0);
    double curPropPrb(0), curActPrb = decrementTree(snt, zrules); 

    tree<node_t>* new_tree = convertDerivation(hg, new_deriv); 
    if(debug_) {
      cerr << endl << "**** Current Derivation ****";
      mnt_->printTree(snt);
      cerr << endl << "**** Sampled Derivation ****";
      mnt_->printTree(new_tree);
    }
    const tree<node_t>* cur_tree = mnt_->trees_[snt];
    tree<node_t>::iterator trit; 
    // get the current derivation's proposal probability 
    trit = cur_tree->begin();
    while(trit != cur_tree->end()) {
      const hieroRule rtr = mnt_->transform(trit->rule, false);
      MultiNT::p0_t p0 = mnt_->get_p0(&rtr, false);
      curPropPrb += mnt_->log_prob(rtr, p0.second);
      ++trit;
    }
    // get the sampled derivation's proposal probability
    trit = new_tree->begin();
    while(trit != new_tree->end()) {
      const hieroRule rtr = mnt_->transform(trit->rule, false);
      MultiNT::p0_t p0 = mnt_->get_p0(&rtr, false);
      smplPropPrb += mnt_->log_prob(rtr, p0.second);
      ++trit;
    }
    // get the sampled derivation's actual probability
    trit = new_tree->begin();
    while(trit != new_tree->end()) {
      const hieroRule rtr = mnt_->transform(trit->rule, false);
      MultiNT::p0_t p0 = mnt_->get_p0(&rtr, false);
      smplActPrb += mnt_->log_prob(rtr, p0.second);
      mnt_->increment(rtr, p0);// increment rule 
      if(zrules.find(rtr) != zrules.end()) zrules.erase(rtr);
      ++trit;
    }
    double acpt = (smplActPrb + curPropPrb) - (curActPrb + smplPropPrb); 
    double rand = log(dRand()); 
    if(debug_) {
      cerr << "Current rules: " << cur_tree->size() << endl;
      cerr << "Current Actual Prob: " << curActPrb << "\tCurrent Proposal Prob: " << curPropPrb << endl;
      cerr << "Sample rules: " << new_tree->size() << endl;
      cerr << "Sample Actual Prob: " << smplActPrb << "\tSample Proposal Prob: " << smplPropPrb << endl;
      cerr << "Accept: " << acpt << "\tRand: " << rand << endl;
    }
    if(acpt >= 0 || (rand <= acpt)) {
      if(debug_)
        cerr << "Incorporating new derivation " << endl;
      // replace current derivation 
      delete cur_tree; 
      mnt_->trees_[snt] = new_tree;
      ++tot_accepted_;
      eraseZeroRules(zrules);
    }
    else { //put the rules back in the PYP
      if(debug_)
        cerr << "Keeping the old derivation " << endl;
      trit = new_tree->begin(); // decrement each rule in PYP
      while(trit != new_tree->end()) {
        mnt_->decrement(mnt_->transform(trit->rule, false));
        ++trit;
      }
      delete new_tree;
      incrementTree(snt);
    }
  }
  void eraseZeroRules(const unordered_set<hieroRule, hieroRule>& zrules) {
    iterate(zrules, zit) {
      stringstream rss;
      grcache_->convertPYPRule(*zit, rss); // convert to cdec format
      TRulePtr rp(new TRule(rss.str())); 
      iterate(rp->f_, sit) {
        if(isTerm(*sit)) { 
          grcache_->ruleCache_[*sit].erase(rp);
        }
      }
      iterate(rp->e_, tit) {
        if(isTerm(*tit)) { 
          grcache_->ruleCache_[*tit].erase(rp);
        }
      }
    }
  }
  hieroRule convertCdecRule(const TRule& crule) {
    //cerr << crule << endl;
    hieroRule hr;
    hr.NT = 0; // singleNT_ condition only 
    const vector<WordID>& src = crule.f_;
    const vector<WordID>& trg = crule.e_; 
    // do source side first
    vector<int> ntIdxs;
    //cerr << crule.lhs_ << " ||| ";
    for(int j=0; j < src.size(); ++j) {
      //cerr << src[j] << " ";
      if(isTerm(src[j])) { // terminal 
        wordID_t sw = srcVcb_->getWordID(TD::Convert(src[j]));
        hr.s.push_back(sw);
      }
      else { // non-terminal
        hr.s.push_back(src[j]);
        ntIdxs.push_back(src[j]);
      }
    }
    //cerr << "||| ";
    for(int j=0; j < trg.size(); ++j) {
      //cerr << trg[j] << " ";
      if(isTerm(trg[j])) { // terminal 
        wordID_t tw = trgVcb_->getWordID(TD::Convert(trg[j]));
        hr.t.push_back(tw);
      }
      else { // non-terminal
        int idx = trg[j] * -1;
        hr.t.push_back(ntIdxs[idx]);
      }
    }
    //cerr << endl;
    return mnt_->transform(hr, false);
  }
  tree<node_t>* convertDerivation(const Hypergraph& hg, const vector<unsigned>& deriv) {
    tree<node_t>* tr = new tree<node_t>;
    int nodeIndex(0);
    map<int, int> parent_lbls;
    for(size_t i = 2; i < deriv.size(); ++i) {
      const TRule& crule = *hg.edges_[deriv[i]].rule_;
      //cerr << crule << endl;
      hieroRule hr;
      hr.NT = 0; // singleNT_ condition only 
      const vector<WordID>& src = crule.e_;
      const vector<WordID>& trg = crule.f_; 
      // since inverted after bi-parse, do target side first
      vector<int> ntIdxs;
      //cerr << crule.lhs_ << " ||| ";
      for(int j=0; j < trg.size(); ++j) {
        //cerr << trg[j] << " ";
        if(trg[j] > 0) { // terminal 
          wordID_t tw = trgVcb_->getWordID(TD::Convert(trg[j]));
          hr.t.push_back(tw);
        }
        else { // non-terminal
          hr.t.push_back(trg[j]);
          ntIdxs.push_back(trg[j]);
          parent_lbls[trg[j]] = crule.lhs_;
        }
      }
      //cerr << "||| ";
      for(int j=0; j < src.size(); ++j) {
        //cerr << src[j] << " ";
        if(src[j] > 0) { // terminal 
          wordID_t sw = srcVcb_->getWordID(TD::Convert(src[j]));
          hr.s.push_back(sw);
        }
        else { // non-terminal
          int idx = src[j] * -1;
          hr.s.push_back(ntIdxs[idx]);
        }
      }
      //cerr << endl;
      node_t node;
      node.rule = hr;
      node.label = crule.lhs_; 
      node.sstr = crule.prev_i;
      node.send = crule.prev_j;
      if(node.label < nodeIndex) nodeIndex = node.label;
      if(tr->size() == 0) {
        tr->set_head(node);
      }
      else {
        map<int, int>::iterator prit = parent_lbls.find(node.label);
        assert(prit != parent_lbls.end());
        tree<node_t>::iterator it;
        it = std::find(tr->begin(), tr->end(), prit->second);
        assert(it != tr->end());
        tr->append_child(it, node);
      }
    }
    tr->begin()->nodeIndex = nodeIndex;
    assert(tr->size() == deriv.size()-2);
    return tr; 
  }
};
#endif
