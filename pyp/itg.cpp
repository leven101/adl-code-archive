#include "vocab.h"
#include "file.h"
#include "types.h"
#include "google/dense_hash_map"
/*TODO::
 * - nodes should be stored in hash map for O(1) access
 * - each node should be a struct with pointers to rules it holds (terms and/or nonterms)
 * - each visit to a node should query grammar class for all rules that derive from children
 * - all rule (terms and nonterms) should be in grammar class
 * - handle log probs better
 */
using namespace std;

vector<vector<int> > src_sents, trg_sents;
map<int, map<int, float> > lex_rules; // move this to grammar class eventually
float monotoneRuleProb(0), invertedRuleProb(0);
Vocab src_vocab(false), trg_vocab(false);
const float _ZERO_ = -numeric_limits<float>::infinity();
const size_t maxSentSize = 10;
const int noIterations = 2;
// sentence specific data structures
map<vector<int>, pair<float, float> > inside_probs;
map<vector<int>, pair<float, float> > outside_probs;
// data structure and variables to store sufficient statistics
map<int, map<int, float> > suff_stats; 
float mrp_ss(_ZERO_), irp_ss(_ZERO_);
// viterbi parse
bool viterbi(false);
map<vector<int>, pair<vector<int>, vector<int> > > back_ptrs; // node -> best left child, best right child

// functions
void readSentences(); 
void updateSuffStats(vector<int>&, vector<int>&, int, int, int, int, 
  const float, bool = false);
vector<int> makeNode(int a, int b, int c, int d) {
  vector<int> node(4);
  node[0] = a;
  node[1] = b;
  node[2] = c;
  node[3] = d;
  return node;
}
float getLexRuleProb(vector<int>& src, vector<int>& trg,int s, int t, int u, 
  int v) {
  assert((t-s <= 1) && (v-u <= 1)); // get prob of D_sSUv
  float prb(0);
  int src_id = src[s], trg_id = trg[u];
  if(s == t) { // get probability of X => <null/trg>
    prb = lex_rules[0].find(trg_id)->second;
  }
  else if(u == v) { // get probability of X => <src/null>
    prb = lex_rules[src_id].find(0)->second;
  }
  else {
    prb = lex_rules.find(src_id)->second.find(trg_id)->second;
  }
  if(prb < 1e-08)  // smooth probabilites
    prb = 1e-08;
  return log(prb);
}
bool isLexRule(vector<int>& v) {
  assert(v.size() == 4);
  return (v[1]-v[0] <= 1) && (v[3]-v[2] <= 1);
}
float addLogProbs(float delta1, float delta2 ) {
  if (delta1 > delta2) {
    return( (float) (delta1 + log( 1 + exp( delta2-delta1) ) ) );   
  } else  {
    return( (float) (delta2 + log( 1 + exp( delta1-delta2) ) ) );   
  }
}
float getInvertedInsideNodeProb(vector<int>& src, vector<int>& trg,int s, int t, int u, 
  int v) {
  vector<int> node = makeNode(s,t,u,v);
  vector<int> bp1, bp2; // only used if doing viterbi parse
  map<vector<int>, pair<float, float> >::iterator it;
  it = inside_probs.find(node);
  if((it != inside_probs.end()) && (it->second.second != _ZERO_)) // if node has already been checked
    return it->second.second; // return its inverted probability
  float total_inv_prob(_ZERO_);
  for(int S=s; S<=t; ++S) {
    for(int U=u; U<=v; ++U) {
      int notBothNull = (((S-s)*(t-S)) + ((U-u)*(v-U)));
      if(notBothNull) {
        float p1(0), p2(0);
        // this is a node
        //cerr << "  cell(" << s << t << u << v << ") =" ;
        // these are rules
        //cerr << "  cell(" << s << S << U << v << ") cell(" << S << t << u << U << ")" << endl;
        // get inverted rule probabilities for node sSuU 
        vector<int> child1 = makeNode(s, S, U, v);
        it = inside_probs.find(child1); 
        if(it != inside_probs.end() && (it->second.second != _ZERO_)) {  // if node already visited
          p1 = it->second.second;  // retrieve its inside probability
        }
        else { // process node in 'chart' 
          // 2 choices -- either lexical rule so span is length 1 OR rules spans more than one X
          if(!(isLexRule(child1))) {
            p1 += getInvertedInsideNodeProb(src, trg, s, S, U, v); // accumulate child node inside probs 
          }
          else {
            cerr << "getsInvertedhere1?\n";
            cerr << "cell(" << s << "," << S << "," << U << "," << v << ")" << endl; 
            p1 = getLexRuleProb(src, trg, s, S, U, v);
          }
        }
        vector<int> child2 = makeNode(S, t, u, U);
        it = inside_probs.find(child2); 
        if(it != inside_probs.end() && (it->second.second != _ZERO_)) {  // if node already visited
          p2 = it->second.second;  // retrieve its inside probability
        }
        else { // process node in 'chart' 
          // 2 choices -- either lexical rule so span is length 1 OR rules spans more than one X
          if(!(isLexRule(child2))) {
            p2 += getInvertedInsideNodeProb(src, trg, S, t, u, U); // accumulate child node inside probs 
          }
          else {
            cerr << "getsInvertedhere2?\n";
            cerr << "cell(" << S << "," << t << "," << u << "," << U << ")" << endl; 
            p2 = getLexRuleProb(src, trg, S, t, u, U);
          }
        }
        float ruleProb = (p1) + (p2) + log(invertedRuleProb);
        if(viterbi) {
          if(ruleProb > total_inv_prob) {
            total_inv_prob = ruleProb;
            // store backpointers
            bp1 = child1;
            bp2 = child2;
          }
        }
        else {
          if(total_inv_prob == _ZERO_) 
            total_inv_prob = ruleProb;
          else
            total_inv_prob = addLogProbs(total_inv_prob, ruleProb);
        }
      } // notBothNull
    } // end U
  } // end S
  it = inside_probs.find(node);
  if(it == inside_probs.end()) {
    inside_probs[node].second = total_inv_prob;
    inside_probs[node].first = _ZERO_;
  }
  else {
    assert(it->second.second == _ZERO_);
    it->second.second = total_inv_prob;
  }
  if(viterbi) {
    // if inverted word order is currently more probable than montone replace back pointers
    if(inside_probs[node].second > inside_probs[node].first) {
      back_ptrs[node] = make_pair(bp1, bp2);
    }
  }
  //cerr << "Total inverted inside prob for cell (" << s << t << u << v << ") = " 
    //<< total_inv_prob << " (" << exp(total_inv_prob) << ")" << endl;
  return total_inv_prob;
}
float getInsideNodeProb(vector<int>& src, vector<int>& trg,int s, int t, int u, 
  int v) {
  map<vector<int>, pair<float, float> >::iterator it;
  vector<int> bp1, bp2; // only used if doing viterbi parse
  vector<int> node = makeNode(s,t,u,v);
  it = inside_probs.find(node);
  if(it != inside_probs.end() && (it->second.first != _ZERO_)) 
    return it->second.first; // node already visited so return its monotone probability
  float total_inside_prob(_ZERO_);
  for(int S=s; S<=t; ++S) {
    for(int U=u; U<=v; ++U) {
      int notBothNull = (((S-s)*(t-S)) + ((U-u)*(v-U)));
      if(notBothNull) {
        float p1(0), p2(0);
        // this is a node
        //cerr << "cell(" << s << "," << t << "," << u << "," << v << ") =" ;
        // these are rules
        //cerr << "  cell(" << s << S << u << U << ") cell(" << S << t << U << v << ")" << endl;
        // get monotone rule probabilities for node sSuU 
        vector<int> child1 = makeNode(s, S, u, U);
        it = inside_probs.find(child1); 
        if(it != inside_probs.end() && (it->second.first != _ZERO_)) {  // if node already visited
          p1 = it->second.first;  // retrieve its inside log probability
        }
        else { // process node in 'chart' 
          // 2 choices -- either lexical rule so span is length 1 OR rules spans more than one X
          if(!(isLexRule(child1))) {
            p1 += getInsideNodeProb(src, trg, s, S, u, U); // accumulate child node inside probs 
          }
          else {
            cerr << "getshere1?\n";
            cerr << "cell(" << s << "," << S << "," << u << "," << U << ")" << endl; 
            p1 = getLexRuleProb(src, trg, s, S, u, U);
          }
        }
        vector<int> child2 = makeNode(S, t, U, v);
        it = inside_probs.find(child2); 
        if(it != inside_probs.end() && (it->second.first != _ZERO_)) {
          p2 = it->second.first;
        }
        else {
          if(!(isLexRule(child2))) { // get probability for node D_StUv
            p2 += getInsideNodeProb(src, trg, S, t, U, v); 
          }
          else {
            cerr << "getshere2?\n";
            cerr << "cell(" << S << "," << t << "," << U << "," << v << ")" << endl; 
            p2 = getLexRuleProb(src, trg, S, t, U, v);
          }
        }
        float ruleProb = (p1) + (p2) + log(monotoneRuleProb);
        if(viterbi) {
          if(ruleProb > total_inside_prob) {
            total_inside_prob = ruleProb;
            // store backpointers
            bp1 = child1;
            bp2 = child2;
          }
        }
        else {
          if(total_inside_prob == _ZERO_)
            total_inside_prob = ruleProb;
          else
            total_inside_prob = addLogProbs(total_inside_prob, ruleProb);
        }
      } // end if notBothNull
    } // end U
  } //end S
  it = inside_probs.find(node);
  if(it == inside_probs.end()) {
    inside_probs[node].first = total_inside_prob;
    inside_probs[node].second = _ZERO_;
  }
  else {
    assert(it->second.first == _ZERO_); // should be set to this
    it->second.first = total_inside_prob; 
  }
  if(viterbi) { // store optimal children for this node
    // if montone word order is currently more probable then replace back pointers
    if(inside_probs[node].first > inside_probs[node].second) {
      back_ptrs[node] = make_pair(bp1, bp2);
    }
  }
  //cerr << "Total monotone inside prob for cell (" << s << t << u << v << ") = " 
    //<< total_inside_prob << " (" << exp(total_inside_prob) << ")" << endl;
  return total_inside_prob;
}
float getOutsideNodeProb(vector<int>& src, vector<int>& trg, int s, int t, int u, 
  int v, vector<int>& parent) {
  /* for a given rule X->(X1 X2) used in building the table,
   * the outside probability for each child (of the rule) is
   * updated using the outside probability of its parent together 
   * with the inside probability of its sibling nonterminal
   */
  vector<int> node = makeNode(s,t,u,v);
  // test for top node 
  if(s == 0 && (t == (int)src.size()) && (u == 0) && (v == (int)trg.size())) {
    //outside_probs.insert(pair<vector<int>, pair<float,float> >(node, pair<float, float>(1,1)));
    outside_probs.insert(pair<vector<int>, pair<float,float> >(node, pair<float, float>(0,0))); // log(1) = 0
    parent = node;
  }
  assert(parent.size());
  map<vector<int>, pair<float, float> >::iterator it, oit;
  // get parent node outside probability
  it = outside_probs.find(parent);
  assert(it != outside_probs.end());
  float parentProb = it->second.first;
  // update outside probability for all child nodes 
  for(int S=s; S<=t; ++S) {
    for(int U=u; U<=v; ++U) {
      int notBothNull = (((S-s)*(t-S)) + ((U-u)*(v-U)));
      if(notBothNull) {
        // this is a node
        //cerr << "cell(" << s << t << u << v << ") =" ;
        // these are rules
        //cerr << "  cell(" << s << S << u << U << ") cell(" << S << t << U << v << ")" << endl;
        vector<int> X1 = makeNode(s,S,u,U), X2 = makeNode(S,t,U,v);
        // update first child
        float ruleProb = !(isLexRule(X1)) ? log(monotoneRuleProb) : 
          getLexRuleProb(src, trg, s, S, u, U); // either lexical or nonterminal rule
        it = inside_probs.find(X2); // sibling rule inside log probability
        assert(it != inside_probs.end());
        //float p = parentProb * ruleProb * it->second.first;
        float p = parentProb + ruleProb + it->second.first;
        oit = outside_probs.find(X1);
        if(oit == outside_probs.end()) {
          outside_probs[X1].first = p; // adds node  
          outside_probs[X1].second = _ZERO_;
        }
        else {
          oit->second.first = addLogProbs(oit->second.first, p);
        }
        // update second child
        ruleProb = !(isLexRule(X2)) ? log(monotoneRuleProb) : 
          getLexRuleProb(src, trg, S, t, U, v);
        it = inside_probs.find(X1);
        assert(it != inside_probs.end());
        p = parentProb + ruleProb + it->second.first;
        oit = outside_probs.find(X2);
        if(oit == outside_probs.end()) {
          outside_probs[X2].first = p; // adds node  
          outside_probs[X2].second = _ZERO_;
        }
        else {
          oit->second.first = addLogProbs(oit->second.first, p);
        }
        //outside_probs[X2].first += p;
      }
    }
  }
  return 1;
}
float getInvertedOutsideNodeProb(vector<int>& src, vector<int>& trg, int s, int t, int u, 
  int v, vector<int>& parent) {
  vector<int> node = makeNode(s,t,u,v);
  // test for top node 
  if(s == 0 && (t == (int)src.size()) && (u == 0) && (v == (int)trg.size())) {
    parent = node;
  }
  assert(parent.size());
  map<vector<int>, pair<float, float> >::iterator it, oit;
  // get parent node outside probability
  it = outside_probs.find(parent);
  assert(it != outside_probs.end());
  float parentProb = it->second.second;
  // update outside probability for all child nodes 
  for(int S=s; S<=t; ++S) {
    for(int U=u; U<=v; ++U) {
      int notBothNull = (((S-s)*(t-S)) + ((U-u)*(v-U)));
      if(notBothNull) {
        vector<int> X1 = makeNode(s,S,U,v), X2 = makeNode(S,t,u,U);
        // update first child
        float ruleProb = !(isLexRule(X1)) ? log(invertedRuleProb) : 
          getLexRuleProb(src, trg, s, S, U, v); // either lexical or nonterminal rule
        it = inside_probs.find(X2); // sibling rule inside probability
        assert(it != inside_probs.end());
        //float p = parentProb * ruleProb * it->second.second;
        float p = parentProb + ruleProb + it->second.second;
        oit = outside_probs.find(X1);
        if(oit == outside_probs.end()) {
          outside_probs[X1].second = p; // adds node  
          outside_probs[X1].first = _ZERO_;
        }
        else {
          oit->second.second = addLogProbs(oit->second.second, p);
        }
        //outside_probs[X1].second += p; // add node if not already there
        // update second child
        ruleProb = !(isLexRule(X2)) ? log(invertedRuleProb) : 
           getLexRuleProb(src, trg, S, t, u, U);
        it = inside_probs.find(X1);
        assert(it != inside_probs.end());
        //p = parentProb * ruleProb * it->second.second;
        p = parentProb + ruleProb + it->second.second;
        oit = outside_probs.find(X2);
        if(oit == outside_probs.end()) {
          outside_probs[X2].second = p; // adds node  
          outside_probs[X2].first = _ZERO_; // adds node  
        }
        else {
          oit->second.second = addLogProbs(oit->second.second, p);
        }
        //outside_probs[X2].second += p;
      }
    }
  }
  return 1;
}
void printLexRule(vector<int>& src, vector<int>& trg, vector<int> node) {
  assert(node.size() == 4);
  if(node[0] == node[1]) {
    cout << "X ==> < NULL/" << trg_vocab.getWord(trg[node[2]]) << " >\n"; 
  }
  else if(node[2] == node[3]) {
    cout << "X ==> < " << src_vocab.getWord(src[node[0]]) << "/NULL >\n";
  }
  else {
    cout << "X ==> < " << src_vocab.getWord(src[node[0]]) << "/" << 
      trg_vocab.getWord(trg[node[2]]) << " >\n";
  }
}
//map<vector<int>, pair<vector<int>, vector<int> > > back_ptrs; // node -> best left child, best right child
void viterbiParse(vector<int>& src, vector<int>& trg, vector<int>& parent) {
  map<vector<int>, pair<vector<int>, vector<int> > >::iterator bpit;
  bpit = back_ptrs.find(parent);
  assert(bpit != back_ptrs.end());
  vector<int>& child1 = bpit->second.first;
  vector<int>& child2 = bpit->second.second;
  if(!isLexRule(child1)) {
    viterbiParse(src, trg, child1);
  }
  else {
    printLexRule(src, trg, child1);
  }
  if(!isLexRule(child2)) {
    viterbiParse(src, trg, child2);
  }
  else {
    printLexRule(src, trg, child2);
  }
}
void parseSentence(vector<int>& src, vector<int>& trg) {
  inside_probs.clear(); // clear any old sentence statistics
  outside_probs.clear();
  // Initialize 'chart' 
  //cerr << "INITIAL LEXICAL PROBS\n";
  for(size_t t=0; t < src.size(); ++t) {
    for(size_t v=0; v < trg.size(); ++v) {
      assert(lex_rules[src[t]].find(trg[v]) != lex_rules[src[t]].end());
      vector<int> node = makeNode(t, t+1, v, v+1);
      float prob = getLexRuleProb(src, trg, t, t+1, v, v+1); 
      inside_probs[node].first = prob;
      inside_probs[node].second = prob;
      //if(viterbi)
        //cerr << "cell(" << t << "," << t+1 << "," << v << "," << v+1 << ") = " << prob << endl; 
    }
  }
  for(size_t t=0; t < src.size(); ++t) {
    for(size_t v=0; v <= trg.size(); ++v) {
      assert(lex_rules[src[t]].find(0) != lex_rules[src[t]].end());
      vector<int> node = makeNode(t, t+1, v, v);
      float prob = getLexRuleProb(src, trg, t, t+1, v, v); 
      inside_probs[node].first = prob; 
      inside_probs[node].second = prob;
      //if(viterbi)
        //cerr << "cell(" << t << "," << t+1 << "," << v << "," << v << ") = " << prob << endl; 
    }
  }
  for(size_t t=0; t <= src.size(); ++t) {
    for(size_t v=0; v < trg.size(); ++v) {
      assert(lex_rules[0].find(trg[v]) != lex_rules[0].end());
      vector<int> node = makeNode(t, t, v, v+1);
      float prob = getLexRuleProb(src, trg, t, t, v, v+1); 
      inside_probs[node].first = prob;
      inside_probs[node].second = prob; 
      //if(viterbi)
        //cerr << "cell(" << t << "," << t << "," << v << "," << v+1 << ") = " << prob << endl; 
    }
  }
  // functions that recursively get every inside rule of this node and all children
  //cerr << "INSIDE PROBS\n";
  for(size_t t=1; t <= src.size(); ++t) {
    for(size_t s=0; s < t; ++s) { 
     for(size_t v=1; v <= trg.size(); ++v) {
        for(size_t u=0; u < v; ++u) {
          if((t-s)+(v-u) <= 2) continue;
          getInsideNodeProb(src, trg, s, t, u, v);
          getInvertedInsideNodeProb(src, trg, s, t, u, v);
  } } } } 
  if(viterbi) {
    vector<int> topNode = makeNode(0, src.size(), 0, trg.size());
    viterbiParse(src, trg, topNode);
    cerr << endl;
    back_ptrs.clear();
  }
  else {
    // get outside counts
    vector<int> parentNode;
    //cerr << "OUTSIDE PROBS\n";
    for(size_t t=src.size(); t > 0; --t) {
      for(size_t s=0; s < t; ++s) { 
        for(size_t v=trg.size(); v > 0; --v) {
          for(size_t u=0; u < v; ++u) {
            if(t-s+v-u <= 2) continue;
            getOutsideNodeProb(src, trg, s, t, u, v, parentNode);
            getInvertedOutsideNodeProb(src, trg, s, t, u, v, parentNode);
            parentNode = makeNode(s,t,u,v);
    } } } } // O(n^6)
    //iterate(outside_probs, oit) {
    //  cerr << exp(oit->second.first) << " -- " << exp(oit->second.second) << endl;
    //}
    //cerr << "COLLECT COUNTS\n";
    vector<int> topNode = makeNode(0, src.size(), 0, trg.size());
    const float tip = inside_probs.find(topNode)->second.first;
    for(size_t t=src.size(); t > 0; --t) {
      for(size_t s=0; s < t; ++s) { 
        for(size_t v=trg.size(); v > 0; --v) {
          for(size_t u=0; u < v; ++u) {
            updateSuffStats(src, trg, s, t, u, v, tip);
    } } } }
    for(size_t t=0; t < src.size(); ++t) {
      for(size_t v=0; v < trg.size(); ++v) {
        updateSuffStats(src, trg, t, t+1, v, v+1, tip, true);
      }
    }
    for(size_t t=0; t < src.size(); ++t) {
      for(size_t v=0; v <= trg.size(); ++v) {
        updateSuffStats(src, trg, t, t+1, v, v, tip, true);
      }
    }
    for(size_t t=0; t <= src.size(); ++t) {
      for(size_t v=0; v < trg.size(); ++v) {
        updateSuffStats(src, trg, t, t, v, v+1, tip, true);
      }
    }
  }
}
void updateSuffStats(vector<int>& src, vector<int>& trg, int s, int t, int u, 
  int v, const float totalInsideProb, bool isLexRule) {
  vector<int> parentNode = makeNode(s,t,u,v);
  map<vector<int>, pair<float, float> >::iterator oit, iit1, iit2;
  map<int, map<int, float> >::iterator ssit; 
  map<int, float>::iterator ssit2; 
  // function to stat counts for rules used to build node (s,t,u,v) 
  if(isLexRule) {
    //cerr << "getting lex_prob for cell(" << s << t << u << v << ") : " ;
    float outsideProb(_ZERO_);
    float p(_ZERO_);
    oit = outside_probs.find(parentNode);
    if(oit != outside_probs.end()) {
      outsideProb = oit->second.first; 
    }
    if(outsideProb != _ZERO_) {
      float ruleProb = getLexRuleProb(src, trg, s, t, u, v); 
      p = (ruleProb - totalInsideProb) + outsideProb;
    }
    if(p == _ZERO_) return;
    //cerr << " p = " << p << " -- " << exp(p) << endl;
    if(s == t) { // then null/trg
      //suff_stats[0][trg[u]] += (ruleProb / totalInsideProb) * outsideProb; 
      ssit = suff_stats.find(0); // 0 has already been inserted 
      if(ssit == suff_stats.end()) {
        suff_stats[0][trg[u]] = p;
      }
      else {
        ssit2 = ssit->second.find(trg[u]); 
        if(ssit2 == ssit->second.end()) {
          (ssit->second)[trg[u]] = p; 
        }
        else {
          ssit2->second = addLogProbs(ssit2->second, p); 
        }
      }
    }
    else if(u == v) { // then src/null
      //suff_stats[src[s]][0] += (ruleProb / totalInsideProb) * outsideProb;
      ssit = suff_stats.find(src[s]);
      if(ssit == suff_stats.end()) {  // src[s] not found
        suff_stats[src[s]][0] = p; // add src[s]/null 
      }
      else {
        ssit2 = ssit->second.find(0);
        if(ssit2 == ssit->second.end()) { // null not found
          (ssit->second)[0] = p;
          //ssit->second.insert(make_pair(0, p));
        }
        else {
          ssit2->second = addLogProbs(ssit2->second, p);
        }
      }
    }
    else {
      //suff_stats[src[s]][trg[u]] += (ruleProb / totalInsideProb) * outsideProb;
      ssit = suff_stats.find(src[s]);
      if(ssit == suff_stats.end()) {
        suff_stats[src[s]][trg[u]] = p; 
      }
      else {
        ssit2 = ssit->second.find(trg[u]);
        if(ssit2 == ssit->second.end()) {
          (ssit->second)[trg[u]] = p;
        }
        else {
          ssit2->second = addLogProbs(ssit2->second, p);
        }
      }
    }
  }
  else {
    for(int S=s; S<=t; ++S) {
      for(int U=u; U<=v; ++U) {
        int notBothNull = (((S-s)*(t-S)) + ((U-u)*(v-U)));
        if(notBothNull) {
          oit = outside_probs.find(parentNode);
          float outsideProb1 = oit->second.first; // monotone
          float outsideProb2 = oit->second.second; // inverted
          vector<int> child1 = makeNode(s,S,u,U), child2 = makeNode(S,t,U,v);
          iit1 = inside_probs.find(child1);
          iit2 = inside_probs.find(child2);
          if((outsideProb1 != _ZERO_) && (iit1->second.first != _ZERO_) 
              && (iit2->second.first != _ZERO_)) {
            float p = (log(monotoneRuleProb) - totalInsideProb) +
              (outsideProb1 + iit1->second.first + iit2->second.first);
            if(mrp_ss == _ZERO_)
              mrp_ss = p;
            else
              mrp_ss = addLogProbs(mrp_ss, p);
          }
          if((outsideProb2 != _ZERO_) && (iit1->second.second != _ZERO_) 
              && (iit2->second.second != _ZERO_)) {
            float p = ((log(invertedRuleProb) - totalInsideProb) + 
              (outsideProb2 + iit1->second.second + iit2->second.second));
            if(irp_ss == _ZERO_)
              irp_ss = p;
            else 
              irp_ss = addLogProbs(irp_ss, p); 
          }
          //mrp_ss += ((monotoneRuleProb / totalInsideProb) * 
          //  (outsideProb1 * inside_probs[child1].first * inside_probs[child2].first)); 
          //irp_ss += ((invertedRuleProb / totalInsideProb) * 
          //  (outsideProb2 * inside_probs[child1].second * inside_probs[child2].second)); 
        } // notBothNull
      } // end U
    } // end S
  }
}
void normalizeProbs(bool print) {
  // count total expectations for each source word
  iterate(suff_stats, ssit) {
    float denom(_ZERO_);
    iterate(ssit->second, tit) {
      //cerr << "old probs [" << src_vocab.getWord(ssit->first) << "/" << 
        //trg_vocab.getWord(tit->first) << "] = " << lex_rules[ssit->first][tit->first] << endl;
      if(denom == _ZERO_)
        denom = tit->second;
      else
        denom = addLogProbs(denom, tit->second);
      //cerr << "tit->second = " << tit->second << " / " << exp(tit->second) << endl;
    }
    iterate(ssit->second, tit) {  // normalize 
      lex_rules[ssit->first][tit->first] = exp(tit->second - denom); // convert back into standard probs 
      if(print) {
        cerr << "new probs [" << src_vocab.getWord(ssit->first) << "/" << 
          trg_vocab.getWord(tit->first) << "] = " << lex_rules[ssit->first][tit->first] << endl;
      }
    }
    if(print)
      cerr << endl;
  }
  // normalize nonterminal rules
  mrp_ss = exp(mrp_ss);
  irp_ss = exp(irp_ss);
  monotoneRuleProb = mrp_ss / (mrp_ss + irp_ss); 
  invertedRuleProb = irp_ss / (mrp_ss + irp_ss); 
  if(print) {
    cerr << "monotoneRuleProb = " << monotoneRuleProb << endl;
    cerr << "invertedRuleProb = " << invertedRuleProb << endl;
  }
  mrp_ss = _ZERO_;
  irp_ss = _ZERO_;
}
void normalizeProbs2() {
  /* normalize as joint conditional probability p(s,t|X) */
  float denom(addLogProbs(mrp_ss,irp_ss));
  // count total expectations for all rules
  iterate(suff_stats, ssit) {
    iterate(ssit->second, tit) {
      denom = addLogProbs(denom, tit->second);
    }
  }
  monotoneRuleProb = exp(mrp_ss - denom);
  invertedRuleProb = exp(irp_ss - denom);
  iterate(suff_stats, ssit) {
    iterate(ssit->second, tit) {
      lex_rules[ssit->first][tit->first] = exp(tit->second - denom); // convert back into standard probs 
    }
  }
  mrp_ss = _ZERO_;
  irp_ss = _ZERO_;
}
void EM(int iterations) {
  assert(src_sents.size() == trg_sents.size());
  for(int itr=0; itr < iterations; ++itr) {
    // E-step is inside-outside
    cerr << "Starting iteration " << (itr+1) << endl;
    for(size_t i=0; i < src_sents.size(); ++i) {
      cerr << "Sentence " << (i+1) << " of " << src_sents.size() << endl;
      parseSentence(src_sents[i], trg_sents[i]);
    }
    // M-step
    cerr << "Normalizing parameters for iteration " << (itr+1) << endl;
    //normalizeProbs((itr + 1 == noIterations));
    normalizeProbs2();
    suff_stats.clear(); // clear old statistics
  }
}

//source and target sentences
int main(int argc, char** argv) { 
  time_t start, finish;
  readSentences();
  time(&start);
  cerr << "\n\nINSIDE-OUTSIDE\n\n";
  EM(noIterations);
  cerr << "\n\nVITERBI PARSE\n\n";
  viterbi = true;
  for(size_t i=0; i < 30; ++i) {
    parseSentence(src_sents[i], trg_sents[i]);
  }
  time(&finish);
  cerr << "time taken was " << difftime(finish, start) << " seconds\n";
  return 1;
}

void readSentences() {
  map<int, float> map2copy;
  lex_rules[0] = map2copy; // setup null/target entry
  FileHandler fsrc("source.text", ios::in);
  FileHandler ftrg("target.text", ios::in);
  string line;
  vector<string> srcline, trgline; 
  while(getline(fsrc, line)) {
    vector<int> srclineIDs, trglineIDs; 
    Utils::splitToStr(line, srcline, " ");
    getline(ftrg, line); // get aligned target sentence  
    Utils::splitToStr(line, trgline, " ");
    if(srcline.size() >= maxSentSize || (trgline.size() >= maxSentSize))
      continue;
    iterate(trgline, w) { 
      trglineIDs.push_back(trg_vocab.getWordID(*w));
    }
    iterate(srcline, w) { // for each source word 
      int srcID = src_vocab.getWordID(*w); // add each word to vocab
      srclineIDs.push_back(srcID);
      if(lex_rules.find(srcID) == lex_rules.end()) // add src word to rule list 
        lex_rules[srcID] = map2copy;
      if(lex_rules[srcID].find(0) == lex_rules[srcID].end())
        lex_rules[srcID][0] = 1.0 / (float)srcline.size();
      // for each source word add target rule 
      iterate(trglineIDs, t) { 
        int trgID = *t;
        if(lex_rules[srcID].find(trgID) == lex_rules[srcID].end()) { 
          lex_rules[srcID][trgID] = 1.0 / float(srcline.size());
        }
      if(lex_rules[0].find(trgID) == lex_rules[0].end())
        lex_rules[0][trgID] = 1.0 / (float)trgline.size();
      }
    }
    // add sentence pair to in memory corpus
    src_sents.push_back(srclineIDs);
    trg_sents.push_back(trglineIDs);
    //if(src_sents.size() == 4) break;
  }
  fsrc.close();
  ftrg.close();
  // initialize 2 grammar rules
  monotoneRuleProb = .5;
  invertedRuleProb = .5;
  // test
  /*iterate(lex_rules, i) { 
    cout << "[" <<  i->first << " / " << src_vocab.getWord(i->first) << "]\t";
    iterate(i->second, ii)
      cout << "[" << ii->first << "/" << trg_vocab.getWord(ii->first) << "]\t"; 
    cout << endl;
  }*/
}
