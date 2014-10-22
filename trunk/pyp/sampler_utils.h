#ifndef wa_sampler_h 
#define wa_sampler_h

#include "multiNTs.h"

vector<int> precomputeSpans(tree<node_t>::iterator trit, bool debug) {
  const hieroRule& hr = trit->rule;
  vector<int> symSpans(hr.s.size() + 1, 0);
  tree<node_t>::iterator it; 
  /*if(debug) {
    cout << "PrecomputeSpans" << endl;
    printTreeNode(trit);
  }*/
  for(size_t j=0; j < symSpans.size(); ++j) { // precompute span of each child node in rule 
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
double dRand() {
  return drand48();
  //return rand() * (1/double(RAND_MAX));
}
int iRand(int modval=1) {
  return lrand48() % modval;
  //return rand() % modval;
}
bool isInsertionRule(const hieroRule& X) {
  bool insRule = false;
  if((X.s.size() == 1) && (!isTerm(X.s[0]))) { //source side is just NT 
    if((X.t.size() == 1) && (!isTerm(X.t[0]))) {
      insRule = true;
      //cerr << "WARNING: " << X << " is an insertion rule. Not adding to parsing grammar." << endl;
    }
    else if(X.t.size() > 1) { // target side has more than one symbol
      iterate(X.t, tit) {
        if(isTerm(*tit) != isTerm(X.t[0])) { // mix of terminals and NTs
          //cerr << "WARNING: " << X << " is an insertion rule. Not adding to parsing grammar." << endl;
          insRule = true;
          break;
        }
      }
    }
  }
  return insRule;
}
wordID_t getWordID(string word, Vocab* vcb) {
  if(word == "NULL")
    return Vocab::kOOVWordID;
  else
    return vcb->getWordID(word);
}
#endif
