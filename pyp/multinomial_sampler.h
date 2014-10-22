#ifndef wa_multinomial_sampler_h 
#define wa_multinomial_sampler_h

#include "multiNTs.h"
#include "sampler_utils.h"

class MultinomialSampler {
private:
  MultiNT* mnt_;
public: 
  MultinomialSampler(void* threadargs) {
    sampleBatchMN(threadargs);
  }
  void sampleBatchMN(void* threadArgs) {
    threadData_t* td = (threadData_t*)threadArgs; 
    mnt_ = td->mnt;
    set<int> srcNTCache, trgNTCache; // thread specific caching
    int lftRuleCache; 
    for(int t = td->sntBegin; t < td->sntEnd; ++t) {
      tree<node_t>::iterator trit = mnt_->trees_[t]->begin();
      bool debug = t == mnt_->snt2debug_;
      while(trit != mnt_->trees_[t]->end()) { // iterate top-down through all legal spans in alignment grid 
        hieroRule prtRule = trit->rule, prtTrsfm = mnt_->transform(prtRule, trit); 
        if((!td->skipDecr) || (trit != mnt_->trees_[t]->begin())) 
          mnt_->decrement(prtTrsfm); // remove current rule from restr configuration
        MultiNT::p0_t prtP0 = mnt_->get_p0(&prtRule, debug);
        double prtProb = mnt_->log_prob(prtTrsfm, prtP0.second) * td->anneal;
        vector<int> symSpans = precomputeSpans(trit, debug); 
        int maxSrcSpan = std::min(trit->send - trit->sstr, td->maxspan);
        for(int srcspan=maxSrcSpan; srcspan > 0; --srcspan) { // top level loop from biggest to smallest 
          for(int srIdxStart=0; srIdxStart < (int)prtRule.s.size(); ++srIdxStart) { // start index of rule lhs 
            bool newSpan(false);
            hieroRule rule;
            bool allNewSrc(true); 
            int sidx=srIdxStart, ss, se;
            int curspan(0);
            bool checkTarget(false);
            do { // build source rule
              rule.s.push_back(prtRule.s[sidx]); //add to the rule
              curspan += (symSpans[sidx+1] - symSpans[sidx]);
              if(curspan == srcspan) {
                ss = symSpans[srIdxStart];
                se = ss + curspan;
                rule.NT = mnt_->NTSymbol(mnt_->src_sents[t][ss], mnt_->src_sents[t][se-1]); //get source side word class
                checkTarget = true; // sample something
                break;
              }
            } while(curspan < srcspan && (++sidx < (int)prtRule.s.size()));
            if(checkTarget) { // build all target side rules
              vector<pair<vector<int>, int> > trgSpans;
              int maxTrgSpan = std::min((int)prtRule.t.size(), td->maxspan);
              for(int trgspan = maxTrgSpan; trgspan > 0; --trgspan) {
                bool allNewTrg(true);
                for(int tIdxStart=0; tIdxStart + trgspan <= (int)prtRule.t.size(); ++tIdxStart) {
                  rule.t.clear();
                  for(int tidx=tIdxStart; tidx-tIdxStart < trgspan; ++tidx) {
                    rule.t.push_back(prtRule.t[tidx]);
                  }
                  if(validRule(rule, prtRule, &allNewSrc, &allNewTrg, 
                    srcNTCache, trgNTCache, &lftRuleCache)) {
                    trgSpans.push_back(make_pair(rule.t, tIdxStart));
                  } // end validRule
                } // end tIdxStart
              } // end trgspan
              if(trgSpans.size()) {
                newSpan = sampleAll(rule, trgSpans, srIdxStart, 
                  prtProb, trit, mnt_->trees_[t], td->anneal, debug, ss, se); 
              }
              if(newSpan) {
                prtRule = trit->rule;
                prtTrsfm = mnt_->transform(prtRule, trit);
                prtP0 = mnt_->get_p0(&prtRule, debug);
                prtProb = mnt_->log_prob(prtTrsfm, prtP0.second) * td->anneal;
                symSpans = precomputeSpans(trit, debug);
                //if(alreadyOn(rule))
                  //--srIdxStart; //keep src start index the same
              }
            }  // end checkTarget
          }  // end srIdxStart 
        }  // end maxSrcSpan
        mnt_->increment(prtTrsfm, prtP0); // add back current rule into restr
        ++trit;
      } // end sentence 
      if(debug) mnt_->printTree(t);
    } // end batch of sentences 
    //pthread_exit((void*) threadArgs);
  }
  bool sampleAll(hieroRule& rule, const vector<pair<vector<int>, int> >& trgSpans, 
    const int sIdx, const double prtSpanPrb, tree<node_t>::iterator node, tree<node_t>* fullTree, 
    const float annealT, const bool debug, const int ss, const int se) {
    bool change(false);
    vector<double> probs;
    map<int, pair<double, double> > mchildPrbs;
    double allChildPrbs(0);
    double currConfigPrb = (prtSpanPrb + allChildPrbs) * annealT; // combine and anneal 
    probs.push_back(currConfigPrb);
    if(debug)
      cout << endl << node->rule << "\t" << probs.back() << endl;
    // for each rule in trgSpan, get probability
    iterate(trgSpans, ts) {
      rule.t = ts->first;
      if(debug) cout << rule << "\t";
      if(alreadyOn(rule)) {
        tree<node_t>::iterator child = std::find(node.begin(), node.end(), rule.s[0]);
        hieroRule ocr = child->rule;
        // combine child with parent
        hieroRule npr = combineChildParent(node->rule, ocr, sIdx, ts->second, debug);
        // get probs for potential rules 
        MultiNT::p0_t nprP0 = mnt_->get_p0(&npr, debug);
        double nprPrb = mnt_->log_prob(mnt_->transform(npr, node), nprP0.second) * annealT; // anneal 
        // get prob of all children minus this one
        double remChildPrbs = allChildPrbs - mchildPrbs[rule.s[0]].first;
        double combPrb = (remChildPrbs + nprPrb) * annealT; // combine and anneal 
        probs.push_back(Log<double>::add(probs.back(), combPrb));
        if(debug) cout << "[" << combPrb << ", ";
      }
      else {
        // split out rule from parent
        int lbl = fullTree->begin()->nodeIndex; // get current label from full tree
        hieroRule npr = newParentRule1(node->rule, rule, sIdx, ts->second, lbl-1); // create new parent rule 
        // get probs for potential rules 
        MultiNT::p0_t nprP0 = mnt_->get_p0(&npr, debug);
        MultiNT::p0_t ncrP0 = mnt_->get_p0(&rule,debug);
        double nprPrb = mnt_->log_prob(mnt_->transform(npr, node, rule.NT), nprP0.second);
        double ncrPrb = mnt_->log_prob(mnt_->transform(rule, node), ncrP0.second); // assume independence for now
        double combPrb = (nprPrb + ncrPrb + allChildPrbs) * annealT; // normalize and anneal 
        probs.push_back(Log<double>::add(probs.back(), combPrb));
        if(debug) cout << "[" << combPrb << ", ";
      }
      if(debug) cout << probs.back() << "]" << endl;
    }
    double rSample = probs.back() + log(dRand());
    vector<double>::iterator pit = std::lower_bound(probs.begin(), probs.end(), rSample);
    int index = int(pit - probs.begin());
    if(debug) cout << "Selected index " << index << endl;
    if(index > 0) { 
      change = true;
      rule.t = trgSpans[index-1].first;
      int tStr = trgSpans[index-1].second;
      if(alreadyOn(rule)) {
        tree<node_t>::iterator child = std::find(node.begin(), node.end(), rule.s[0]);
        hieroRule ocr = child->rule;
        mnt_->decrement(mnt_->transform(ocr, child));
        // combine child with parent
        hieroRule npr = combineChildParent(node->rule, ocr, sIdx, tStr, debug);
        // restructure tree
        node->rule = npr; 
        fullTree->reparent(node, child); // move any children of 'child' to parent 
        fullTree->erase(child); 
      }
      else {
        int& lbl = fullTree->begin()->nodeIndex; // get current label from full tree
        hieroRule npr = newParentRule1(node->rule, rule, sIdx, tStr, lbl-1); // create new parent rule 
        mnt_->updateTree(node, fullTree, npr, rule, --lbl, ss, se);
        mnt_->increment(mnt_->transform(rule, node), mnt_->get_p0(&rule, debug)); // add new child rule to restaurant
      }
    }
    return change;
  }

};
#endif 

