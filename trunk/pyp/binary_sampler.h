#ifndef wa_binary_sampler_h 
#define wa_binary_sampler_h

#include "multiNTs.h"
#include "sampler_utils.h"


class BinarySampler {
private:
  MultiNT* mnt_;
  bool longSnt_;
public: 
  BinarySampler(void* threadArgs) {
    sampleBatch(threadArgs);
  }
  bool block_allow(const hieroRule& rule) {
    if(!longSnt_) {
      if(isInsertionRule(rule)) return false;
      if(rule.arity() > 2) return false; //DEFICIENT!!!
    }
    return true;
  }
  void sampleBatch(void* threadArgs) {
    threadData_t* td = (threadData_t*)threadArgs; 
    mnt_ = td->mnt; 
    set<int> srcNTCache, trgNTCache; // thread specific caching
    int lftRuleCache; 
    for(int t = td->sntBegin; t < td->sntEnd; ++t) {
      longSnt_ = ((mnt_->src_sents[t].size() >= 40) || (mnt_->trg_sents[t].size() >= 40));
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
            bool newPrtSpan(false);
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
                /*if(debug) {
                  cout << "checking: " << rule << endl;
                  cout << "ss: " << ss << "\tse: " << se << endl;
                }*/
                checkTarget = true; // sample something
                break;
              }
            } while(curspan < srcspan && (++sidx < (int)prtRule.s.size()));
            if((!longSnt_) && (rule.arity() > 2)) checkTarget = false; 
            if(checkTarget) {
              int maxTrgSpan = std::min((int)prtRule.t.size(), td->maxspan);
              for(int trgspan = maxTrgSpan; trgspan > 0; --trgspan) {
                bool allNewTrg(true);
                for(int tIdxStart=0; tIdxStart + trgspan <= (int)prtRule.t.size(); ++tIdxStart) {
                  rule.t.clear();
                  for(int tidx=tIdxStart; tidx-tIdxStart < trgspan; ++tidx) {
                    rule.t.push_back(prtRule.t[tidx]); // build target side of the rule
                  }
                  if(debug) {
                    cout << "Current parent rule : " << prtRule << endl;
                    cout << "Checking validity: " << rule << endl;
                  }
                  if(validRule(rule, prtRule, &allNewSrc, &allNewTrg, 
                    srcNTCache, trgNTCache, &lftRuleCache)) {
                    if(debug) { 
                      cout << "rule is valid" << endl;
                    }
                    if(alreadyOn(rule)) {
                      newPrtSpan = sampleChildOff(rule.s[0], srIdxStart, tIdxStart, prtP0.second,
                        prtProb, trit, mnt_->trees_[t], td->anneal, debug, t);
                      if(newPrtSpan) {
                        --srIdxStart; //keep src start index the same
                      }
                    }
                    else {
                      newPrtSpan = sampleChildOn(rule, srIdxStart, tIdxStart, prtP0.second, 
                        prtProb, trit, mnt_->trees_[t], td->anneal, debug, ss, se);
                    }
                    if(newPrtSpan) { // get new parent's stats
                      prtRule = trit->rule;
                      prtTrsfm = mnt_->transform(prtRule, trit);
                      prtP0 = mnt_->get_p0(&prtRule, debug);
                      prtProb = mnt_->log_prob(prtTrsfm, prtP0.second) * td->anneal;
                      trgspan = -1; // kill outer target loop
                      symSpans = precomputeSpans(trit, debug);
                      break;
                    }
                  } // end validRule
                  else if(debug) cout << "rule is invalid" << endl;
                } // end tIdxStart
              } // end trgspan
            }  // end checkTarget
          }  // end srIdxStart 
        }  // end maxSrcSpan
        mnt_->increment(prtTrsfm, prtP0); // add back current rule into restr
        ++trit;
      } // end sentence 
      if(debug) mnt_->printTree(t);
    } // end batch of sentences 
    mnt_ = NULL; 
  }
  bool sampleChildOff(const int nodeLbl, const int sIdx, const int tIdx, 
    const double prtSpanP0, const double  prtSpanPrb, 
    tree<node_t>::iterator node, tree<node_t>* fullTree, const float annealT,
    const bool debug, int treeNo) {
    bool newSample(false); 
    hieroRule prtSpan = node->rule;
    // find child span to merge
    tree<node_t>::iterator child = std::find(node.begin(), node.end(), nodeLbl);
    assert(child != node.end());
    const hieroRule& ocr = child->rule; // old child rule
    const hieroRule ocrTr = mnt_->transform(ocr, child);
    // have child span so merge with parent span 
    hieroRule npr = prtSpan;
    npr.s.erase(npr.s.begin() + sIdx);
    npr.s.insert(npr.s.begin() + sIdx, ocr.s.begin(),
      ocr.s.end());
    npr.t.erase(npr.t.begin() + tIdx);
    npr.t.insert(npr.t.begin() + tIdx, ocr.t.begin(),
      ocr.t.end());
    if(!block_allow(npr)) return false;
    const hieroRule nprTr = mnt_->transform(npr, node);
    mnt_->decrement(ocrTr);
    // get probs for potential rules 
    MultiNT::p0_t nprP0 = mnt_->get_p0(&npr, debug);
    MultiNT::p0_t ocrP0 = mnt_->get_p0(&ocr, debug);
    double nprPrb = mnt_->log_prob(nprTr, nprP0.second) * annealT; // anneal 
    double ocrPrb = mnt_->log_prob(ocrTr, ocrP0.second); 
    double currRulesProb = (prtSpanPrb + ocrPrb) * annealT; // combine and anneal 
    double log_normalizer = Log<double>::add(currRulesProb, nprPrb);
    double rSample = log_normalizer + log(dRand());
    if(debug) {
      cout << "sampling off...\n";
      cout << "currparentRule : " << prtSpan << " [" << prtSpanP0 << ", " 
        << mnt_->restr_->at(prtSpan.NT).count(mnt_->transform(prtSpan, node)) << ", " << prtSpanPrb << "]" << endl;
      cout << "currchildRule : " << ocr << " [" << ocrP0.second << ", " 
        << mnt_->restr_->at(ocr.NT).count(ocrTr) << ", " << ocrPrb << "]" << endl;
      cout << "current Rules Prob (combined and annealed) : " << currRulesProb << endl;
      cout << "newParentRule : " << npr << " [" << nprP0.second << ", " << 
        mnt_->restr_->at(npr.NT).count(nprTr) << ", " << nprPrb << "]" << endl;
      cout << "normalizer : " << log_normalizer << "\tsampled point : " << rSample << endl;
    }
    if(rSample > currRulesProb) {
      if(debug) {
        cout << " sample accepted\n\n";
      }
      node->rule = npr; // restructure tree
      fullTree->reparent(node, child); // move any children of 'child' to parent 
      fullTree->erase(child); 
      newSample = true;
    }
    else {
      if(debug) {
        cout << " sample rejected\n\n";
      }
      mnt_->increment(ocrTr, ocrP0); // put child rule back into restr
    }
    return newSample;
  }
  bool sampleChildOn(const hieroRule& ncr, const int sIdx, const int tIdx, 
    const double prtSpanP0, const double prtSpanPrb,
    tree<node_t>::iterator node, tree<node_t>* fullTree, const float annealT,
    const bool debug, const int ss, const int se) {
    bool newSample(false);
    hieroRule prtSpan = node->rule;
    int& lbl = fullTree->begin()->nodeIndex; // get current label from full tree
    hieroRule npr = newParentRule1(prtSpan, ncr, sIdx, tIdx, lbl-1); // create new parent rule 
    if(!block_allow(npr)) return newSample;
    if(!block_allow(ncr)) return newSample;
    // get probs for potential rules 
    MultiNT::p0_t nprP0 = mnt_->get_p0(&npr, debug);
    MultiNT::p0_t ncrP0 = mnt_->get_p0(&ncr, debug);
    const hieroRule nprTr(mnt_->transform(npr, node, ncr.NT));
    const hieroRule ncrTr(mnt_->transform(ncr, node));
    double nprPrb = mnt_->log_prob(nprTr, nprP0.second); 
    double ncrPrb = mnt_->log_prob(ncrTr, ncrP0.second); // assume independence for now
    double newRulesProb = (nprPrb + ncrPrb) * annealT; // normalize and anneal 
    double log_normalizer = Log<double>::add(prtSpanPrb, newRulesProb);
    double rSample = log_normalizer + log(dRand());
    if(debug) {
      cout << "sampling on...\n";
      cout << "oldparentRule : " << prtSpan << " [" << prtSpanP0 << ", " 
        << mnt_->restr_->at(prtSpan.NT).count(mnt_->transform(prtSpan, node)) << ", " << prtSpanPrb << "]" << endl;
      cout << "newParentRule : " << npr << " [" << nprP0.second << ", " 
        << mnt_->restr_->at(npr.NT).count(nprTr) << ", " << nprPrb << "]" << endl;
      cout << "newChildRule : " << ncr << " [" << ncrP0.second << ", " 
        << mnt_->restr_->at(ncr.NT).count(ncrTr) << ", " << ncrPrb << "]" << endl;
      cout << "newRulesProb (combined and annealed) : " << newRulesProb << endl;
      cout << "normalizer : " << log_normalizer << "\tsampled point : " << rSample << endl;
    }
    if(rSample > prtSpanPrb) {
      if(debug) {
        cout << " sample accepted\n\n";
      }
      mnt_->updateTree(node, fullTree, npr, ncr, --lbl, ss, se);
      mnt_->increment(ncrTr, ncrP0); // add new child rule to restaurant
      newSample = true;
    }
    else {
      if(debug) {
        cout << " sample rejected\n\n";
      }
    }
    return newSample;
  }
};

#endif
