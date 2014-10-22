#include "hmm.h"
#include "params.h"
#include "model1.cpp"
int curSent;
vector<float> LL_;
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"source", "./srctoy", "src", Parameters::kStringValue, ""},
  {"target", "./trgtoy", "trg", Parameters::kStringValue, ""},
  {"EMIterations", "1", "em", Parameters::kIntValue, ""},
};

int main(int argc, char**argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  HMM hmm(params.getParam("source"), params.getParam("target"));
  hmm.EM(atoi(params.getParam("EMIterations").c_str()));
  return 1;
}
/* notes - for alpha probs you include current trellis cell.
 *         for beta probs you don't include current cell */
void HMM::clearCounts() {
  counts_->clear();
  totTrgPrb_.clear();
  if(alignCnts_) delete alignCnts_;
  alignCnts_ = new vector<ttable_t>(alignPrb_->size());
}
int HMM::EM(int iterations) {
  // runs the Baum-Welch Algorithm 
  assert(src_sents.size() == trg_sents.size());
  std::vector<wordID_t> *srcSnt, *trgSnt;
  int noSents = src_sents.size();
  // for each iteration
  for(int iter = 1; iter <= iterations; ++iter) {
    clearCounts();
    LL_.push_back(0);
    for(int snt = 0; snt < noSents; ++snt) { // for all sentences
      curSent = snt + 1;
      cerr << "NEW sentence number: " << (snt + 1) << endl;
      trgSnt = &trg_sents[snt]; int trgLen = trgSnt->size();
      srcSnt = &src_sents[snt]; int srcLen = srcSnt->size(); 
      newTrellis(srcLen, trgLen); // instantitate a new trellis
      // run the forward-backward algorithm and gather stats
      float sent_prb = forward(*srcSnt, *trgSnt);
      LL_.back() += log(sent_prb) * (curSent==1?10:20);
      backward(*srcSnt, *trgSnt);
      count(*srcSnt, *trgSnt, sent_prb);
    } // end all sentences
    // update alignment and translation probabilities
    normalize();
  } // end EM iteration
  for(int ll = 0; ll < LL_.size(); ++ll)
    cerr << "Log Likelihood " << ll << " = " << LL_[ll] << endl;
  return 1;
}
float HMM::forward(vector<wordID_t>& srcSnt, vector<wordID_t>& trgSnt) {
  assert(alpha_ && beta_);
  //translation prob * alignment prob * sum_(rows)(trellis_[j-1])
  int trgLen = trgSnt.size();
  int srcLen = srcSnt.size(); 
  ttable_t::iterator trprb, alprb;
  for(int j = 0; j < srcLen; ++j) {  // trellis columns
    for(int i = 0; i < trgLen; ++i) {  // trellis rows
      cerr << "alphatrellis["<<i<<"]["<<j<<"] = " <<trgSnt[i] <<","<<srcSnt[j] ;
      trprb = wrdProbs_->find(wrdPair_t(srcSnt[j], trgSnt[i]));
      assert(trprb != wrdProbs_->end());
      if(j == 0) { // start case
        alprb = (*alignPrb_)[0].find(wrdPair_t(0,trgSnt[i])); 
        assert(alprb != (*alignPrb_)[0].end());
        alpha_[i][j] = alprb->second * trprb->second;
      }
      else {  // current node is sum of all possible prior states
        for(int row = 0; row < trgLen; ++row) {
          // state of [row][col] before going to this [row][col]
          alprb = (*alignPrb_)[0].find(wrdPair_t(trgSnt[row], trgSnt[i]));
          assert(alprb != (*alignPrb_)[0].end());
          alpha_[i][j] += alpha_[row][j-1] * // previous node
            alprb->second * trprb->second; // jump prb, translation prb 
        }
      }
      cerr << " = " << alpha_[i][j] << endl;
    }
  }
  // get full sentence probability
  float sent_prb(0);
  for(int i = 0; i < trgLen; ++i) sent_prb += alpha_[i][srcLen-1];
  return sent_prb;
}
float HMM::backward(vector<wordID_t>& srcSnt, vector<wordID_t>& trgSnt) {
  assert(alpha_ && beta_);
  int trgLen = trgSnt.size();
  int srcLen = srcSnt.size(); 
  ttable_t::iterator trprb, alprb;
  for(int i = 0; i < trgLen; ++i) beta_[i][srcLen-1] = 1.0; // start case
  for(int j = srcLen-2; j >= 0; --j) {  // trellis columns (start from end)
    for(int i = 0; i < trgLen; ++i) {  // trellis rows
      cerr << "betatrellis["<<i<<"]["<<j<<"] = " <<trgSnt[i] <<","<<srcSnt[j] ;
      for(int row = 0; row < trgLen; ++row) { // for each state
        // need prob of jump state producing next word 
        trprb = wrdProbs_->find(wrdPair_t(srcSnt[j+1], trgSnt[row]));
        assert(trprb != wrdProbs_->end());
        // get probability of this state transitioning to next col state
        alprb = (*alignPrb_)[0].find(wrdPair_t(trgSnt[i], trgSnt[row]));
        assert(alprb != (*alignPrb_)[0].end());
        beta_[i][j] += beta_[row][j+1] * // next node 
          alprb->second * trprb->second;
      }
      cerr << " = " << beta_[i][j] << endl;
    }
  } 
  return 0;
}
/* count() is inefficient since we can easily move functionality in count() to backward() */
float HMM::count(vector<wordID_t>& srcSnt, vector<wordID_t>& trgSnt, const float sent_prb) {
  int trgLen = trgSnt.size();
  int srcLen = srcSnt.size(); 
  ttable_t::iterator trprb, alprb;
  vector<float> gamma, delta;
  for(int j = 0; j < srcLen; ++j) { 
    for(int i = 0; i < trgLen; ++i) {  // for each generating state
      float a = alpha_[i][j]; // get alpha prob for state
      delta.push_back(0.0); 
      if(j < srcLen - 1) { // for all but last column
        for(int row = 0; row < trgLen; ++row) { // to every other state
          // need prob of jump state producing next word 
          trprb = wrdProbs_->find(wrdPair_t(srcSnt[j+1], trgSnt[row]));
          assert(trprb != wrdProbs_->end());
          // get probability of state transition
          alprb = (*alignPrb_)[0].find(wrdPair_t(trgSnt[i], trgSnt[row]));
          assert(alprb != (*alignPrb_)[0].end());
          // need backward probability 
          float b = (a * beta_[row][j+1] * alprb->second * trprb->second);
          gamma.push_back(b / sent_prb);
          if(curSent == 1) gamma.back() *= 10;
          else gamma.back() *= 20;
          delta.back() += gamma.back();
          // increment all alignment counts
          (*alignCnts_)[0][wrdPair_t(trgSnt[i], trgSnt[row])] += gamma.back();
        }
      }
      else delta.back() = a / sent_prb; 
      // increment prob counts  (add if not in already)
      (*counts_)[wrdPair_t(srcSnt[j], trgSnt[i])] += delta.back();  // * number of matching sentences in corpus
      totTrgPrb_[trgSnt[i]] += delta.back();
      if(j == 0) { // increment start state probs for 0->current trg ID
        (*alignCnts_)[0][wrdPair_t(0, trgSnt[i])] += delta.back();
      }
    }
  }
  assert(int(delta.size()) == srcLen*trgLen);
  cerr << "Gamma:\n";
  iterate(gamma, itr) cerr << *itr << endl;
  cerr << "Delta:\n";
  iterate(delta, itr) cerr << *itr << endl;
  cerr << "totTrgPrb_:\n";
  iterate(totTrgPrb_, itr) 
    cerr << itr->first << " = " << itr->second << endl;
  return 0;
}
float HMM::normalize() {
  int srcVcbSize = srcVcb->size(); 
  int trgVcbSize = trgVcb->size();
  cerr << "srcVcbSize = " << srcVcbSize << "\ttrgVcbSize = " << trgVcbSize << endl;
  for(int sID = 1; sID <= srcVcbSize; ++sID) 
    for(int tID = 1; tID <= trgVcbSize; ++tID) {
      wrdPair_t wp(sID, tID);
      ttable_t::iterator trprb = wrdProbs_->find(wp);
      // update probabilities with counts 
      if(trprb != wrdProbs_->end()) {
        ttable_t::iterator trcnt = counts_->find(wp);
        assert(trcnt != counts_->end()); //????
        trprb->second = trcnt->second / totTrgPrb_[tID];  
      }
    }
  cerr << "Probs:\n";
  piterate(wrdProbs_, li)
    cerr << srcVcb->getWord(li->first.first) << "|" << trgVcb->getWord(li->first.second) << "=" << li->second << endl;
  // normalize transition probs 
  for(int l = 0; l <= 0/*maxSntLen_*/; ++l) {
    if((!alignCnts_->at(l).size())||(!alignPrb_->at(l).size())) continue;
    std::map<count_t, float> alTots;
    // get normalization denominator
    iterate(alignCnts_->at(l), li) { 
      cerr << "li-first = " << li->first.first << "/" << li->first.second << " " << li->second << endl;
      alTots[li->first.first] += li->second;
    }
    cerr << "align Totals:\n";
    iterate(alTots, li) 
      cerr << li->first << " " << li->second << endl;
    iterate(alignPrb_->at(l), ap) {
      ttable_t::iterator li = alignCnts_->at(l).find(ap->first);
      assert(li != alignCnts_->at(l).end());
      ap->second = li->second / alTots[li->first.first];
    }
  }
  cerr << "align Probs:\n";
  iterate(alignPrb_->at(0), li) 
    cerr << trgVcb->getWord(li->first.first) << "|" << trgVcb->getWord(li->first.second) << "=" << li->second << endl;
  return 0;
}
void HMM::newTrellis(const int cols, const int rows) {
  // instantiates a new row x col trellis
  if(alpha_ || beta_) freeTrellis(rows); 
  alpha_ = new float*[rows];
  beta_ = new float*[rows];
  for(int i = 0; i < rows; ++i) {
    alpha_[i] = new float[cols];
    std::fill(alpha_[i], alpha_[i] + cols, 0);
    beta_[i] = new float[cols];
    std::fill(beta_[i], beta_[i] + cols, 0);
  }
}
void HMM::freeTrellis(const int rows) {
  for(int i=0; i < rows; ++i) {
    delete[] alpha_[i];
    delete[] beta_[i];
  }
  delete[] alpha_;
  delete[] beta_;
  alpha_ = 0;
  beta_ = 0;
}
