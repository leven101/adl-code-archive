#include "naiveBayes.h"
void NaiveBayes::train() {
  iterate(train_set_, trItr) { // for each training website 
    LABEL_T label = trItr->second; 
    if(binary_) binaryHistograms(trItr->first, label);
    else normalizedHistograms(trItr->first, label);
  }
  cerr << "total spam labels=" << noTrSpamLbls_ << "\ttotal nonspam labels=" << noTrNonspamLbls_ << endl;
}
void NaiveBayes::normalizedHistograms(int hostID, LABEL_T label) {
  if(label == SPAM) ++noTrSpamLbls_;  // collect label counts for priors
  else if(label == NONSPAM) ++noTrNonspamLbls_;
  // for each feature, count value given label 
  for(size_t f=0; f < data_.at(hostID).size(); ++f) { // for each feature
    float value = data_[hostID][f];
    if(label == SPAM) ++spamFeats_[f][value];
    else if(label == NONSPAM) ++nonspamFeats_[f][value];
  }
  newstuff_ = true;
}
void NaiveBayes::computeStats() {
  for(size_t f=0; f < spamFeats_.size(); ++f) {
    float samples(0), value(0);
    iterate(spamFeats_[f], itr) {
      samples += itr->second;
      value += itr->first * itr->second;
    }
    float mu = value / samples;
    value = 0;
    iterate(spamFeats_[f], itr) {
      value += (pow(itr->first - mu, 2) * itr->second); 
    }
    float sigma = value / samples;
    spamStats_[f] = pair<float,float>(mu, sigma);
  }
  for(size_t f=0; f < nonspamFeats_.size(); ++f) {
    float samples(0), value(0);
    iterate(nonspamFeats_[f], itr) {
      samples += itr->second;
      value += itr->first * itr->second;
    }
    float mu = value / samples;
    value = 0;
    iterate(nonspamFeats_[f], itr) {
      value += (pow(itr->first - mu, 2) * itr->second); 
    }
    float sigma = value / samples;
    nonspamStats_[f] = pair<float,float>(mu, sigma);
  }
}
void NaiveBayes::binaryHistograms(int id, LABEL_T label) {
  if(label == SPAM) ++noTrSpamLbls_;  // collect label counts for priors
  else if(label == NONSPAM) ++noTrNonspamLbls_;
  // binary features - did feature occur given label 
  for(size_t f=0; f < data_.at(id).size(); ++f) { // for each feature
    if(data_[id][f] != 0) { // if feature occurs then count it given class
      if(label == SPAM) ++binSpamFeats_[f];
      else if(label == NONSPAM) ++binNonspamFeats_[f];
    }
  }
}
LABEL_T NaiveBayes::predict(int node, float* ps, float* pns, bool doit) {
  // get priors
  //float noExamples = noTrSpamLbls_ + noTrNonspamLbls_,
  //  priorSpam(noTrSpamLbls_ / noExamples), 
  //  priorNonspam(noTrNonspamLbls_ / noExamples);
  float priorSpam(0.1), 
    priorNonspam(0.9);
  if(binary_) {
    binaryPredict(node, ps, pns); 
  }
  else {
    //normalPredict(node, ps, pns, doit); 
    gaussPredict(node, ps, pns);
  }
  *ps += log(priorSpam);
  *pns += log(priorNonspam);
  if(doit) cerr << "ps = " << *ps << "\tpns = " << *pns << endl;
  return UNKNOWN; 
}
void NaiveBayes::gaussPredict(int nodeID, float* ps, float* pns) {
  *ps = 0;
  *pns = 0;
  if(newstuff_) {
    computeStats();
    newstuff_ = false;
  }
  map<float, float>::iterator m_it;
  for(size_t idx=0; idx < data_[nodeID].size(); ++idx) { // for each feature
    float probSpam(0), probNonspam(0);
    float featVal = data_[nodeID][idx];
    probSpam = gaussian(spamStats_[idx].first, spamStats_[idx].second, featVal); 
    probNonspam = gaussian(nonspamStats_[idx].first, nonspamStats_[idx].second, featVal); 
    *ps += log(probSpam);
    *pns += log(probNonspam);
  }
}
float NaiveBayes::normalPredict(int testPage, float* ps, float* pns, bool doit) {
  if(doit) cerr << "spam labels = " << noTrSpamLbls_ << "\t nonspam lables = " << noTrNonspamLbls_ << endl;
  *ps = 0;
  *pns = 0;
  map<float, float>::iterator m_it;
  for(size_t idx=0; idx < data_[testPage].size(); ++idx) { // for each feature
    float probSpam(0), probNonspam(0);
    if(doit) { 
      cerr << "Num. of spamFeats = " << spamFeats_[idx].size();
      cerr << "\t Num. of nonspamFeats = " << nonspamFeats_[idx].size() << endl;
    }
    float featVal = data_[testPage][idx];
    m_it = spamFeats_[idx].find(featVal); 
    if(m_it != spamFeats_[idx].end()) { // if this value was in training data with SPAM label
      probSpam = m_it->second;
    }
    if(doit) cerr << "feature " << idx << ": count spamFeats_= " << probSpam; 
    probSpam /= noTrSpamLbls_;
    if(doit) cerr << "  prob = " << probSpam;
    m_it = nonspamFeats_[idx].find(featVal); 
    if(m_it != nonspamFeats_[idx].end()) { // if this value was in training data with NONSPAM label
      probNonspam = m_it->second;
    }
    if(doit) cerr << "\t count nonspamFeats " << probNonspam; 
    probNonspam /= noTrNonspamLbls_; 
    if(doit) cerr << "  prob = " << probNonspam << endl;
    if(probSpam > 0 && probNonspam > 0) {
      *ps += log(probSpam);
      *pns += log(probNonspam);
    }
    else if(probSpam > 0 && probNonspam == 0) {
      *ps += log(probSpam);
      *pns += log(probSpam) - 3;
    }
    else if(probSpam == 0 && probNonspam > 0) {
      *ps += log(probNonspam) - 3;
      *pns += log(probNonspam);
    }
  } // end
  return 0;
}
float NaiveBayes::binaryPredict(int id, float* ps, float* pns) {
  *ps = 0;
  *pns = 0;
  for(size_t idx=0; idx < data_[id].size(); ++idx) {// get prob of each feature occurring
    float probSpam(0), probNonspam(0);
    if(data_[id][idx] != 0) {  
      probSpam = (float)binSpamFeats_[idx];
      probNonspam = (float)binNonspamFeats_[idx];
    }
    else { // or prob of it not occurring
      probSpam = noTrSpamLbls_ - (float)binSpamFeats_[idx];
      probNonspam = noTrNonspamLbls_ - (float)binNonspamFeats_[idx];
    }
    probSpam /= noTrSpamLbls_;
    probNonspam /= noTrNonspamLbls_;
    *ps += probSpam > 0 ? log(probSpam) : 0;
    *pns += probNonspam > 0 ? log(probNonspam) : 0;
  } // end for each feature
  return 0;
}
