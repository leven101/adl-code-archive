#include "base.h"

void Base::printTrainTestFeatures() {
  cerr << "Printing all trianing and testing features to file...\n";
  FileHandler trInst("./train.norm", std::ios::out, false);
  iterate(train_labels_, trItr) { 
    //int fno(1);
    int lbl = trItr->second == SPAM ? 1 : 0;
    trInst << lbl << " ";
    vector<float>& feats = data_.at(trItr->first);  
    iterate(feats, iitr) 
      trInst << *iitr << " ";
      //trInst << fno++ << ":" << *iitr << " ";
    trInst << endl;
  }
  trInst.close();
  FileHandler tstInst("./test.norm", std::ios::out, false);
  iterate(test_labels_, tstItr) { 
    int lbl = tstItr->second == SPAM ? 1 : 0;
    tstInst << lbl << " ";
    vector<float>& feats = data_.at(tstItr->first);  
    iterate(feats, iitr) 
      tstInst << *iitr << " ";
    tstInst << endl;
  }
  tstInst.close();
}
void Base::selectFeatSet() {
  numTotNodes_ = features_->linkFeatures_.size();
  assert((int)features_->contentFeatures_.size() == numTotNodes_); 
  int totCFeats = features_->contentFeatures_[0].size(); 
  int totLFeats = features_->linkFeatures_[0].size(); 
  int totHGFeats = graph_->hostFeatures_.size() ? graph_->hostFeatures_[0].size() : 0;
  if(features_->lmFeatures_.size() == 0) {
    for(int i=0; i < numTotNodes_; ++i) {
      features_->streamLMFeatures(i);
    }
  }
  int totLMFeats = features_->lmFeatures_[0].size(); 
  cerr << "Total no. of features = " << totCFeats << " +" 
    << totLFeats << " + " << totHGFeats <<  " + " << totLMFeats << endl;
  // filter/order features based on fisher ratio
  //vector<int> contFeatIdx = fisherDiscriminant(features_->contentFeatures_);
  //vector<int> linkFeatIdx = fisherDiscriminant(features_->linkFeatures_);
  //vector<int> hgFeatIdx = fisherDiscriminant(graph_->hostFeatures_);
  int useContFeats = totCFeats; 
  int useLinkFeats = 0;//totLFeats; 
  int useHGFeats = totHGFeats;
  int useLMFeats = 0; //totLMFeats;
  numTotFeats_ = useLinkFeats + useContFeats + useHGFeats + useLMFeats; 
  cerr << "Number of features used = " << numTotFeats_ << endl;
  assert(data_.empty());
  vector<float> initFeats(numTotFeats_, 0);
  for(int i=0; i < numTotNodes_; ++i) {
    data_.push_back(initFeats);
    for(int f=0; f < useContFeats; ++f) // for selected content features
      data_[i][f] = features_->contentFeatures_.at(i).at(f); //contFeatIdx[f]);
    for(int f=0; f < useLinkFeats; ++f)  // for selected link features
      data_[i][f+useContFeats] = features_->linkFeatures_.at(i).at(f); //linkFeatIdx[f]);
    for(int f=0; f < useHGFeats; ++f)
      data_[i][f+useContFeats+useLinkFeats] = graph_->hostFeatures_.at(i).at(f); //hgFeatIdx[f]); 
    for(int f=0; f < useLMFeats; ++f)
      data_[i][f+useContFeats+useLinkFeats+useHGFeats] = features_->lmFeatures_.at(i).at(f); // 
  }
}
void Base::loadLabels() {
  cerr << "Loading testing and training labels...\n";
  if(!features_->hostsMissingFeats_.size()) 
    features_->loadHostsMissingFeats();
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year +
    "/training/webspam-uk" + year + "-set1-labels.txt"; 
  vector<string> tmpV;
  string line;
  FileHandler fin(infile, std::ios::in);
  while(getline(fin, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int webPageID = atoi(tmpV[0].c_str());
    //if(webPageID >= numTotNodes_) break; 
    // don't add training hosts missing all features
    if(features_->hostsMissingFeats_.find(webPageID) != 
        features_->hostsMissingFeats_.end()) continue;
    //if(features_->smallLMSet.find(webPageID) == features_->smallLMSet.end()) continue;
    Utils::trim(tmpV[1]);
    if(tmpV[1] == "undecided") continue;
    LABEL_T label = tmpV[1] == "spam" ? SPAM : (tmpV[1] == "nonspam" ? NONSPAM : UNKNOWN); 
    train_labels_[webPageID] = label;
  }
  fin.close();
  // load test set
  infile = params_->getParamValue("top-dir") + year + 
    "/testing/webspam-uk" + year + "-set2-labels.txt"; 
  FileHandler fin2(infile, std::ios::in);
  while(getline(fin2, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int webPageID = atoi(tmpV[0].c_str());
    //if(webPageID >= numTotNodes_) break; 
    if(features_->hostsMissingFeats_.find(webPageID) != 
        features_->hostsMissingFeats_.end()) {
      continue;
    }
    //if(features_->smallLMSet.find(webPageID) == features_->smallLMSet.end()) continue;
    Utils::trim(tmpV[1]);
    if(tmpV[1] == "undecided") continue;
    //LABEL_T label = tmpV[1] == "spam" ? SPAM : (tmpV[1] == "nonspam" ? NONSPAM : UNKNOWN); 
    //test_labels_[webPageID] = label;
    test_labels_[webPageID] = UNKNOWN; // do this now to make sure not training on test;
  }
  fin2.close();
  cerr << "training Istances = " << train_labels_.size() << endl;
  cerr << "testing Istances = " << test_labels_.size() << endl;
}
void Base::loadTestLabels() {
  vector<string> tmpV;
  string line, year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year + 
    "/testing/webspam-uk" + year + "-set2-labels.txt"; 
  map<int, LABEL_T>::iterator tit;
  FileHandler fin2(infile, std::ios::in);
  while(getline(fin2, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int webPageID = atoi(tmpV[0].c_str());
    tit = test_labels_.find(webPageID); 
    if(tit != test_labels_.end()) { 
      Utils::trim(tmpV[1]);
      LABEL_T label = tmpV[1] == "spam" ? SPAM : (tmpV[1] == "nonspam" ? NONSPAM : UNKNOWN); 
      tit->second = label;
    }
  }
  fin2.close();
}
vector<int> Base::fisherDiscriminant(vector<vector<float> >& X) {
  map<int, LABEL_T> labeled_feats_ = train_labels_;
  assert(X.size() >= labeled_feats_.size());
  int noFeats = X[0].size();
  vector<float> FR(noFeats);
  float noSpamLbls(0), noNonspamLbls(0), minstds = 1;
  float sMu, nsMu; // holds means for spam/nonspam
  iterate(labeled_feats_, trItr) { // for each training website 
    LABEL_T label = trItr->second; 
    if(label == SPAM) 
      ++noSpamLbls;  // collect label counts
    else if(label == NONSPAM) 
      ++noNonspamLbls;
  }
  //cerr << noSpamLbls << " spam labels and " << noNonspamLbls << " nonspam labels\n";
  assert(noSpamLbls + noNonspamLbls == labeled_feats_.size());
  for(int f=0; f < noFeats; ++f) {
    float totSpam(0), totNonspam(0); 
    iterate(labeled_feats_, trItr) { // this loops gets mean of each feature
      float val = X[trItr->first][f]; 
      LABEL_T label = trItr->second; 
      if(label == SPAM) {
        totSpam += val; 
      }
      else if(label == NONSPAM) {
        totNonspam += val; 
      }
    }
    sMu = totSpam / noSpamLbls;
    nsMu = totNonspam / noNonspamLbls;
    //cerr << "Feature " << f << " has a average of " << sMu; 
    //cerr << " for spam and " << nsMu << " for nonspam\n"; 
    // compute sample-based variance for each feature
    totSpam = 0; 
    totNonspam = 0;
    iterate(labeled_feats_, trItr) {
      float val = X[trItr->first][f]; 
      LABEL_T label = trItr->second; 
      if(label == SPAM) 
        totSpam += pow(val - sMu, 2);
      else if(label == NONSPAM) 
        totNonspam += pow(val - nsMu, 2);
    } 
    float std2 = (totSpam / noSpamLbls) + (totNonspam / noNonspamLbls);
    // calc Fisher Discriminant for this feature
    if(std2 == 0) {
      if(sMu != nsMu) 
        FR[f] = -pow(sMu - nsMu, 2); //in case of the stds is zero but the average not same, mark it with negative value to deal with later
      else 
        FR[f] = 0;  //if stds and average difference are both zero, set Fisher Ratio to 0
    }
    else {
      FR[f] = pow(sMu - nsMu, 2) / std2;  //calculate Fisher Ratio value
      if(std2 < minstds) 
        minstds = std2;  //renew the minimum non-zero std value, if poissable
    }
    //sMu = abs(sMu);
    //nsMu = abs(nsMu);
    //if((sMu < 1e-6) && (nsMu < 1e-6))
      //FR[f] = 0;
  } // end of feature loop
  vector<int> featIdx(noFeats);
  for(int f=0; f < noFeats; ++f) {
    featIdx[f] = f;
    if(FR[f] < 0)  //to avoid same very small noise signal's effect,set those features Fisher Ratio value zero
      FR[f] = -FR[f] / minstds;
    //cerr << "FR[" << f << "] = " << FR[f] << endl;
  }

  double xp;//xp is a temporary variable for Fisher Ratio exchange
  int xi;   //xi is a temporary variable for feature index exchange
  //sort feature index based on Fisher Ratio
  for(size_t ii=0;ii<FR.size();ii++) {
    for (size_t jj=ii;jj<FR.size();jj++) {
      if(FR[jj]>FR[ii]){
        xp=FR[ii];
        xi=featIdx[ii];
        FR[ii]=FR[jj];
        featIdx[ii]=featIdx[jj];
        FR[jj]=xp;
        featIdx[jj]=xi;
      }
    }
  }
  /*for(int i=0; i < (int)FR.size(); ++i)
    cerr << i << " " << featIdx[i] << " " << FR[i] << endl;
  cerr << "*********************************\n\n";
  */
  return featIdx;
}

