#include "goldlabels.h" 

GoldLabels::GoldLabels(Parameters& p) {
  params_ = &p;
  noTestEp_ = atoi(params_->getParam("test-epochs").c_str());
  ss_ = params_->getParam("no-ss") == Parameters::kFalseValue;
  load();
  buildClasses();
  if(ss_) {
    subsample();
  }
  else {
    lbls_ = org_;
  }
  noClasses_ = TS::unique(&lbls_[0], lbls_.size());
  // delay trends, ie, 'test tomorrow's trends'
  delay_ = atoi(params_->getParam("delay").c_str());
  cerr << "Delaying gold labels by " << delay_ << " epochs" << endl;
  lbls_.erase(lbls_.begin(), lbls_.begin()+delay_);
  stats();
}
void GoldLabels::load() {
  // gold labels are from last column in file 
  FileHandler fin(params_->getParam("ts-file"), std::ios::in);
  string line;
  int curline(0), skip = atoi(params_->getParam("skip").c_str());
  vector<string> vline;
  const int tsidx = atoi(params_->getParam("ts-index").c_str());
  const bool open = params_->getParam("open") == Parameters::kTrueValue;
  while(getline(fin, line)) {
    if(++curline <= skip) continue; 
    Utils::splitToStr(line, vline, " ");
    assert(tsidx < (int)vline.size());
    string feats = tsidx > -1 ? vline[tsidx] : vline.back();
    if(feats.find("/") == string::npos) {
      org_.push_back(atof(feats.c_str()));
    }
    else { // contains multiple features
      Utils::splitToStr(feats, vline, "/"); // split features
      org_.push_back(atof(vline[(open?0:1)].c_str()));
    }
  }
  fin.close();
  cerr << "Loaded " << org_.size() << " gold labels." << endl;
}
void GoldLabels::subsample() {
  cerr << "Subsampling data..." << endl;
  if(noClasses_ > 2) {
    cerr << "WARNING: Unable to subsample data with more than 2 classes\n";
    return;
  }
  srand48(123);
  size_t zeros = std::count(org_.begin(), org_.end(), 0); 
  vector<int> oneIdxs;
  for(size_t i=0; i < org_.size(); ++i) {
    if(org_[i] == 1) oneIdxs.push_back(i);
  }
  std::random_shuffle(oneIdxs.begin(), oneIdxs.end(), Utils::iRand);
  int todel = oneIdxs.size() - zeros; 
  for(int i=0; i < todel-noTestEp_; ++i) {
    int ridx = Utils::iRand(oneIdxs.size()); // choose random epoch
    eps2Skip_.insert(oneIdxs[ridx]);
    oneIdxs.erase(oneIdxs.begin()+ridx);
  }
  cerr << "Skipping " << eps2Skip_.size() << " epochs: "; 
  iterate(eps2Skip_, it) cerr << *it << ", ";
  cerr << endl;
  cerr << "Keeping epochs ";
  for(size_t i=0; i < org_.size(); ++i) {
    if(eps2Skip_.find(i) == eps2Skip_.end()) { 
      cerr << i << ", ";
      lbls_.push_back(org_[i]);
    }
  }
  cerr << endl;
  org_ = lbls_;
}
void GoldLabels::stats() {
  // set noTrainEp_, tstEpochs_, trEpochs_
  cerr << "Have " << noClasses_ << " classes with " << totEps();
  cerr << " total epochs."<< endl;
  noTrainEp_ = totEps() - noTestEp_;
  vector<int> trLbls(noClasses_, 0), tstLbls(noClasses_, 0);
  for(int i=0; i < noTrainEp_; ++i) {
    trEpochs_.insert(i);
    ++trLbls[lbls_[i]];
  }
  cerr << "Number of train epochs: " << noTrainEp_ << "\t";
  for(size_t i=0; i < noClasses_; ++i) {
    cerr << trLbls[i] << " of class " << i << ", ";
  }
  cerr << endl;
  for(int i=noTrainEp_; i < totEps(); ++i) {
    tstEpochs_.insert(i);
    ++tstLbls[lbls_[i]];
  }
  cerr << "Number of test epochs:" << noTestEp_ << "\t";
  for(size_t i=0; i < noClasses_; ++i) {
    cerr << tstLbls[i] << " of class " << i << ", ";
  }
  cerr << endl;
}
void GoldLabels::buildClasses() {
  //org_ = TS::binary_band(org_);
  //return;
  iterate(org_, it) {
    if(*it < -0.2) {
      *it = 0;
    }
    else if((*it >= -0.2) && (*it < -0.1)) {
      *it = 1;
    }
    else if((*it >= -0.1) && (*it < 0)) {
      *it = 2;
    }
    else if((*it >= 0) && (*it < 0.1)) {
      *it = 3;
    }
    else if((*it >= 0.1) && (*it < 0.2)) {
      *it = 4;
    }
    else { // *it >= 0.2 
      *it = 5;
    }
  }
}
