#include "features.h"
#include "html.h"

void Features::normalizeContentFeatures() {
  cerr << "Normalizing all content features... " << endl;
  const int numFeatures(contentFeatures_[0].size());
  string outfile = params_->getParamValue("top-dir") + params_->getParamValue("year");
  outfile += "/features/normalized.content";
  normalize(numFeatures, contentFeatures_, outfile);
}
void Features::normalizeLinkFeatures() {
  cerr << "Normalizing all link-based features... " << endl;
  const size_t numFeatures(linkFeatures_[0].size());
  string outfile = params_->getParamValue("top-dir") + params_->getParamValue("year");
  outfile += "/features/normalized.links";
  normalize(numFeatures, linkFeatures_, outfile);
}
void Features::normalize(int numFeatures, vector<vector<float> >& features,
  string outfile) {
  for(int i=0; i < numFeatures; ++i) { // for each feature
    // replace all values by number of websites less than value
    std::map<float, int> featCounts;
    // count
    for(size_t j=0; j < features.size(); ++j) {
      assert((int)features[j].size() == numFeatures);
      featCounts[features[j][i]]++;
    }
    // count counts
    std::map<float, int> lessThan;
    int numLessThan(0);
    iterate(featCounts, fitr) {
      lessThan[fitr->first] = numLessThan;
      numLessThan += fitr->second;
    }
    // replace 
    for(size_t j=0; j < features.size(); ++j) {
      float tmp = features[j][i];
      features[j][i] = lessThan[tmp];  // number of instances having a value of less than x_ij for feature j
      features[j][i] /= (float)features.size();
    }
  }
  // save normalized values to file 
  if(!outfile.empty()) {
    FileHandler fnormal(outfile, std::ios::out, false);
    for(size_t i=0; i < features.size(); ++i) {
      iterate(features[i], itr)
        fnormal << *itr << ",";
      fnormal << endl;
    }
  }
}
void Features::loadNormalizedContentFeatures() {
  cerr << "Loading normalized content-based features...\n";
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year + "/features/normalized.content";
  if(!Utils::fileExists(infile)) {
    loadRawContentFeatures();
    normalizeContentFeatures();
  }
  else {
    string line;
    int TEST_NUM = atoi(params_->getParamValue("test-num").c_str()), cnt(0);
    FileHandler fin1(infile, std::ios::in);
    while(getline(fin1, line)) {
      vector<float> tmpV;
      if(line[0] == '#') continue;
      Utils::splitToFloat(line, tmpV);
      contentFeatures_.push_back(tmpV);
      if(++cnt >= TEST_NUM) break;
    }
    fin1.close();
  }
}
void Features::loadNormalizedLinkFeatures() {
  cerr << "Loading normalized link features...\n";
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year + "/features/normalized.links";
  if(!Utils::fileExists(infile)) {
    loadRawLinkFeatures();
    normalizeLinkFeatures();
  }
  else {
    string line;
    int TEST_NUM = atoi(params_->getParamValue("test-num").c_str()), cnt(0);
    FileHandler fin1(infile, std::ios::in);
    while(getline(fin1, line)) {
      vector<float> tmpV;
      if(line[0] == '#') continue;
      Utils::splitToFloat(line, tmpV);
      linkFeatures_.push_back(tmpV);
      if(++cnt >= TEST_NUM) break;
    }
    fin1.close();
  }
}
void Features::loadRawLinkFeatures() {
  cerr << "Loading raw link-based features...\n";
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year + 
    "/features/uk-" + year + "-05.link_based_features.csv"; //_transformed.csv";
  string line;
  int TEST_NUM = atoi(params_->getParamValue("test-num").c_str());
  FileHandler flinks(infile, std::ios::in);
  vector<float> tmpVF;
  int website(0);
  while(getline(flinks, line)) {
    if(line[0] == '#') continue;
    Utils::splitToFloat(line, tmpVF);
    vector<float> features(tmpVF.begin()+3, tmpVF.end());
    linkFeatures_.push_back(features);
    if(++website >= TEST_NUM) break;
  }
  flinks.close();
}
void Features::loadRawContentFeatures() {
  cerr << "Loading raw content-based features...\n";
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year + 
    "/features/uk-" + year + "-05.obvious_features.csv";
  string line;
  int TEST_NUM = atoi(params_->getParamValue("test-num").c_str());
  FileHandler fin1(infile, std::ios::in);
  vector<float> tmpVF;
  int website(0);
  while(getline(fin1, line)) {
    if(line[0] == '#') continue;
    Utils::splitToFloat(line, tmpVF);
    vector<float> features(tmpVF.begin()+2, tmpVF.end());
    contentFeatures_.push_back(features);
    if(++website >= TEST_NUM) break;
  }
  fin1.close();
  website=0;
  infile = params_->getParamValue("top-dir") + year +
    "/features/uk-" + year + "-05.content_based_features.csv";
  FileHandler fin2(infile, std::ios::in);
  while(getline(fin2, line)) {
    if(line[0] == '#') continue;
    Utils::splitToFloat(line, tmpVF);
    int targetSite = (int)tmpVF[0];
    assert(tmpVF.size() == 98);
    while(website < targetSite) { // fill missing sites with zeroes for all features
      //cout << website << " versus " << targetSite << endl;
      for(int i=0; i < 96; ++i)
        contentFeatures_[website].push_back(0);
      hostsMissingFeats_.insert(website);
      if(++website >= TEST_NUM) break;
    }
    for(size_t i=2; i < tmpVF.size(); ++i)
      contentFeatures_[website].push_back(tmpVF.at(i));
    if(++website >= TEST_NUM) break;
  }
  fin2.close();
  while(website < (int)contentFeatures_.size()) { // handling missing end cases
    for(int i=0; i < 96; ++i)
      contentFeatures_[website].push_back(0);
    hostsMissingFeats_.insert(website);
    ++website;
  }
  // save missing hosts to file 
  infile = params_->getParamValue("top-dir") + year + 
    "/training/hostsMissingFeats.txt";
  FileHandler fout(infile, std::ios::out, false);
  iterate(hostsMissingFeats_, sitr) fout << *sitr << endl;
  fout.close();
}
void Features::loadHostsMissingFeats() {
  cerr << "Loading hosts missing content-based features...\n";
  string infile = params_->getParamValue("top-dir") + 
    params_->getParamValue("year") + "/training/hostsMissingFeats.txt";
  string line;
  FileHandler fin(infile, std::ios::in);
  while(getline(fin, line)) {
    Utils::trim(line);
    hostsMissingFeats_.insert(atoi(line.c_str()));
  }
  fin.close();
  assert(hostsMissingFeats_.size() == 5168); // verify all is okay
}
void Features::printStuff() {
  cout << "size of contentFeatures = " << contentFeatures_.size() << endl;
  iterate(contentFeatures_, vitr) {
    piterate(vitr, itr)
      cout << *itr << "  ";
    cout << endl;
  }
}
//TODO: test this properly
void Features::streamLMFeatures(int id) {
  static const vector<float> dummy(15,0);
  vector<float> feats;
  processWarcFile(id, feats);
  if(feats.size() == 0) 
    lmFeatures_[id] = dummy;
  else {
    for(size_t i=0; i < feats.size(); ++i) {
      if(feats[i] == -0) feats[i] = 0; 
      lmFeatures_[id].push_back(logisticF(feats[i])); // normalize
    }
  }
}
void Features::loadLMFeatures() {
  string infile = "../data/2007/features/newFeats.csv";
  FileHandler fin(infile, std::ios::in);
  string line;
  int TEST_NUM = atoi(params_->getParamValue("test-num").c_str()), cnt(0);
  while(getline(fin, line)) {
    vector<float> tmpV;
    if(line[0] == '#') continue;
    Utils::splitToFloat(line, tmpV);
    int id = tmpV[0]; 
    smallLMSet.insert(id);
    for(size_t i=1; i < tmpV.size(); ++i) {
      if(tmpV[i] == -0) tmpV[i] = 0; 
      lmFeatures_[id].push_back(logisticF(tmpV[i]));
    }
    if(++cnt >= TEST_NUM) break;
  }
  vector<float> dummy(15,0);
  for(size_t i=0; i < lmFeatures_.size(); ++i) {
    if(lmFeatures_[i].size() == 0) 
      lmFeatures_[i] = dummy;
  }
  fin.close();
}
