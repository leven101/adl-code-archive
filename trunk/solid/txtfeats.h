#ifndef txtfeats_solid_h
#define txtfeats_solid_h

#include <unistd.h>
#include "header.h"
#include "curl_theysay.h"

class TFIDF;
class DocFeats {
public:
  friend class TFIDF;
  typedef vector<word_t> snt_t;
  DocFeats(const string, bool, Parameters&);
  ~DocFeats() { 
  }
  vector<vf_t> txtFeats_;
  map<wordID_t, int> wrdCnts_;
  float noSnts() { return txtFeats_.size(); }
  float noPos() { return totPos_; }
  float noNeg() { return totNeg_; }
  float pcntPos() { return noSnts() ? totPos_/noSnts() : 0; }
  float pcntNeg() { return noSnts() ? totNeg_/noSnts() : 0; }
  int noNeut() { return totNeut_; }
  string source() { return fname_; }
  void saveFeatures(const string& fname) {
    cerr << "saving features to file " << fname << endl;
    FileHandler fout(fname, std::ios::out, false);
    fout.precision(3);
    for(size_t i=0; i < txtFeats_.size(); ++i) {
      for(size_t j=0; j < txtFeats_[i].size(); ++j) {
        fout << txtFeats_[i][j] << " ";
      }
      fout << endl;
    }
    fout.close();
  }
private:
  int totPos_, totNeg_, totNeut_, totFeats_;
  int curSnt_, processEvery_;
  FileHandler* fsnts_, *ffeats_;
  const string fname_;
  time_t rlStart_, rlEnd_;
  Parameters* params_;
  void theySayFeatsBatch(const snt_t&);
  bool parseSentiment(const string&);
  void loadFeatures(const string);
  void computeFeatures(const string);
  void checkRateLimit(); 
  void sentimentStats();
};
DocFeats::DocFeats(const string fname, bool precomputed, 
  Parameters& params): fname_(fname), params_(&params) {
  if(precomputed) {
    loadFeatures(fname);
    sentimentStats();
  }
  else {
    curSnt_=0;
    computeFeatures(fname);
  }
}
void DocFeats::sentimentStats() {
  // how many negative, positive, neutral sentences 
  totPos_ = 0;
  totNeg_ = 0;
  totNeut_ = 0;
  for(size_t i=0; i < txtFeats_.size(); ++i) {
    if(txtFeats_[i][0] == -1)
      ++totNeg_;
    else if(txtFeats_[i][0] == 1)
      ++totPos_;
    else ++totNeut_;
  }
}
void DocFeats::computeFeatures(const string fname) {
  FileHandler fin(fname, std::ios::in);
  string line;
  vector<string> vtmp;
  fsnts_ = new FileHandler(fname + ".sntcs", std::ios::out, false);
  ffeats_ = new FileHandler(fname + ".feats", std::ios::out, false);
  ffeats_->precision(3);
  //processEvery_ = 1;
  int totSnts(0), sntsPruned(0);
  time(&rlStart_);
  while(getline(fin, line)) {
    Utils::trim(line);
    if(line.empty()) continue;
    if(++totSnts==1) continue;
    Utils::trim(line);
    if(line.find("AN Document") != string::npos) {
      //cout << line << endl;
      *fsnts_ << endl;
      *ffeats_ << endl;
    }
    Utils::splitToStr(line, vtmp, " ");
    if(StrUtils::pruneSnts(line, vtmp)) {
      ++sntsPruned;
      continue;
    }
    ++curSnt_;
    snt_t snt;
    iterate(vtmp, it) {
      if(!StrUtils::pruneWord(*it)) {
        snt.push_back(*it);
      }
    }
    theySayFeatsBatch(snt);
  }
  cerr << "File contined a total of " << totSnts << " sentences." << endl;
  cerr << "Removed " << sntsPruned << " unsavory sentences." << endl;
  cerr << "Processed a total of " << curSnt_ << " good sentences." << endl;
  fin.close();
  fsnts_->close();
  ffeats_->close();
}
void DocFeats::loadFeatures(const string path) {
  // 1 file is 1 epoch's features
  //cerr << "Loading features from file: " << path << endl;
  FileHandler fin(path, std::ios::in);
  string line;
  vector<float> vtmp;
  while(getline(fin, line)) {
    if(line.empty()) continue;
    Utils::splitToFloat(line, vtmp, " ");
    //if(vtmp[0] == 0) continue; // skip neutral sentiment
    vf_t feats;
    for(size_t i=0; i < vtmp.size(); ++i) {
      feats.push_back(vtmp[i]);
    }
    txtFeats_.push_back(feats);
  }
  fin.close();
}
void DocFeats::checkRateLimit() {
  if(curSnt_ % 300 == 0) {
    cerr << curSnt_ << "...";
    time(&rlEnd_);
    float secs = difftime(rlEnd_, rlStart_);
    cerr << "Took " << secs << " to process 300 sentences..." << endl;
    if(secs < 60) {
      cerr << "Taking a TheySayAPI nap..." << endl;
      sleep(60);
    }
    time(&rlStart_);
  }
}
void DocFeats::theySayFeatsBatch(const snt_t& snt) {
  // get sentiment features for this sentence 
  // static string batch; // for batch mode
  string batch;
  iterate(snt, wrdp) { // rebuild sentences
    batch += *wrdp + " "; 
  }
  checkRateLimit();
  //if(batch.size() > processEvery_) {
    assert(batch.size() < 20000);
    string data = theySayAPI(batch);  // only sentiment right now
    if(data.find("{\"errors\"") == 0) {
      cerr << data << endl;
      cerr << "SENT: " << batch << endl;
      if(data.find("exceeded your rate limit") != string::npos) {
        cerr << "Taking a TheySayAPI nap..." << endl;
        sleep(60);
        time(&rlStart_); //restart rateLimit timer
      }
      else if(data.find("used up your quota") != string::npos) {
        assert(false);
      }
    }
    else if(data.size()) {
      if(parseSentiment(data))
        *fsnts_ << batch << endl;
      // batch.clear(); // for batch mode
    }
  //}
}
bool DocFeats::parseSentiment(const string& data) {
  vector<string> vBatch, vSntc, tmp;
  Utils::splitToStrMD(data, vBatch, "}, {");
  iterate(vBatch, it) {
    vf_t feats;
    Utils::splitToStr(*it, vSntc, ",");
    if(vSntc[0].find("POSITIVE") != string::npos) {
      feats.push_back(1); 
    }
    else if(vSntc[0].find("NEGATIVE") != string::npos) {
      feats.push_back(-1); 
    }
    else {
      feats.push_back(0); // neutral  
    }
    Utils::splitToStr(vSntc[1], tmp, ":");
    feats.push_back(atof(tmp[1].c_str()));  // prob(positive)
    Utils::splitToStr(vSntc[2], tmp, ":");
    feats.push_back(atof(tmp[1].c_str()));  // prob(negative)
    Utils::splitToStr(vSntc[3], tmp, ":");
    feats.push_back(atof(tmp[1].c_str()));  // prob(neutral)
    Utils::splitToStr(vSntc[4], tmp, ":");
    feats.push_back(atof(tmp[1].c_str()));  // word count or confidence score 
    //if(feats[1] == 0 && feats[2] == 0 && feats[3] == 0) 
      //return false;
    //txtFeats_.push_back(feats);
    iterate(feats, it) { // save to file
      *ffeats_ << *it << " ";
    }
    *ffeats_ << endl;
    // save sentences to file 
    //Utils::splitToStrMD(*it, vSntc, "\"text\":");
    //Utils::trim(vSntc.back(), " }]\n\"");
    //*fsnts_ << vSntc.back() << endl;
  }
  return true;
}
class TFIDF {
// TF-IDF function from http://lucene.apache.org/core/3_6_1/api/all/org/apache/lucene/search/Similarity.html
public: 
  static vector<vf_t> calculate(const vector<DocFeats*>& docs) {
    map<wordID_t, float> docFreq;
    // calculate number of documents containing terms
    for(size_t d=0; d < docs.size(); ++d) {
      iterate(docs[d]->wrdCnts_, wit) {
        docFreq[wit->first]++;
      }
    }
    // for each document calculate its TF-IDF
    vector<vf_t> tfidfFeats(docs.size());
    for(size_t d=0; d < docs.size(); ++d) {
      vf_t snt_tfidf(docs[d]->wrdCnts_.size());
      iterate(docs[d]->wrdCnts_, wit) {
        wordID_t id = wit->first;
        float df = docFreq.find(id)->second; 
        //cout << "df: " << df << endl;
        float idf = 1.0 + (log((float)docs.size() / df+1));
        //cout << "idf: " << idf << endl;
        float tfidf = pow(wit->second, 0.5) * idf;
        //cout << "tidf: " << tfidf << endl;
        snt_tfidf[id-1] = tfidf;
      }
      tfidfFeats[d] = snt_tfidf;
    }
    return tfidfFeats;
  }
};
#endif
