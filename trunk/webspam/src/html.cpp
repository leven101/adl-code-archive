#include "html.h"

map<string, int> hostnames_;
set<int> train_set_, test_set_;
static void loadHostnames();
static bool extractTHostFromWarcData(string);
static void sortRawHtml(string&);
void findMissingSites();
void buildSeedLM(OnlineRLM<unsigned>*);
static void testLM(OnlineRLM<unsigned>*, string, bool = false);

int main(int argc, char** argv) {
  loadHostnames();
  spLM = new OnlineRLM<unsigned>(50, 20, 50, 3, &vcb_, 8); // LM for nonspam   
  nspLM = new OnlineRLM<unsigned>(300, 20, 50, 3, &vcb_, 8); // LM for nonspam   
  FileHandler fseed("../data/2007/lms/seedLM.orlm.gz", ios::in|ios::binary);
  seedLM = new OnlineRLM<unsigned>(&fseed, order_);
  fseed.close();
  DirList dir("../data/2007/raw/sites/", "*.warc.gz");
  int cnt(0), begin = atoi(argv[1]), end = atoi(argv[2]);
  iterate(dir.cFileList_, dit) {
    int id = atoi(dit->substr(dit->find_last_of("/") + 1, dit->find(".warc.gz")).c_str());
    ++cnt;
    if(cnt < begin) continue;
    else if(cnt >= end) break;
    cerr << "processing node " << id << " (" << cnt << ")" << endl;
    processWarcFile(id, features_[id]);
    FileHandler fFeats("./newFeats.csv", ios::out|ios::app);
    fFeats << id << ",";
    iterate(features_[id], vvit) {
      fFeats << *vvit << ","; 
    }
    fFeats << endl;
    fFeats.close();
  } // end dir iterate
  // save extracted features and LMs
  cerr << "Saving LMs and features...\n";
  FileHandler fSpLM("../data/2007/lms/spLM.orlm.gz", ios::out|ios::binary);
  FileHandler fNspLM("../data/2007/lms/nspLM.orlm.gz", ios::out|ios::binary);
  fSpLM.close();
  fNspLM.close();
  // free memory
  cerr << "Freeing memory...\n";
  delete spLM;
  delete nspLM;
  delete seedLM;
  return 1;
}
static bool extractTHostFromWarcData(string swarc) {
//warc/0.9 10757 response http://www.mattenltd.co.uk/ 20060920234350 message/http uuid:c6f7927d-aaea-4e53-b121-c4a594218d8a
  vector<string> tmpV;
  Utils::splitToStr(swarc, tmpV, " ");
  assert(tmpV[3].substr(0,7) == "http://");
  string name = tmpV[3].substr(7);
  name = name.substr(0, name.find("/")); 
  cerr << name << endl; 
  map<string, int>::iterator mit = hostnames_.find(name);
  if(mit == hostnames_.end()) {
    //cerr << "EffedUpIsh. Couldn't find " << name << endl;
    return false;
  }
  int hostID = mit->second; 
  if((train_set_.find(hostID) != train_set_.end()) || 
    (test_set_.find(hostID) != test_set_.end())) {
    return true;
  }
  else 
    return false;
}
static void loadHostnames() {
  vector<string> tmpV;
  string line;
  string infile = "../data/2007/training/WEBSPAM-UK2007-hostnames.txt"; 
  FileHandler fin(infile, ios::in);
  while(getline(fin, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int id = atoi(tmpV[0].c_str());
    Utils::trim(tmpV[1]);
    hostnames_[tmpV[1]] = id;
  }
  fin.close();
  infile = "../data/2007/training/webspam-uk2007-set1-labels.txt"; 
  FileHandler f2(infile, ios::in);
  while(getline(f2, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int id = atoi(tmpV[0].c_str());
    train_set_.insert(id);
  }
  f2.close();
  infile = "../data/2007/testing/webspam-uk2007-set2-labels.txt"; 
  FileHandler f3(infile, ios::in);
  while(getline(f3, line)) {
    Utils::splitToStr(line, tmpV, " ");
    int id = atoi(tmpV[0].c_str());
    test_set_.insert(id);
  }
  f3.close();
  features_.resize(hostnames_.size());
}
static void sortRawHtml(string& infile) {
  cerr << infile << endl;
  FileHandler** fhandlers;
  fhandlers = new FileHandler*[hostnames_.size()];
  for(size_t i=0; i < hostnames_.size(); ++i)
    fhandlers[i] = 0;
  FileHandler fin(infile, ios::in);
  string line;
  int id = -1;
  vector<string> tmpV;
  while(getline(fin, line)) {
    Utils::trim(line);
    if(line.empty()) continue;
    // seperate by website 
    if(COMPARE("warc/0.9", line.substr(0,8).c_str())) { 
      Utils::splitToStr(line, tmpV, " ");
      assert(tmpV[3].substr(0,7) == "http://");
      string name = tmpV[3].substr(7);
      name = name.substr(0, name.find("/")); 
      map<string, int>::iterator mit = hostnames_.find(name);
      assert(mit != hostnames_.end());
      id = mit->second;
      assert((train_set_.find(id) != train_set_.end()) ||
          (test_set_.find(id) != test_set_.end()));
      if(fhandlers[id] == 0) {  // for every website create filehandler
        string sid = "../data/2007/raw/sites/" + Utils::IntToStr(id) + ".warc";
        fhandlers[id] = new FileHandler(sid, ios::out|ios::app, false);
      }
    }
    *fhandlers[id] << line << endl; 
  }
  // free it all
  cerr << "freeing it all\n";
  for(size_t i=0; i < hostnames_.size(); ++i) {
    if(fhandlers[i]) {
      fhandlers[i]->close();
      delete fhandlers[i];
    }
  }
  delete[] fhandlers;
}
void findMissingSites() {
  DirList dir("../data/2007/raw/sites/", "*.warc.gz");
  set<int> here;
  iterate(dir.cFileList_, ditr) {
    int id = atoi(ditr->substr(ditr->find_last_of("/") + 1, ditr->find(".warc.gz")).c_str());
    here.insert(id);
  }
  iterate(train_set_, ts)
    if(here.find(*ts) == here.end())
      cerr << *ts << endl; 
  iterate(test_set_, ts)
    if(here.find(*ts) == here.end())
      cerr << *ts << endl; 
}
void buildSeedLM(OnlineRLM<unsigned>* seedLM) {
  FileHandler fin("../data/2007/raw/12-2007-text.pages", ios::in);
  string line;
  vector<string> vs;
  map<string, int> allgrams;
  int cnt(0);
  while(getline(fin, line)) {
    Utils::trim(line);
    if(line == "==P=>>>>=i===<<<<=T===>=A===<=!Junghoo!==>")
      continue;
    if(!line.empty()) {
      if(++cnt % 5000000 == 0) { 
        cerr << "Processed " << cnt << " lines...\n";
        cerr << "adding " << allgrams.size() << " unique grams to seed LM\n";
        iterate(allgrams, git) {
            Utils::splitToStr(git->first, vs, " ");
            seedLM->update(vs, git->second);
        }
        allgrams.clear();
      }
      map<string, int> grams = countStreamNgrams(line, order_);
      if(grams.size() == 0) continue;
      iterate(grams, git) {
        allgrams[git->first] += git->second;
      }
    }
  }
  fin.close();
  FileHandler fout("../data/2007/lms/seedLM.orlm.gz", ios::out|ios::binary);
  seedLM->save(&fout);
  seedLM->analyze();
  fout.close();
}
static void testLM(OnlineRLM<unsigned>* orlm, string fname, bool gc) {
  orlm->analyze();
  orlm->corpusSize();
  FileHandler ftest(fname, std::ios::in);
  string line;
  std::vector<string> ngram;
  while(getline(ftest, line)) {
    Utils::trim(line);
    if(Utils::splitToStr(line, ngram, " ") > 0) { // split tokens in ngram
      int len = ngram.size();
      wordID_t wrdIDs[len];
      for(int i=0; i < len; i++) {
        wrdIDs[i] = orlm->vocab_->getWordID(ngram[i]);
      }
      if(gc == true) {
        int vals[len];
        int fnd = orlm->sbsqQuery(wrdIDs, len, vals, false);
        cout << (fnd == len ? vals[0] : 0) << endl;
      }
      else
        cout <<  orlm->getProb(wrdIDs, len, NULL) << endl;
    }
  }
  ftest.close();
}
