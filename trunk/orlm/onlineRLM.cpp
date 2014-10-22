#include "onlineRLM.h"
#include "file.h"
#include "params.h"
#include "dirList.h"
#include "stats.h"

// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"memory", "2", "m", Parameters::kIntValue, "Total megabytes for filter"},
  {"width", "20", "w", Parameters::kIntValue, "Width in bits for each ngram"},
  {"bucket-size", "100", "b", Parameters::kIntValue, "Number of cells in each bucket"},
  {"order", "5", "o", Parameters::kIntValue, "N-gram order"},
  {"qbase", "8", "qb", Parameters::kFloatValue, "Quantization base (used as 2^(1/qbase))"},
  {"top-dir", "_NO_DIR_", "dir", Parameters::kStringValue, "Directory with training data."},
  {"data-ext", "*", "ext", Parameters::kStringValue, "Data file's extensions."},
  {"write", "orlm.gz", "", Parameters::kStringValue, "Path to save trained ORLM to file"},
  {"orlm", "_NO_ORLM_", "", Parameters::kStringValue, "Path to read in trained ORLM from file"},
  {"go", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal training"}, 
  {"query", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal querying"},
  {"adapt", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal adaption"},
  {"stats", Parameters::kFalseValue, "", Parameters::kBoolValue, "Signal stat collection"},
  {"analyze", Parameters::kFalseValue, "", Parameters::kBoolValue, "Analyze a ORLM"},
  {"input", "_NO_FILE_", "d", Parameters::kStringValue, ""},
  {"get-counts", Parameters::kFalseValue, "gc", Parameters::kBoolValue, "Return counts instead of logprobs"},
  {"adapt-severe", Parameters::kFalseValue, "", Parameters::kBoolValue, "Remove all unused events before adapting"},
  {"xx", "300000000", "", Parameters::kIntValue, "Keys to target while testing"}
};

typedef std::map<std::vector<string>, count_t> map_t;
Vocab vocab;

int getPeriod(std::string str) {
  const date_t period = "stream";
  return atoi(str.substr(str.rfind(period) 
    + period.size()).c_str());
}
void adapt(Parameters& params) {
  count_t seen(0), unseen(0);
  time_t start, finish;
  int order = atoi(params.getParamValue("order").c_str());
  // load orlm
  string orlm = params.getParamValue("orlm");
  FileHandler fin(orlm, std::ios::in|std::ios::binary);
  OnlineRLM<count_t> ph(&fin, order);
  fin.close();
  //ph.analyze();
  //ph.countHits();
  //ph.countPrefixes();
  //ph.removeNonMarked(); // for SEVERE deletions
  time(&start);
  //map_t* unseenGrams = new map_t;
  map_t unseenGrams;
  DirList dir_(params.getParamValue("top-dir"),    // which directory
               params.getParamValue("data-ext"));  // which file extension
  iterate(dir_.cFileList_, itr) { // open each file
    // only for this adaption period 
    cerr << "Processing " << *itr << endl;
    // open each file
    FileHandler gramsFile(*itr, std::ios::in);
    string line; 
    std::vector<string> ngram, vtmp;
    while(getline(gramsFile, line) && (seen + unseen < atoi(params.getParam("xx").c_str()))) {
      // split up ngram and associated count
      if(Utils::splitToStr(line, vtmp) != 2) continue;
      Utils::splitToStr(vtmp[0], ngram, " "); // split tokens in ngram
      int value = atoi(vtmp[1].c_str());
      int vals[ngram.size()];
      for(int i = 0; i < (int)ngram.size(); ++i) vals[i] = 0;
      bool binserted(false);
      int fnd = ph.sbsqQuery(ngram, vals, true);
      if(fnd == (int)ngram.size()) { // if sequence found matches this ngram
        ++seen;
        //value = std::max(vals[0], value);// keep largest value
        value = vals[0] + value;
        binserted = ph.update(ngram, value); // update to mark recent
      }
      if(!binserted) { // if ngram not handled above
        unseenGrams[ngram] += value; // store in buffer 
        ++unseen;
      }
      if((seen + unseen) % 1000000 == 0) 
        cerr << "Processed " << seen << " (seen) + " << unseen << " (new) ngrams\n";
      if(unseenGrams.size() >= 10000000) {
        //ph.heurDelete(unseenGrams.size() * 3);  // CONSERVATIVE DELETE 
        // insert all items in unseenGrams
        ph.vocab_->makeOpen();
        cerr << "inserting " << unseenGrams.size() << endl;
        for(map_t::iterator it = unseenGrams.begin();
            it != unseenGrams.end(); ++it) {
          ph.insert(it->first, it->second); // new counts
        }
        ph.vocab_->makeClosed();
        unseenGrams.clear();// clear unseenGrams
      }
    }
    gramsFile.close();
  } // end file iteration
  // insert final items in unseenGrams
  ph.vocab_->makeOpen();
  cerr << "inserting " << unseenGrams.size() << endl;
  //ph.heurDelete(unseenGrams.size() * 3);  // CONSERVATIVE DELETE 
  for(map_t::iterator it = unseenGrams.begin(); it != unseenGrams.end(); ++it) {
    ph.insert(it->first, it->second); // new counts
  }
  ph.vocab_->makeClosed();
  unseenGrams.clear();// clear unseenGrams
  // clear all hit markings
  ph.clearMarkings();
  time(&finish);
  cerr << "time taken was " << difftime(finish, start) << " seconds\n";
  //save orlm
  FileHandler fout(orlm + ".adapted", std::ios::out|std::ios::binary);
  ph.save(&fout);
  fout.close();
}
void load(Parameters& params) {
  cerr << "Querying ORLM\n";
  string fname = params.getParamValue("orlm"); 
  FileHandler fin(fname, std::ios::in|std::ios::binary);
  OnlineRLM<count_t> ph(&fin, 5);
  //ph.analyze();
  fin.close();
  FileHandler ftest(params.getParamValue("input"), std::ios::in);
  string line;
  std::vector<string> ngram, vtmp;
  while(getline(ftest, line)) {
    Utils::trim(line);
    if(Utils::splitToStr(line, ngram, " ") > 0) { // split tokens in ngram
      int len = ngram.size();
      wordID_t wrdIDs[len];
      for(int i=0; i < len; i++) {
        wrdIDs[i] = ph.vocab_->getWordID(ngram[i]);
      }
      if(params.getParamValue("get-counts") == Parameters::kTrueValue) {
        int vals[len];
        int fnd = ph.sbsqQuery(wrdIDs, len, vals, false);
        cout << (fnd == len ? vals[0] : 0) << endl;
      }
      else
        cout <<  ph.getProb(wrdIDs, len, NULL) << endl;
    }
  }
  ftest.close();
  //save LM with markings
  //Utils::rtrim(fname, ".gz");
  //FileHandler fout(fname + ".marked.gz", std::ios::out|std::ios::binary, false);
  //ph.save(&fout);
  //fout.close();
}
void train(Parameters& params) {
  time_t start, finish;
  time(&start);
  cerr << "Initializing Dynamic Bloomier filter...\n";
  srand(time(NULL));
  count_t mem = atoi(params.getParamValue("memory").c_str()),
    target = atoi(params.getParamValue("xx").c_str()),
    width = atoi(params.getParamValue("width").c_str()),
    order = atoi(params.getParamValue("order").c_str()),
    bucketRange = atoi(params.getParamValue("bucket-size").c_str());
  float qBase = atof(params.getParamValue("qbase").c_str());
  OnlineRLM<count_t> ph(mem, width, bucketRange, order, &vocab, qBase);
  if(params.getParamValue("go") != Parameters::kTrueValue) {
    ph.analyze();
    return;
  }
  else cerr << "Storing ngrams...\n";
  count_t keys = 0, modVal = (count_t)std::min(10000000.0, target *.1);
  FileHandler fin(params.getParamValue("input"), std::ios::in);
  string line;
  std::vector<string> vline, vwords;
  while(getline(fin, line) /*&& keys < target*/) {
    /*Utils::splitToStr(line, vline);
    std::vector<string> vwords(vline.size() - 1);
    for(int i=1; i < vline.size(); ++i) {
      vwords[i-1] = vline[i];
      Utils::trim(vwords[i-1]);
    }
    ph.insert(vwords, atoi(vline[0].c_str()));*/
    if(Utils::splitToStr(line, vline) >= 2) { // split count from ngram
      Utils::trim(vline[0]);
      Utils::splitToStr(vline[0], vwords, " "); // split tokens in ngram
      ph.insert(vwords, atoi(vline[1].c_str()));
      if(++keys % modVal == 0) cerr << keys << endl;
    }
  }
  fin.close();
  time(&finish);
  ph.analyze();
  ph.countPrefixes();
  cerr << "corpus size=" << ph.corpusSize() << endl;
  cerr << "total keys=" << keys << endl;
  cerr << "vocab size=" << vocab.size() << endl;
  FileHandler fout(params.getParamValue("write"), std::ios::out|std::ios::binary);
  ph.save(&fout);
  fout.close();
  cerr << "time taken was " << difftime(finish, start) << " seconds\n";
  return;
}
void analyze(Parameters& params) {
  string fname = params.getParamValue("orlm"); 
  count_t order = atoi(params.getParamValue("order").c_str());  
  FileHandler fin(fname, std::ios::in|std::ios::binary);
  OnlineRLM<count_t> ph(&fin, order);
  ph.analyze();
  ph.countHits();
  ph.countPrefixes();
  fin.close();
}
bool insertWithProb(const int stream, int order) {
  // return probability of insert with inverse proportion
  // of stream throughput 
  float prob;
  switch(stream) {
    case 2:
      prob = 1.0; 
      break;
    default:
      prob = order < 5 ? 1.0 : .3;
  }
  return RAND_MAX * prob >= rand();
}
void multiStream(Parameters& params) {
/*implements a max-size adaptation scheme*/
  OnlineRLM<count_t>* lm = 0;
  Stats stats(3,5);
  srand(time(NULL));
  int order = atoi(params.getParamValue("order").c_str());
  bool seed(false);
  count_t mem(0);
  DirList dir_(params.getParamValue("top-dir"), params.getParamValue("data-ext"));
  // load lm if one already exists
  if(params.getParamValue("orlm") != "_NO_ORLM_") { 
    FileHandler fin(params.getParamValue("orlm"), std::ios::in|std::ios::binary);
    lm = new OnlineRLM<count_t>(&fin, order);
    lm->vocab_->makeOpen();
    fin.close();
    lm->removeNonMarked();  // clear old data here
    lm->clearMarkings(); // clear all hit markings
    mem = lm->bucketsMemUse();
  }
  else { // else create empty LM
    seed = true;
    cerr << "Initializing Dynamic Bloomier filter...\n";
    srand(time(NULL));
    mem = atoi(params.getParamValue("memory").c_str());
    count_t  width = atoi(params.getParamValue("width").c_str()),
      order = atoi(params.getParamValue("order").c_str()),
      bucketRange = atoi(params.getParamValue("bucket-size").c_str());
    float qBase = atof(params.getParamValue("qbase").c_str());
    lm = new OnlineRLM<count_t>(mem, width, bucketRange, order, &vocab, qBase);
    if(params.getParamValue("go") != Parameters::kTrueValue) {
      lm->analyze();
      delete lm; 
      return;
    }
  }
  assert(lm);
  // open all streams 
  short streamNo = dir_.cFileList_.size(); 
  FileHandler** fstreams = new FileHandler*[streamNo]; 
  for(int i=0; i < streamNo; ++i) { 
    fstreams[i] = new FileHandler(dir_.cFileList_[i], std::ios::in);
  }
  // now round robin -- get an ngram from each stream until memory constraint is reached
  int curStream(-1), seen(0), unseen(0);
  string line;
  std::vector<string> ngram, vline;
  cerr << "Adding new ngrams...\n";
  while((mem * .01) >= lm->hpDictMemUse()) { // allow some percentage for hpDict
    curStream = (curStream + 1) % streamNo;
    if(getline(*(fstreams[curStream]), line)) {
      //insert or update...
      assert(Utils::splitToStr(line, vline) >= 2);// split count from ngram
      Utils::trim(vline[0]);
      Utils::splitToStr(vline[0], ngram, " "); // split tokens in ngram
      int value = atoi(vline[1].c_str());
      if((curStream == 0) && seed) {
        ++unseen;
        lm->insert(ngram, value);
      }
      else {
        int vals[ngram.size()];
        for(int i = 0; i < (int)ngram.size(); ++i) vals[i] = 0;
        bool binserted(false);
        int fnd = lm->sbsqQuery(ngram, vals, true);
        if(fnd == (int)ngram.size()) { // if sequence found matches this ngram
          ++seen;
          //value = std::max(vals[0], value);// keep largest value
          value = vals[0] + value;
          binserted = lm->update(ngram, value); // update to mark recent
        }
        else {
          ++unseen;
          lm->insert(ngram, value);
        }
      }
      stats.add(curStream+1, ngram.size());
      if((seen + unseen) % 10000000 == 0) 
        cerr << "Processed " << seen << " (seen) + " << unseen << " (new) ngrams\n";
    }
    else if(streamNo == 1) break;
  }
  lm->clearMarkings(); // clear all hit markings
  lm->analyze();
  FileHandler fout(params.getParamValue("write"), std::ios::out|std::ios::binary);
  lm->save(&fout);
  fout.close();
  for(int i=0; i < streamNo; ++i) { // close file streams
    fstreams[i]->close();
    delete fstreams[i];
  }
  delete[] fstreams;
  delete lm;
  stats.print();
}
void getStreamStats(Parameters& params) {
  Stats stats(1,5);
  string fin = params.getParamValue("input");
  //DirList dir_(params.getParamValue("top-dir"), params.getParamValue("data-ext"), NULL, true);
  //iterate(dir_.cFileList_, itr) { 
    //int cur_stream = getPeriod(*itr); 
    //cerr << *itr << "-->" << cur_stream << endl;
    //FileHandler gramsFile(*itr, std::ios::in);
    FileHandler gramsFile(fin, std::ios::in);
    string line; 
    std::vector<string> ngram, vline;
    while(getline(gramsFile, line)) {
      /*assert(Utils::splitToStr(line, vline) >= 2);// split count from ngram
      Utils::trim(vline[0]);
      Utils::splitToStr(vline[0], ngram, " "); // split tokens in ngram
      stats.add(cur_stream, ngram.size());
      */
      Utils::splitToStr(line, vline);
      std::vector<string> vwords(vline.size() - 1);
      for(int i=1; i < vline.size(); ++i) {
        vwords[i-1] = vline[i];
        Utils::trim(vwords[i-1]);
      }
      stats.add(1, vwords.size());
    }
  //}
  stats.print();
}
int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  if(params.getParamValue("query") == Parameters::kTrueValue) load(params);
  else if(params.getParamValue("adapt") == Parameters::kTrueValue) adapt(params);
  else if(params.getParamValue("analyze") == Parameters::kTrueValue) analyze(params);
  else if(params.getParamValue("stats") == Parameters::kTrueValue) getStreamStats(params);
  else train(params);
  /*else multiStream(params);*/
  return 1;
}
