// g++ manipData.cpp -I ~/work/src/include/ ~/work/src/file.o 
#include "file.h"
#include "utils.h"
#include "types.h"
#include "dirList.h"

using namespace std;

set<string> primeSrcs_;
void loadSources() {
  primeSrcs_.insert("AP");  //associated press
  primeSrcs_.insert("DJ");  //dow jones
  primeSrcs_.insert("LB");  //reuters
  primeSrcs_.insert("MA");  //market news international
  primeSrcs_.insert("WS");  //wall street journal
  cerr << "Have " << primeSrcs_.size() << " prime data sources." << endl;
}
void rename(char** argv) {
  DirList dir(argv[1]); 
  int month(0);
  vector<string> vs;
  iterate(dir.cFileList_, itr) {
    ++month;
    Utils::splitToStr(*itr, vs, "/");
    string newDir("/");
     for(int i=0; i < vs.size()-1; ++i) {
      newDir += vs[i] + "/";
     }
     string smth = Utils::IntToStr(month);
     while(smth.size() < 2) {
       smth = "0" + smth;
     }
     newDir += smth; 
     string cmd = "mv " + *itr + " " + newDir;
     //cerr << cmd << endl;
     system(cmd.c_str());
     if(month == 12) month = 0;
  }
}
void splitIntoSources() {
  const string headDir = "/Users/ablev/work/solid/data/text/nonfarm";
  loadSources();
  DirList dir(headDir, "nonfarm-*/*/*/*.txt"); 
  vector<string> vs, varticle;
  iterate(dir.cFileList_, itr) {
    cerr << "Processing directory " << *itr << "..." << endl;
    Utils::splitToStr(*itr, vs, "/");
    string year = vs[vs.size()-3];
    string month = vs[vs.size()-2];
    FileHandler fin(*itr, ios::in); 
    string line;
    while(getline(fin, line)) {
      Utils::trim(line);
      if(line == "AN") {
        varticle.push_back(line);
        string source;
        iterate(varticle, it) { // save to source specific directory
          if(it->find("SC ") == 0) {
            Utils::splitToStr(*it, vs, " ");
            assert(vs.size() == 2);
            Utils::trim(vs[1]);
            source = Utils::uppercase(vs[1].substr(0,2));
            if(source == "J") source = "WS";
            if(source == "AW") source = "WS";
            if(source == "NY") source = "WS";
            if(source == "CM") source = "DJ";
            if(source == "FW") source = "DJ";
            break;
          }
        }
        if(primeSrcs_.size() && (primeSrcs_.find(source) == primeSrcs_.end())) { // if source is not a prime source
          source = "OT";
        }
        string dirPath = headDir + "/sources/" + source + "/" + year + "/" + month + "/";
        string cmd = "mkdir -p " + dirPath;
        system(cmd.c_str());
        string fname = dirPath + "/raw.txt"; 
        FileHandler fout(fname, ios::out|ios::app);
        iterate(varticle, it) {
          fout << *it << endl;
        }
        fout.close();
        varticle.clear();
      }
      else {
        varticle.push_back(line);
      }
    }
  }
}
void filterByKeywords(char** argv) {
  // filter sentences by keywords
  set<string> keywords_;
  keywords_.insert("nonfarm");
  keywords_.insert("forecast");
  DirList dir(argv[1], "*.sntcs"); 
  vector<string> vfn;
  string sline, fline;
  iterate(dir.cFileList_, itr) {
    string sntname = *itr;
    // build new file names
    Utils::splitToStr(sntname, vfn, ".");
    string ftname(""), kwsntcs(""), kwfeats("");
    int noOld(0), noNew(0);
    for(int i=0; i < vfn.size()-1; ++i) {
      ftname += vfn[i] + ".";
      kwsntcs += vfn[i] + ".";
      kwfeats += vfn[i] + ".";
    }
    ftname += "feats";
    kwsntcs += "kwsntcs";
    kwfeats += "kwfeats";
    // open and save data subset
    FileHandler fsnt(sntname, ios::in);
    FileHandler ffeat(ftname, ios::in);
    FileHandler fkwsnt(kwsntcs, ios::out);
    FileHandler fkwfeats(kwfeats, ios::out);
    while(getline(fsnt, sline)) {
      getline(ffeat, fline);
      ++noOld;
      Utils::trim(sline);
      Utils::trim(fline);
      iterate(keywords_, kit) {
        if(sline.find(*kit) != string::npos) {
          fkwsnt << sline << endl;
          fkwfeats << fline << endl;
          ++noNew;
          break;
        }
      }
    }
    cerr << "\t" << noOld << " " << sntname << endl; 
    cerr << "\t" << noNew << " " << kwsntcs << endl; 
    fsnt.close();
    ffeat.close();
    fkwsnt.close();
    fkwfeats.close();
  }
}
int main(int argc, char** argv) {
  if(argc == 1) {
    splitIntoSources();
  }
  else
    //filterByKeywords(argv);
    rename(argv);
  return 1;
}
