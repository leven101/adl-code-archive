#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "file.h"
#include "utils.h"
#include "dirList.h"
#include "types.h"

using namespace std;

int main(int argc, char** argv) {
  string dataPath_ = argv[1];
  string weeks[60] = {"0"};
  {
    FileHandler weeksfile("weeks.data", ios::in);
    int i = 0;
    while( getline(weeksfile, weeks[i]) ) i++;    
  }
  DirList dirlist(dataPath_, "*.sents.gz");
  int cnt=0;
  string thisweek = weeks[cnt], nextweek = weeks[cnt+1],
         weeksfiles = "";
  iterate(dirlist.cFileList_, itr) {
    string tmp = itr->substr(itr->find("/199"));
    Utils::trim(tmp, "/.sents.gz");
    //cout << "tmp=" << tmp << endl;
   
    if( tmp >= thisweek && tmp < nextweek ){
      weeksfiles += " " + *itr;
    }
    else {
      string cmd = "zcat " + weeksfiles + 
                   "| gzip -f > week" +  
                   Utils::IntToStr(cnt+1) + 
                   ".sents.gz";
      cout << cmd << endl;
      //int r = system(cmd.c_str());
      weeksfiles = " " + *itr;
      thisweek = weeks[++cnt]; 
      nextweek = weeks[cnt+1];
    }
  }
}// main

/*
  if( argc < 4 ) {
    cerr << "Needs args:\n(1)\tNumber of ngrams\n";
    cerr << "(2)\tcount file paths\n(3)\t{1,2,3)-grams\n";
    cerr << "(4)\tstarting line\n";
    exit(1);
  }
  size_t N_ = atoi(argv[1]), 
    startline_ = argc == 5 ? argc : 0;
  string ngrams_[N_], line, gram_ = argv[3];
  string dataPath_ = argv[2];
  {
    string mainFilePath = dataPath_ + "reuters." 
                          + gram_ + ".sorted.gz";
    vector<string> tmp;
    FileHandler mainFile(mainFilePath, ios::in); 
    for( size_t i=0; i< startline_; i++ )
      getline(mainFile, line);
    for( size_t i = 0; i < N_; i++ ) {
      getline(mainFile, line);
      Utils::tokenizeToStr(line, tmp, "\t");
      Utils::trim(tmp[0]);
      ngrams_[i] = tmp[0];
      tmp.clear();
    }
  }
  {
    FileHandler weeksfile("weeks.data", ios::in);
    int i = 0;
    while( getline(weeksfile, weeks[i]) ) i++;    
  }
  string filepath_ = dataPath_ + "count-ngram/";
  DirList dirlist(filepath_, "*.unigrams.sorted.gz");
  int cnt=0;
  string thisweek = weeks[cnt], nextweek = weeks[cnt+1],
         weeksfiles = "";
  iterate(dirlist.cFileList_, itr) {
    string tmp = itr->substr(itr->find("/199"));
    Utils::trim(tmp, "/.unigrams.sorted.gz");
    //cout << "tmp=" << tmp << endl;
    
    if( tmp >= thisweek && tmp < nextweek ){
      weeksfiles += " " + *itr;
    }
    else {
      //do all the stuff
      string fullcmd_, cmd_ = "gunzip -fc " + 
        weeksfiles + " | grep \"\t";
      // bug bug bug bug hack hack hack
      if(gram_ != "unigrams" ) cmd_ += " ";
      string outFilePath = "../output/week" +
                           Utils::IntToStr(cnt+1) + 
                           "." + gram_ + ".distr.gz";
      FileHandler fout(outFilePath, ios::out, false);
      for( size_t j = 0; j < N_; j++ ) {
        fullcmd_ = cmd_ + ngrams_[j] + "\"$";
        //cout << fullcmd_ << endl;
        FILE *fp = popen(fullcmd_.c_str(), "r");
        if( !fp ) { 
          cerr << "Error in \"" + fullcmd_ + "\"\n";
          return 0;
        }
        int gramCount = 0;
        char psBuffer[1024];
        while( !feof(fp) ) {
          if( fgets(psBuffer, sizeof(psBuffer), fp ) != NULL ) {
            if( psBuffer[strlen(psBuffer)-1] == '\n' )
              psBuffer[strlen(psBuffer)-1] = 0;
            gramCount += atoi(psBuffer); 
          }
        }
        fout << gramCount << "\t" << ngrams_[j] << endl;
        pclose(fp);
      } //ngrams
      fout.flush();
      weeksfiles = " " + *itr;
      thisweek = weeks[++cnt]; 
      nextweek = weeks[cnt+1];
    }
  }
 */ 
