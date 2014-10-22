#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "file.h"
#include "utils.h"

using namespace std;

int main(int argc, char** argv) {
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
  //for( size_t i=0; i < N_; i++ )
  //  cout << ngrams_[i] << endl;
  string period_ = "", year_;
  for( size_t i = 1; i < 14; i++ ) {
    period_ = (i < 10) ? "0" : "";
    period_ +=  Utils::IntToStr(i);
    if( i <= 8 ) {
      year_ = "1997-";
    } 
    else if( i > 8 ) {
      year_ = "1996-";
      if( i == 13 ) period_ = "08";
    }
    string fullcmd_, cmd_ = "gunzip -fc " + 
      dataPath_ + "count-ngram/" + 
      year_ + period_ + "-*." + gram_ + 
      ".sorted.gz | grep \"\t";
    // bug bug bug bug hack hack hack
    if(gram_ != "unigrams" ) cmd_ += " ";
    string outFilePath = year_ + period_ + "."
                         + gram_ + ".distr.gz";
    FileHandler fout(outFilePath, ios::out, false);
    for( size_t j = 0; j < N_; j++ ) {
      fullcmd_ = cmd_ + ngrams_[j] + "\"$";
      //cout << fullcmd_ << endl;
      //return 0;
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
  } // months
}// main
