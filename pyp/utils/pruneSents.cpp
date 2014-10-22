//g++ pruneSents.cpp -I ~/src/include/ ~/src/file.o
#include "file.h"
#include "utils.h"
#include "types.h"

using namespace std;

int main(int argc, char** argv) {
  FileHandler fin1(argv[1], ios::in);
  FileHandler fin2(argv[2], ios::in);
  const int max = atoi(argv[3]);
  bool goldAlgs = argc == 5; 
  FileHandler fout1("./chinese.txt", ios::out, false);
  FileHandler fout2("./english.txt", ios::out, false);
  FileHandler* fin3(0), *fout3(0);
  if(goldAlgs) {
    fin3 = new FileHandler(argv[4], ios::in);
    fout3 = new FileHandler("./algs.txt", ios::out, false);
  }
  string line1, line2, line3;
  vector<string> v;
  int src_above_max(0), trg_above_max(0);
  while(getline(fin1, line1)) {
    getline(fin2, line2);
    if(goldAlgs && (!fin3->eof())) 
      getline(*fin3, line3);
    Utils::splitToStr(line1, v, " ");
    if(v.size() >= max) {
      ++src_above_max;
      continue;
    }
    Utils::splitToStr(line2, v, " ");
    if(v.size() >= max) {
      ++trg_above_max;
      continue;
    }
    fout1 << line1 << endl;
    fout2 << line2 << endl;
    if(goldAlgs&& (!fin3->eof())) 
      *fout3 << line3 << endl;
  }
  cerr << "Src sentences above max: " << src_above_max << endl;
  cerr << "Trg sentences above max: " << trg_above_max << endl;
  fin1.close();
  fin2.close();
  if(goldAlgs) {
    fin3->close();
    delete fin3;
  }
  fout1.close();
  fout2.close();
  if(goldAlgs) {
    fout3->close();
    delete fout3;
  }
  return 1;
}

