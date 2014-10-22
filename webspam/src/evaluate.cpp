#include "evaluate.h"

int main(int argc, char** argv) {
  if(argc < 5) {
    cerr << "Not enuf arguements\n";
    exit(-1);
  }
  Evaluate::Fmeasure(atof(argv[1]), atof(argv[2]), 
    atof(argv[3]), atof(argv[4]));
  return 1;
}
