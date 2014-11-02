#include "ibcc.h"

using namespace std;
// paramter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"input", "./wcOutput", "i", Parameters::kStringValue, "Input file with weak classifier distributions"},
  {"test-epochs", "1", "t", Parameters::kIntValue, "Number of instances to test"},
  {"skip", "0", "", Parameters::kIntValue, "First number of lines in input file to skip"}, 
  {"vb-itr", "5", "vb", Parameters::kIntValue, "Number of Variational Bayes iterations to run"},
  {"unsupervised", Parameters::kFalseValue, "unsup", Parameters::kBoolValue, "Run IBCC in an unsupervised manner"},
  {"nu0", "", "", Parameters::kFloatValue, "Prior values for nu (use as -nu0 \"class0_val class1_val ...\" for each class"},
  {"alp0", "", "", Parameters::kFloatValue, "Prior values for alpha (use as -alp0 \"[matrix of space delimited vals]\""},
};
int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  cerr << endl << "Starting IBCC..." << endl;
  //CIBCC ibcc(params, weakClfs, gldLbls); // send data straight from models
  CIBCC ibcc(params);
  ibcc.train();
  ibcc.backtest();
  return 1;
}
