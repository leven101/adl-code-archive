#include "classify.h"

// parameter definitions
const ParamDefs paramdefs[] = {
  //name, default value, abbrev, type, msg
  {"c-spam", ".1", "s", Parameters::kFloatValue, "spam class slack parameter C_s"},
  {"c-nonspam", ".1", "ns", Parameters::kFloatValue, "nonspam class slack parameter C_ns"},
  {"classifiers", "1", "cc", Parameters::kFloatValue, "number of classifiers to instantiate"},
  {"epochs", "1", "", Parameters::kIntValue, "epoch iterations"},
  {"epsilon", ".5", "e", Parameters::kFloatValue, "epsilon for regression"},
  {"gamma", "1", "g", Parameters::kFloatValue, "edge weights gamma"},
  {"gr", Parameters::kFalseValue, "gr", Parameters::kBoolValue, "indicate to GR or not"},
  {"lambda", "1", "l", Parameters::kFloatValue, "graph regularizer weight lambda"},
  {"load-raw-feats", Parameters::kFalseValue, "rf", Parameters::kBoolValue, "load raw features (instead of normalized ones)"},
  {"load-hg-feats", Parameters::kFalseValue, "hgf", Parameters::kBoolValue, "load hostgraph features"},
  {"percent-data", "100", "%", Parameters::kIntValue, "percent of data to use"},
  {"test-num", "200000", "x", Parameters::kIntValue, "test amount to load (debugging)"},
  {"top-dir", "/Users/abby/work/webspam/data/", "dir", Parameters::kStringValue, "top level directory for data"},
  {"year", "2007", "", Parameters::kStringValue, "year of data"},
  {"slack", ".01", "C", Parameters::kFloatValue, "class independent slack parameter"},
  {"start", "", "", Parameters::kIntValue, ""},
  {"end", "", "", Parameters::kIntValue, ""},
};
int main(int argc, char** argv) {
  Parameters params(argc, argv, paramdefs, NumOfParams(paramdefs));
  Classifier c(params, atoi(params.getParam("classifiers").c_str()));
  c.testBagging();
  return 1;
}

