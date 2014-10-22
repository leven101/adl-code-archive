#include "hostGraph.h"

const map<int, int>* HostGraph::getOutEdges(const int nodej) {
  map<int, map<int, int> >::const_iterator graph_itr;
  graph_itr = out_links_.find(nodej);
  if(graph_itr != out_links_.end()) // if node_j is in hostgraph
    return &graph_itr->second;
  return 0;
}
const map<int, int>* HostGraph::getInEdges(const int nodej) {
  map<int, map<int, int> >::const_iterator graph_itr;
  graph_itr = in_links_.find(nodej);
  if(graph_itr != in_links_.end()) // if node_j is in hostgraph
    return &graph_itr->second;
  return 0;
}
int HostGraph::findEdges(const int nodej, const int nodei) {
  int found(0);
  map<int, map<int, int> >::const_iterator graph_itr;
  graph_itr = out_links_.find(nodej);
  if(graph_itr != out_links_.end()) {  // if node_j is in hostgraph
    const map<int, int>& edges = graph_itr->second;
    map<int, int>::const_iterator edge_itr = edges.find(nodei);
    if(edge_itr != edges.end()) // if node_i connects to node_j
      found = edge_itr->second;
  }
  return found;
}
void HostGraph::loadHostGraph() {
  cerr << "Loading raw host graph...";
  string infile = params_->getParamValue("top-dir") + 
    params_->getParamValue("year") + "/features/hostgraph.txt";
  string line;
  vector<string> svec;
  vector<int> ivec;
  int TEST_NUM = atoi(params_->getParamValue("test-num").c_str());
  int website(-1);
  FileHandler fin(infile, std::ios::in);
  while(getline(fin, line)) {
    Utils::trim(line);
    if(!(line.empty() || (website == -1))) {
      Utils::splitToStr(line, svec, " ");
      iterate(svec, itr) {
        Utils::splitToInt(*itr, ivec, ":");
        assert(ivec.size() == 2);
        out_links_[website][ivec[0]] = ivec[1];
        in_links_[ivec[0]][website] = ivec[1];
        //cerr << website << "-->" << ivec[0] << "::" << ivec[1] << endl;
      }
    }
    if(++website >= TEST_NUM) break;
  }
  fin.close();
  cerr << "Loaded " << out_links_.size() << " graph nodes.\n";
  /*iterate(train_labels_, itr) {
    if(out_links_.find(itr->first) == out_links_.end() &&
      (itr->second == SPAM))
        cerr << "not missing " << itr->first << endl;
  }*/
}
const vector<vector<float> >* HostGraph::loadHostFeatures() {
  string year = params_->getParamValue("year");
  string infile = params_->getParamValue("top-dir") + year+ "/features/graphFeatures.";
  infile += params_->getBoolValue("load-raw-feats") ? "raw" : "normalized";
  if(!Utils::fileExists(infile)) {
    hostFeatures_.resize(f_->linkFeatures_.size());
    const int numFeats(10);
    const int featIdx[numFeats] = {16, 17, 38, 39, 30, 31, 6, 7, 14, 15}; 
    int start = atoi(params_->getParam("start").c_str());
    int end = atoi(params_->getParam("end").c_str());
    assert(end <= (int)f_->linkFeatures_.size());
    for(int n=start; n < end; ++n) { // for each website
      if(n % 5000 == 0) cerr << "Processed " << n << " nodes...\n";
      EDGE_PT outLinks = getOutEdges(n);
      EDGE_PT inLinks = getInEdges(n);
      for(int i=0; i < numFeats; ++i) { // get each of the feature values
        //cerr << "webiste " << n << "  feature " << featIdx[i] << " has value " << measure(n, featIdx[i]) << endl;
        hostFeatures_[n].push_back(measure(n, featIdx[i]));
        if(inLinks) {
          hostFeatures_[n].push_back(gengF_2(featIdx[i], inLinks));
          hostFeatures_[n].push_back(gengF_4And5(featIdx[i], inLinks));
          hostFeatures_[n].push_back(gengF_6And7(featIdx[i], inLinks));
        }
        else {
          for(int z=0; z < 3; ++z)
            hostFeatures_[n].push_back(0);
        }
        if(outLinks) {
          hostFeatures_[n].push_back(gengF_3(featIdx[i], outLinks));
          hostFeatures_[n].push_back(gengF_4And5(featIdx[i], outLinks));
          hostFeatures_[n].push_back(gengF_6And7(featIdx[i], outLinks));
        }
        else {
          for(int z=0; z < 3; ++z)
            hostFeatures_[n].push_back(0);
        }
      }
    }
    infile += params_->getParam("start") + "-" + params_->getParam("end");
    FileHandler fout(infile, std::ios::out, false);
    //for(size_t i=0; i < hostFeatures_.size(); ++i) {
    for(int i=start; i < end; ++i) {
      iterate(hostFeatures_[i], hf)
        fout << *hf << ",";
      fout << endl;
    }
    fout.close();
  }
  else {
    cerr << "Loading host graph features from " << infile << endl;
    string line;
    int TEST_NUM = atoi(params_->getParamValue("test-num").c_str());
    FileHandler fin1(infile, std::ios::in);
    vector<float> tmpVF;
    int website(0);
    while(getline(fin1, line)) {
      if(line[0] == '#') continue;
      Utils::splitToFloat(line, tmpVF);
      hostFeatures_.push_back(tmpVF);
      if(++website >= TEST_NUM) break;
    }
    fin1.close();
  }
  cerr << "hostFeatures.size() = " << hostFeatures_.size() << endl;
  //normalizeFeatures();
  //exit(1);
  return &hostFeatures_;
}
