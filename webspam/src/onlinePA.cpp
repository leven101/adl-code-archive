#include "onlinePA.h"

bool operator+=(vector<float>& x, vector<float>& y) {
  int size = std::min(x.size(), y.size());
  for(int i=0; i < size; ++i) 
    x[i] += y[i];
  return true; 
}
// for batch testing
map<int, float> V, scores, oldScores, denoms, stableScores;
int outerLoops, innerLoops;  
vector<float> old_weights;

void OnlinePAAlg::train() {
  int epochs = atoi(params_->getParamValue("epochs").c_str());
  time_t start, end;
  time(&start);
  // run initial weight training 
  int numTrain = numTrainExamples(), modv = numTrain * .1; 
  int counter(0);
  for(int i=0; i < epochs; ++i) {
    counter = 0;
    iterate(train_set_, tr) {
      optimizeWeights(tr->first, tr->second, false);
      if(++counter >= numTrain) break;
    }
  }
  if(params_->getBoolValue("gr")) {
    counter = 0;
    iterate(test_set_, tr) {
      graphRegularizer(tr->first);
      optimizeWeights(tr->first, scores[tr->first], true);
      if(++counter >= numTrain) break;
    }
    /* // To use all data (labeled and unlabeled)
    for(int i=0; i < numTrain; ++i) {
      if(i % modv == 0)
        cerr << i << endl;
      graphRegularizer(i);
      optimizeWeights(i, scores[i], true);
    }*/
  }
  time(&end);
  double elapTime = difftime(end, start);
  cerr << "Time taken was " << elapTime << " seconds\n";
}
void OnlinePAAlg::graphRegularizer(int nodeI) {
  stableScores = scores; // copy for convergence guarantee
  innerLoops=0;
  float v_i = computeV_i(nodeI); 
  do {
    scores[nodeI] = v_i + computeMSProd(nodeI);
  } while(!innerConvergence());
}
void OnlinePAAlg::optimizeWeights(int node, float corr, bool regress) {
  //bool useSlacks(false);
  // find correct label in train or not
  //map<int, LABEL_T>::const_iterator lit = train_set_.find(node); 
  //if(lit != train_set_.end()) { 
    //useSlacks = true;
  //}
  if(!regress) {
    binaryUpdate(node, sign(corr));
  }
  else {
    regressionUpdate(node, corr, true);
  }
}
void OnlinePAAlg::binaryUpdate(int node, LABEL_T corrLbl) {
  const vector<float>& nodeFeats = data_.at(node);  
  float scalar = dotProduct(nodeFeats);
  float loss = std::max(0.01f, 1.0f - ((float)corrLbl * scalar));
  //cerr << "scalar = " << scalar << "\tloss = " << loss << "\tcorrect = " << corrLbl << endl;
  if(loss > 0.01) updateWeights(nodeFeats, loss, corrLbl, true);
}
void OnlinePAAlg::regressionUpdate(int node, float corrScr, bool useSlacks) {
  static const float epsilon = atof(params_->getParamValue("epsilon").c_str()); 
  const vector<float>& nodeFeats = data_.at(node);  
  float scalar = dotProduct(nodeFeats);
  float absVal = fabs(scalar - corrScr); 
  float loss = (absVal <= epsilon) ? 0 : absVal - epsilon; 
  if(loss > 0) updateWeights(nodeFeats, loss, sign(corrScr - scalar), useSlacks);
}
void OnlinePAAlg::updateWeights(const vector<float>& x, float loss, 
  LABEL_T label, bool useSlacks) {
  static const float c = C_;  //remember initial value
  assert(weights_.size() == x.size());
  if(useSlacks) {
    if(label == NONSPAM) C_ = C_ns;
    else if(label == SPAM) C_ = C_s;
  }
  else C_ = c;
  float tau = getTau(loss, x);
  //cerr << "Tau = " << tau << endl;
  for(int i=0; i < (int)x.size(); ++i) {
    weights_[i] += (tau * label * x[i]); 
  }
  ++numUpdates_;
  cumm_weights_ += weights_;
}
LABEL_T OnlinePAAlg::predict(int node, float* scalar) {
  vector<float>& nodeFeats = data_.at(node);
  if(params_->getBoolValue("gr")) {
    graphRegularizer(node);
    map<int, float>::iterator it = scores.find(node);
    assert(it != scores.end());
    optimizeWeights(node, it->second, true);
  }
  *scalar = dotProduct2(nodeFeats);
  return sign(*scalar);
}
float OnlinePAAlg::computeV_i(int nodeID) {
  map<int, LABEL_T>::const_iterator lbl = train_set_.find(nodeID);
  float l = train_set_.size();
  float b_i = lbl != train_set_.end() ? 1 / l : 0; // b_i term only if labeled node 
  float numer = (b_i * lbl->second) + (lambda_ * dotProduct(data_[nodeID])); 
  float edgeWeights = biDirLinkWeights(nodeID);
  float denom = b_i + lambda_ + (gamma_ * edgeWeights);  
  denoms[nodeID] = denom;
  return numer / denom;
}
float OnlinePAAlg::computeMSProd(int nodeI) {
  map<int, float>::iterator sitr;
  const map<int, int>* edges_ij = graph_->getOutEdges(nodeI);
  float s_i(0);
  if(edges_ij) {
    piterate(edges_ij, eit) {
      int nodeJ = eit->first;
      sitr = oldScores.find(nodeJ);
      if(sitr != oldScores.end()) {  // needed while testing subset 
        float M_ij = computeM_ij(nodeI, nodeJ, eit->second); 
        s_i += sitr->second * M_ij; 
      }
    }
  }
  return s_i;
}
float OnlinePAAlg::computeM_ij(int nodeI, int nodeJ, int edgeCnt_ij) {
  float a_ij = graph_->edgeWeights(edgeCnt_ij); 
  float a_ji = graph_->edgeWeights(graph_->findEdges(nodeJ, nodeI)); 
  bool reg = getScore(nodeI) < getScore(nodeJ);
  float edgeWeight = reg ? 0.1 * (a_ij + a_ji) : a_ij + a_ji; 
  float numer = gamma_ * edgeWeight;
  return numer / denoms[nodeI];
}
float OnlinePAAlg::biDirLinkWeights(const int node) {
  float edgesum(0); // get a_ij + a_ji
  const map<int, int>* edges = graph_->getOutEdges(node); 
  if(edges) { 
    float scr_i = getScore(node);
    // for each edge get weight in both directions
    piterate(edges, edge_itr) {
      float a_ij = graph_->edgeWeights(edge_itr->second); 
      float a_ji = graph_->edgeWeights(graph_->findEdges(edge_itr->first, node));
      float scr_j = getScore(edge_itr->first);
      edgesum += scr_i < scr_j ? 0.1 * (a_ij + a_ji) : a_ij + a_ji;
    }
  }
  return edgesum;
}
float OnlinePAAlg::getScore(int node) {
  map<int, float>::const_iterator sit = stableScores.find(node);
  if(sit != stableScores.end() /*&& sit->second != 0*/)
    return sit->second;
  else if(node < (int)data_.size())
    return dotProduct2(data_[node]);  
  else
    return 0;
}
bool OnlinePAAlg::outerConvergence() {
  assert(old_weights.size() == weights_.size());
  for(size_t i=0; i < weights_.size(); ++i) {
    int old = int(old_weights[i]  * 1e6);
    int wei = int(weights_[i] * 1e6);
    //if(weights_[i] != old_weights[i]) {
      //cerr << i << ": " << weights_[i] << "  vs  " << old_weights[i] << endl;
    if(old != wei) {
      //cerr << i << ": " << wei << "  vs  " << old << endl;
      old_weights = weights_;
      return false;
    }
  }
  return true;
}
bool OnlinePAAlg::innerConvergence() {
  // check if scores have converged
  int max = 0;
  if(innerLoops++ == 0) {
    oldScores = scores;
    return false; 
  }
  //cerr << "Been in convered() " << innerLoops << " times\n";
  if(innerLoops > 400) { cerr << "NOT CONVERGING...EXITING\n"; exit(1); } 
  map<int, float>::iterator oldItr;
  iterate(scores, newItr) {
    //if(newItr->first == 5820) continue;
    oldItr = oldScores.find(newItr->first);
    assert(oldItr != oldScores.end());
    if(oldItr->second != newItr->second) {
      //cerr << std::setprecision(9) << "oldScore[" << oldItr->first << "] = " << oldItr->second;
      //cerr << "\tnewScore[" << newItr->first << "] = " << newItr->second << endl;
      oldScores = scores;
      return false;
    }
    if(++max == 200) { break; } 
  }
  oldScores = scores;
  return true;
}
void OnlinePAAlg::batchTrain() {
  cerr << "C_spam = " << C_s << "\tC_nonspam = " << C_ns << endl;
  time_t start, end;
  time(&start);
  batchOptimizeWeights(false); // get initial weights
  outerLoops=0;
  old_weights = weights_; // remember initial weights 
  do {
    ++outerLoops;
    if(outerLoops % 100 == 0) 
      cerr << "Starting outer convergence loop # " << outerLoops << endl;
    graphRegularizer((outerLoops==1));
    batchOptimizeWeights(true);
  } while(!outerConvergence());
  time(&end);
  double elapTime = difftime(end, start);
  cerr << "Time taken was " << elapTime << " seconds\n";
}
void OnlinePAAlg::batchOptimizeWeights(bool regress) {
  if(!regress) { // binary PA Alg
    iterate(train_set_, trItr) { // for each training instance
      binaryUpdate(trItr->first, trItr->second);
    }
  }
  else { // use regression PA alg
    if(scores.size()) {
      iterate(scores, sit) {
        regressionUpdate(sit->first, scores[sit->first], false);
      }
    }
    else {
      iterate(train_set_, it)
        regressionUpdate(it->first, (float)it->second, true);
    }
  }
}
void OnlinePAAlg::graphRegularizer(bool init) {
  V.clear();
  innerLoops=0;
  static int numTrain = numTrainExamples(); 
  int counter(0);
  iterate(train_set_, tr) {
    ++counter;
    V[tr->first] = computeV_i(tr->first);
    if(init) 
      scores[tr->first] = 0;
    if(counter >= numTrain) 
      break;
  }
  do {
    iterate(V, v) {
      assert(scores.find(v->first) != scores.end());
      scores[v->first] = v->second + computeMSProd(v->first); 
    }
  } while(!innerConvergence());
  //cerr << "\tInner Convergence in " << innerLoops << " iterations\n";
  stableScores = scores; // copy for convergence guarantee
}
