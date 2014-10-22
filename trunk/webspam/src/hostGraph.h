#ifndef INC_WEBSPAM_HOSTGRAPH_H
#define INC_WEBSPAM_HOSTGRAPH_H
#include "features.h"


typedef const map<int, int>* EDGE_PT;
class HostGraph {
public:
  HostGraph(Parameters& params, Features& feats): params_(&params),
  f_(&feats) {
    loadHostGraph();
    if(params_->getBoolValue("load-hg-feats"))
      loadHostFeatures();
    else 
      cerr << "Not using hostgraph features!!!!!!!!!!\n";
  }
  ~HostGraph() {
  }
  int findEdges(const int nodej, const int); 
  const map<int, int>* getOutEdges(const int);
  const map<int, int>* getInEdges(const int);
  float edgeWeights(int numEdges) {
    return log(1 + numEdges);
    //return sqrt(numEdges);
    //return 1;
  }
  float edgeWeights2(int numEdges) {
    return 1;
  }
  vector<vector< float> > hostFeatures_;
private:
  map<int, map<int, int> > out_links_, in_links_;
  Parameters* const params_;
  Features* const f_;
  void loadHostGraph(); 
  const vector<vector<float> >* loadHostFeatures();
  float measure(const int node, int idx) {
   if(node < (int)f_->linkFeatures_.size())
     return f_->linkFeatures_.at(node).at(idx); 
   return 0;
  }
  float avoidNullVals(float numer, float denom) {
    float res(0);
    if(denom == 0) 
      res = numer != 0 ? 0 : 1;
    else 
      res = numer /denom;
    return res; 
  }
  float gengF_2(int featIdx, EDGE_PT inEdges) {
    float numer(0), denom(0); 
    int totEdges(0);
    //piterate(inEdges, i) { // for each inlink
      //totEdges += i->second;
    //}
    piterate(inEdges, i) { // for each inlink
      // take measure(h) * weight(h, H) = weight on number of links between nodes
      numer += measure(i->first, featIdx) * edgeWeights2(totEdges);
      EDGE_PT outEdges = getOutEdges(i->first);
      if(outEdges) {
        piterate(outEdges, o) {
          denom += edgeWeights2(o->second); // sum of weight(h,g) of all outlinks of h
        }
      }
    }
    return avoidNullVals(numer, denom);
  }
  float gengF_3(int featIdx, EDGE_PT outEdges) {
    float numer(0), denom(0); 
    int totEdges(0);
    //piterate(outEdges, o) { // for each inlink
      //totEdges += o->second;
    //}
    piterate(outEdges, o) { // for each inlink
      // take measure(h) * weight(h, H) = weight on number of links between nodes
      numer += measure(o->first, featIdx) * edgeWeights2(totEdges);
      EDGE_PT inEdges = getInEdges(o->first);
      if(inEdges) {
        piterate(inEdges, i) {
          denom += edgeWeights2(i->second); // sum of weight(h,g) of all outlinks of h
        }
      }
    }
    return avoidNullVals(numer, denom);
  }
  float gengF_4And5(int featIdx, EDGE_PT edges) {
    float numer(0), denom(0); 
    // get the inlinks of the edges 
    piterate(edges, e) {
      EDGE_PT inEdges = getInEdges(e->first);
      if(inEdges) {
        piterate(inEdges, ii) {
          numer += measure(ii->first, featIdx);
          denom += ii->second; // gives |Inlink(Inlink(H))| 
        }
      }
    }
    return avoidNullVals(numer, denom);
  }
  float gengF_6And7(int featIdx, EDGE_PT edges) {
    float numer(0), denom(0); 
    // get the outlinks of the edges 
    piterate(edges, e) {
      EDGE_PT outEdges = getOutEdges(e->first);
      if(outEdges) {
        piterate(outEdges, oo) {
          numer += measure(oo->first, featIdx);
          denom += oo->second; // gives |Outlink(Outlink(H))| 
        }
      }
    }
    return avoidNullVals(numer, denom);
  }
  void normalizeFeatures() {
  cerr << "Normalizing hostGraph features... " << endl;
  const size_t numFeatures(hostFeatures_[0].size());
  string outfile = params_->getParamValue("top-dir") + params_->getParamValue("year");
  outfile += "/features/graphFeatures.normalized";
  f_->normalize(numFeatures, hostFeatures_, outfile);
}
};
#endif
