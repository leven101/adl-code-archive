#ifndef timeseries_solid_h
#define timeseries_solid_h

#include <numeric> 
#include "header.h"

class TS {
public:
  struct Stats {
    f_t mean, var, stddev;
    Stats(const f_t m, const f_t v, const f_t sd) {
      mean = m;
      var = v;
      stddev = sd;
    }
  };
  static Stats sd(const f_t* x, const f_t sz) { 
    // mean
    f_t mean=0;
    for(int i=0; i < sz; ++i) {
      mean += x[i];
    }
    mean /= sz;
    // variance
    f_t var=0;
    for(int i=0; i < sz; ++i) {
      var += pow(x[i]-mean, 2.); 
    }
    var /= (sz-1.0);
    // standard deviation
    f_t sd = pow(var, .5); 
    Stats stats(mean, var, sd);
    return stats;
  }
  static vf_t sd_normalize(const vf_t& x) {
    Stats stats = sd(&x[0], x.size());
    vf_t xnorm(x.size());
    for(size_t i=0; i < x.size(); ++i) {
      xnorm[i] = (x[i] - stats.mean) / stats.stddev;
    }
    return xnorm;
  }
  static f_t* sd_normalize(const f_t* x, const int sz) {
    /* Allocates memory. Need to free later!!! */
    Stats stats = sd(x, sz);
    f_t* xnorm = new f_t[sz]; 
    for(int i=0; i < sz; ++i) {
      xnorm[i] = (x[i] - stats.mean) / stats.stddev;
    }
    return xnorm;
  }
  static vf_t unit_normalize(const vf_t& x) {
    float min = *std::min_element(x.begin(), x.end());
    float max = *std::max_element(x.begin(), x.end());
    vf_t xnorm(x.size());
    for(size_t i=0; i < x.size(); ++i) {
      xnorm[i] = (x[i] - min) / max;
    }
    return xnorm;
  }
  static f_t* unit_normalize(const f_t* x, const int sz) {
    /* Allocates memory. Need to free later!!! */
    f_t* xnorm = new f_t[sz]; 
    f_t min = *std::min_element(x, x + sz);
    f_t max = *std::max_element(x, x + sz);
    for(int i=0; i < sz; ++i) {
      xnorm[i] = (x[i] - min) / max;
    }
    return xnorm;
  }
  static vf_t diff(const vf_t& x, const int lag=1,
    const int levels=1) {
    vf_t xdiff(x.size());
    for(int i=0; i < lag; ++i) {
      xdiff[i]=0;
    }
    for(size_t t=lag; t < x.size(); ++t) {
      xdiff[t] = x[t] - x[t-lag];
    }
    if(levels > 1) {
      xdiff = diff(xdiff, lag, levels-1);
    }
    return xdiff;
  }
  static f_t* diff(const f_t*x, const int sz, const int lag=1,
    const int levels=1) {
    /* Allocates memory. Need to free later!!! */
    const int newsize = sz-lag;
    f_t* xdiff = new f_t[newsize];
    for(int t=lag; t < sz; ++t) {
      xdiff[t-lag] = x[t] - x[t-lag];
    }
    if(levels > 1) {
      //xdiff = diff(xdiff, newsize, lag, levels-1);
      f_t* tmp = diff(xdiff, newsize, lag, levels-1);
      delete [] xdiff;
      //xdiff = new f_t[newsize-lag];
      xdiff = tmp;
    }
    return xdiff;
  }
  static void testDiff(const int lag=1, const int levels=1) {
    vector<f_t> f;
    int csum(0);
    for(int i=1; i <= 10; ++i) {
      csum += i;
      f.push_back(csum);
    }
    for(int i=1; i <= 10; ++i) {
      f[i] = f[i] + f[i-1];
    }
    for(int i=0; i < 10; ++i) 
      cout << f[i] << "  ";
    cout << endl;
    f_t* d = diff(&f[0], f.size(), lag, levels);
    int newsz = f.size()-levels;
    for(int i=0; i < newsz; ++i) 
      cout << d[i] << "  ";
    cout << endl;
    delete [] d;
  }
  static int unique(const f_t* x, const int sz) {
    set<f_t> uniqvals(x, x+sz);
    return uniqvals.size();
  }
  static vf_t binary_band(const vf_t& feats) {
    vf_t tmp(feats.size());
    tmp[0] = 0;
    for(size_t i=1; i < feats.size(); ++i) {
      tmp[i] =  feats[i] <= feats[i-1] ? 0 : 1;
    }
    return tmp;
  }
  static vf_t sma(const vf_t& v, const float n) {
    assert(n < v.size());
    vf_t sma;
    for(int i=0; i < n-1; ++i) {
      sma.push_back(0);
    }
    for(int i=0; i <= v.size()-n; ++i) {
      float sum = std::accumulate(v.begin() + i, v.begin() + (i+n), 0.0f);
      float s = sum / n;
      sma.push_back(s);
    } 
    assert(sma.size() == v.size());
    return sma;
  }
};
#endif
