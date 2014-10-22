#ifndef dist_solid_h 
#define dist_solid_h

//#include <boost/math/special_functions/digamma.hpp>
#include "header.h"
#include <numeric> 

class Dist {
public:
  static double digamma(double x) {
    double result = 0, xx, xx2, xx4;
    assert(x > 0);
    for ( ; x < 7; ++x)
    result -= 1/x;
    x -= 1.0/2.0;
    xx = 1.0/x;
    xx2 = xx*xx;
    xx4 = xx2*xx2;
    result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
    return result;
  }
  static double psi(double x) {
    return Dist::digamma(x);
    //return boost::math::digamma(x);
  }
  static double ldirichlet(vector<double>& x, vector<double>& a) {
    assert(x.size() == a.size());
    assert(accumulate(x.begin(), x.end(), 0.0) == 1);
    for(unsigned i=0; i < x.size(); ++i) {
      assert((x[i] > 0) && (x[i] < 1));
      assert(a[i] > 0); 
    }
    double num(0);
    for(unsigned i=0; i < a.size(); ++i) {
      num += lgamma(a[i]);
    }
    double denom = lgamma(accumulate(a.begin(), a.end(), 0.0));
    double z = -(num - denom);
    double prod(0);
    for(unsigned i=0; i < x.size(); ++i) {
      prod += log(x[i]) * (a[i] - 1);
    }
    return prod + z; 
  }
  /* lbetadist() returns the log probability density of x under a Beta(alpha,beta)
   * distribution. - copied from Mark Johnson's gammadist.c */
  static double lbetadist(double x, double alpha, double beta) {
    assert(x > 0);
    assert(x < 1);
    assert(alpha > 0);
    assert(beta > 0);
    return (alpha-1)*log(x)+(beta-1)*log(1-x)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta);
  }
  static double lgammadist(long double x, long double alpha, long double beta) {
    assert(alpha > 0);
    assert(beta > 0);
    if(x==0)
      return (alpha-1) - alpha*log(beta) - lgamma(alpha);
    return (alpha-1)*log(x) - alpha*log(beta) - x/beta - lgamma(alpha);
  }
  static double log_factorial(size_t n) { 
    // return log of factorial
    return lgamma(n + 1);
    //return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }
  static double geomDist(float k, float p) {
    //k = num getting prob for, p = prob of success
    return pow(1.0f-p, k - 1.0f) * (p);
  }
  static double log_geomDist(float k, float p) {
    return log(1.0f-p) * (k-1) + log(p); 
  }
  static double poissDist(size_t k, float lambda) {
    // k = number getting probability for, lambda = mean of distribution
    double num = pow(lambda, k) * exp(-lambda);
    double denom = exp(log_factorial(k)); 
    return num / denom;
  }
  static double log_poissDist(size_t k, float lambda) {
    // k = number getting probability for, lambda = mean of distribution
    double num = (log(lambda) * k)  - lambda;
    double denom = log_factorial(k); 
    return num - denom;
  }
};
#endif
