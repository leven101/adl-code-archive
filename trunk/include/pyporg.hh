#ifndef _pyporg_hh
#define _pyporg_hh

#include <math.h>
#include <map>
#include <tr1/unordered_map>
#include <tr1/tuple>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "log_add.h"
#include "mt19937ar.h"
#include "slice_sampler.h"

using std::tr1::hash;
using std::tr1::unordered_map;
using std::tr1::tuple;

//
// Pitman-Yor process with customer and table tracking
//
template <typename Dish, typename Hash=hash<Dish> >
class PYPORG : protected unordered_map<Dish, int, Hash>
{
public:
  using unordered_map<Dish,int>::const_iterator;
  using unordered_map<Dish,int>::iterator;
  using unordered_map<Dish,int>::begin;
  using unordered_map<Dish,int>::end;

  PYPORG(double a, double b, unsigned long seed = 0, Hash hash=Hash());
  ~PYPORG() {}

  int increment(Dish d, double p0);
  int decrement(Dish d);

  double pseudo_increment(Dish d, double p0, double count=1.0);
  double pseudo_prob(Dish d, double p0);
  void clear_expected_counts();

  // lookup functions
  int count(Dish d) const;
  double prob(Dish dish, double p0) const;
  double prob(Dish dish, double dcd, double dca, 
              double dtd, double dta, double p0) const;
  double unnormalised_prob(Dish dish, double p0) const;
  tuple<double,double,double> 
    component_probs(Dish dish, double p0, double dcd=0.0, double dca=0.0, 
                    double dtd=0.0, double dta=0.0) const;

  int num_customers() const { return _total_customers; }
  int num_types() const { return unordered_map<Dish,int>::size(); }
  bool empty() const { return _total_customers == 0; }

  double log_prob(Dish dish, double log_p0) const;
  // nb. d* are NOT logs
  double log_prob(Dish dish, double dcd, double dca, 
                       double dtd, double dta, double log_p0) const;

  int num_tables(Dish dish) const;
  int num_tables() const;

  double a() const { return _a; }
  void set_a(double a) { _a = a; }

  double b() const { return _b; }
  void set_b(double b) { _b = b; }

  void clear();
  std::ostream& debug_info(std::ostream& os) const;

  double log_restaurant_prob() const;
  double log_prior() const;
  static double log_prior_a(double a, double beta_a, double beta_b);
  static double log_prior_b(double b, double gamma_c, double gamma_s);

  template <typename Uniform01>
    void resample_prior(Uniform01& rnd);
  template <typename Uniform01>
    void resample_prior_a(Uniform01& rnd);
  template <typename Uniform01>
    void resample_prior_b(Uniform01& rnd);

protected:
  double _a, _b; // parameters of the Pitman-Yor distribution
  double _a_beta_a, _a_beta_b; // parameters of Beta prior on a
  double _b_gamma_s, _b_gamma_c; // parameters of Gamma prior on b

  struct TableCounter {
    TableCounter() : tables(0) {};
    int tables;
    std::map<int, int> table_histogram; // num customers at table -> number tables
  };
  typedef unordered_map<Dish, TableCounter, Hash> DishTableType;
  DishTableType _dish_tables;
  int _total_customers, _total_tables;

  typedef boost::mt19937 base_generator_type;
  typedef boost::uniform_real<> uni_dist_type;
  typedef boost::variate_generator<base_generator_type&, uni_dist_type> gen_type;

  double m_expected_customer_counts;
  double m_expected_table_counts;
  std::map<Dish, double> m_expected_dish_counts;
  std::map<Dish, double> m_expected_dish_table_counts;

//  uni_dist_type uni_dist;
//  base_generator_type rng; //this gets the seed
//  gen_type rnd; //instantiate: rnd(rng, uni_dist)
                //call: rnd() generates uniform on [0,1)

  // Function objects for calculating the parts of the log_prob for 
  // the parameters a and b
  struct resample_a_type {
    int n, m; double b, a_beta_a, a_beta_b;
    const DishTableType& dish_tables;
    resample_a_type(int n, int m, double b, double a_beta_a, 
                    double a_beta_b, const DishTableType& dish_tables)
      : n(n), m(m), b(b), a_beta_a(a_beta_a), a_beta_b(a_beta_b), dish_tables(dish_tables) {}

    double operator() (double proposed_a) const {
      double log_prior = log_prior_a(proposed_a, a_beta_a, a_beta_b);
      double log_prob = 0.0;
      double lgamma1a = lgamma(1.0 - proposed_a);
      for (typename DishTableType::const_iterator dish_it=dish_tables.begin(); dish_it != dish_tables.end(); ++dish_it) 
        for (std::map<int, int>::const_iterator table_it=dish_it->second.table_histogram.begin(); 
             table_it !=dish_it->second.table_histogram.end(); ++table_it) 
          log_prob += (table_it->second * (lgamma(table_it->first - proposed_a) - lgamma1a));

      log_prob += (proposed_a == 0.0 ? (m-1.0)*log(b) 
                   : ((m-1.0)*log(proposed_a) + lgamma((m-1.0) + b/proposed_a) - lgamma(b/proposed_a)));
      assert(std::isfinite(log_prob));
      return log_prob + log_prior;
    }
  };

  struct resample_b_type {
    int n, m; double a, b_gamma_c, b_gamma_s;
    resample_b_type(int n, int m, double a, double b_gamma_c, double b_gamma_s)
      : n(n), m(m), a(a), b_gamma_c(b_gamma_c), b_gamma_s(b_gamma_s) {}

    double operator() (double proposed_b) const {
      double log_prior = log_prior_b(proposed_b, b_gamma_c, b_gamma_s);
      double log_prob = 0.0;
      log_prob += (a == 0.0  ? (m-1.0)*log(proposed_b) 
                  : ((m-1.0)*log(a) + lgamma((m-1.0) + proposed_b/a) - lgamma(proposed_b/a)));
      log_prob += (lgamma(1.0+proposed_b) - lgamma(n+proposed_b));
      return log_prob + log_prior;
    }
  };
   
  /* lbetadist() returns the log probability density of x under a Beta(alpha,beta)
   * distribution. - copied from Mark Johnson's gammadist.c
   */
  static long double lbetadist(long double x, long double alpha, long double beta);

  /* lgammadist() returns the log probability density of x under a Gamma(alpha,beta)
   * distribution - copied from Mark Johnson's gammadist.c
   */
  static long double lgammadist(long double x, long double alpha, long double beta);

};

template <typename Dish, typename Hash>
PYPORG<Dish,Hash>::PYPORG(double a, double b, unsigned long seed, Hash)
: unordered_map<Dish, int, Hash>(10), _a(a), _b(b), 
  //_a_beta_a(1), _a_beta_b(1), _b_gamma_s(1), _b_gamma_c(1),
  _a_beta_a(1), _a_beta_b(1), _b_gamma_s(10), _b_gamma_c(0.1),
  _total_customers(0), _total_tables(0)//,
  //uni_dist(0,1), rng(seed == 0 ? (unsigned long)this : seed), rnd(rng, uni_dist)
{
  m_expected_customer_counts=0.0; 
  m_expected_table_counts=0.0;
//  std::cerr << "\t##PYPORG<Dish,Hash>::PYPORG(a=" << _a << ",b=" << _b << ")" << std::endl;
  //set_deleted_key(-std::numeric_limits<Dish>::max());
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::prob(Dish dish, double p0) const
{
  assert (p0 >= 0);
  assert (p0 <= 1);
  int c = count(dish), t = num_tables(dish), nc = num_customers();
  double r = num_tables() * _a + _b;
  //std::cerr << "\t\t\t\tPYPORG<Dish,Hash>::prob(" << dish << "," << p0 << ") c=" << c << " r=" << r << std::endl;
  if (nc == 0) return p0;
  else if (c > 0)
    return (c - _a * t + r * p0) / (nc + _b);
  else
    return r * p0 / (nc + _b);
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::unnormalised_prob(Dish dish, double p0) const
{
  int c = count(dish), t = num_tables(dish), nc=num_customers();
  double r = num_tables() * _a + _b;
  if (nc == 0)    return p0;
  else if (c > 0) return (c - _a * t + r * p0);
  else            return r * p0;
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::prob(Dish dish, double dcd, double dca, 
                     double dtd, double dta, double p0) const {
  double c = count(dish) + dcd, t = num_tables(dish) + dtd, nc=num_customers() + dca;
  assert(c >= t);
  double r = (num_tables() + dta) * _a + _b;
  if (nc == 0) return p0;
  else if (c > 0)
    return (c - _a * t + r * p0) / (nc + _b);
  else
    return r * p0 / (nc + _b);
}

template <typename Dish, typename Hash>
tuple<double,double,double>
PYPORG<Dish,Hash>::component_probs(Dish dish, double p0,
                               double dcd, double dca, 
                               double dtd, double dta) const {
  double c = count(dish) + dcd, t = num_tables(dish) + dtd, nc=num_customers() + dca;
  assert(c >= t);

  double r = (num_tables() + dta) * _a + _b;
  double cache = (c == 0.0 ? 0.0 : (c - _a * t)); 
  double new_table = r * p0;
  double denominator = nc + _b;

  if (cache + new_table > denominator) {
    std::cout << std::endl;
    std::cout << count(dish) << " " << dcd << std::endl;
    std::cout << num_tables(dish) << " " << dtd << std::endl;
    std::cout << num_customers() << " " << dca << std::endl;
    std::cout << num_tables() << " " << dta << std::endl;
  }
  
  return std::tr1::make_tuple(cache,new_table,denominator);
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_prob(Dish dish, double log_p0) const
{
  using std::log;
  int c = count(dish), t = num_tables(dish), nc=num_customers();
  double r = log(num_tables() * _a + _b);
  if (nc == 0) { 
    return log_p0;
  }
  else if (c > 0) {
    return Log<double>::add(log(c - _a * t), r + log_p0)
      - log(num_customers() + _b);
  }
  else {
    return r + log_p0 - log(num_customers() + _b);
  }
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_prob(Dish dish, double dcd, double dca, 
                         double dtd, double dta, double log_p0)
const
{
  using std::log;
  int c = count(dish) + dcd, t = num_tables(dish) + dtd, nc=num_customers();
  double r = log((num_tables() + dta) * _a + _b);
  if (nc == 0) return log_p0;
  else if (c > 0)
    return Log<double>::add(log(c - _a * t), r + log_p0)
      - log(num_customers() + dca + _b);
  else
    return r + log_p0 - log(num_customers() + dca + _b);
}

template <typename Dish, typename Hash>
int 
PYPORG<Dish,Hash>::increment(Dish dish, double p0) {
  assert (p0 >= 0);
  assert (p0 <= 1);

  int delta = 0;
  TableCounter &tc = _dish_tables[dish];

  // seated on a new or existing table?
  int c = count(dish), t = num_tables(dish), T = num_tables();
  assert(c >= t);
  double pshare = (c > 0) ? (c - _a*t) : 0.0;
  double pnew = (c > 0) ? (_b + _a*T) * p0 : 1.0;
  assert (pshare >= 0.0);
  //assert (pnew > 0.0);

  //if (rnd() < pnew / (pshare + pnew)) {
  if (mt_genrand_res53() < pnew / (pshare + pnew)) {
    // assign to a new table
    tc.tables += 1;
    tc.table_histogram[1] += 1;
    _total_tables += 1;
    delta = 1;
  }
  else {
    // randomly assign to an existing table
    // remove constant denominator from inner loop
    //double r = rnd() * (c - _a*t);
    double r = mt_genrand_res53() * (c - _a*t);
    for (std::map<int,int>::iterator
         hit = tc.table_histogram.begin();
         hit != tc.table_histogram.end(); ++hit) {
      r -= ((hit->first - _a) * hit->second);
      if (r <= 0) {
        tc.table_histogram[hit->first+1] += 1;
        hit->second -= 1;
        if (hit->second == 0)
          tc.table_histogram.erase(hit);
        break;
      }
    }
    if (r > 0) {
      std::cerr << r << " " << c << " " << _a << " " << t << " " << T << " " << p0 <<  std::endl;
      assert(false);
    }
    delta = 0;
  }

  unordered_map<Dish,int,Hash>::operator[](dish) += 1;
  _total_customers += 1;

  return delta;
}

// Increment and probability calculations based on expected table counts.
// Use in conjunction with clear_expected_counts to approxinate probs
// for multiple dependent increments.
template <typename Dish, typename Hash>
double PYPORG<Dish,Hash>::pseudo_prob(Dish dish, double p0) {
  double dcd = m_expected_dish_counts[dish]; 
  double dtd = m_expected_dish_table_counts[dish]; 

  double c = count(dish) + dcd, t = num_tables(dish) + dtd, 
         nc=num_customers() + m_expected_customer_counts;
  assert(c >= t);

  double r = (num_tables() + m_expected_table_counts) * _a + _b;
  double cache = (c == 0.0 ? 0.0 : (c - _a * t)); 
  double new_table = r * p0;
  double denominator = nc + _b;
  double numerator = cache+new_table;

  // get the components of the probability computation for this dish
//  std::tuple<double,double,double> components 
//    = component_probs(dish, p0,
//                      m_expected_dish_counts[dish], 
//                      m_expected_customer_counts,
//                      m_expected_dish_table_counts[dish],
//                      m_expected_table_counts);

  // return the probability of the dish before the increment
//  if (numerator > denominator)
//    std::cout << numerator << " " << denominator << " " << prob(dish, p0) << " " << p0 << std::endl;
//  assert (numerator <= denominator);
//  assert (numerator > 0.0);
//  assert (denominator > 0.0 );
  return numerator / denominator;
}
template <typename Dish, typename Hash>
double PYPORG<Dish,Hash>::pseudo_increment(Dish dish, double p0, double inc_count) {
  double& dcd = m_expected_dish_counts[dish]; 
  double& dtd = m_expected_dish_table_counts[dish]; 
  assert(inc_count > 0.0 && inc_count <= 1.0);

  double c = count(dish) + dcd, t = num_tables(dish) + dtd;
  assert(c >= t);

  double r = (num_tables() + m_expected_table_counts) * _a + _b;
  double cache = (c == 0.0 ? 0.0 : (c - _a * t)); 
  double new_table = r * p0;
  //double denominator = nc + _b;
  double numerator = cache+new_table;
  // get the components of the probability computation for this dish
  //std::tuple<double,double,double> components 
  //  = component_probs(dish, p0,
  //                    m_expected_dish_counts[dish], 
  //                    m_expected_customer_counts,
  //                    m_expected_dish_table_counts[dish],
  //                    m_expected_table_counts);

  //double numerator = std::get<0>(components)+std::get<1>(components);
  //double denominator = std::get<2>(components);
  
  // calculate the expected probability that this increment created a new table,
  // scaled by the increment count.
  //double new_table_prob = count * (std::get<1>(components) / numerator);
  double new_table_prob = inc_count * (new_table / numerator);
  //assert(new_table_prob > 0.0 && new_table_prob <= 1.0);
  //m_expected_dish_counts[dish] += count;
  //m_expected_dish_table_counts[dish] += new_table_prob;
  dcd += inc_count;
  dtd += new_table_prob;
  m_expected_customer_counts += inc_count;
  m_expected_table_counts += new_table_prob;
  //assert(m_expected_dish_counts[dish] <= m_expected_customer_counts);
  //assert(m_expected_dish_table_counts[dish] <= m_expected_table_counts);

  // return the probability of the dish before the increment
  //assert (numerator <= denominator);
  //assert (numerator > 0.0);
  //assert (denominator > 0.0 );
  return new_table_prob;
}

template <typename Dish, typename Hash>
void PYPORG<Dish,Hash>::clear_expected_counts() {
  m_expected_customer_counts = 0.0;
  m_expected_table_counts = 0.0;
  m_expected_dish_counts.clear();
  m_expected_dish_table_counts.clear();
}

template <typename Dish, typename Hash>
int 
PYPORG<Dish,Hash>::count(Dish dish) const
{
  typename unordered_map<Dish, int>::const_iterator 
    dcit = find(dish);
  if (dcit != end())
    return dcit->second;
  else
    return 0;
}

template <typename Dish, typename Hash>
int 
PYPORG<Dish,Hash>::decrement(Dish dish)
{
  typename unordered_map<Dish, int>::iterator dcit = find(dish);
  if (dcit == end()) {
    std::cerr << dish << std::endl;
    assert(false);
  } 

  int delta = 0;

  typename unordered_map<Dish, TableCounter>::iterator dtit = _dish_tables.find(dish);
  if (dtit == _dish_tables.end()) {
    std::cerr << dish << std::endl;
    assert(false);
  } 
  TableCounter &tc = dtit->second;

  //std::cerr << "\tdecrement for " << dish << "\n";
  //std::cerr << "\tBEFORE histogram: " << tc.table_histogram << " ";
  //std::cerr << "count: " << count(dish) << " ";
  //std::cerr << "tables: " << tc.tables << "\n";

  //double r = rnd() * count(dish);
  double r = mt_genrand_res53() * count(dish);
  for (std::map<int,int>::iterator hit = tc.table_histogram.begin();
       hit != tc.table_histogram.end(); ++hit)
  {
    r -= (hit->first) * hit->second;
    if (r <= 0)
    {
      if (hit->first > 1)
        tc.table_histogram[hit->first-1] += 1;
      else
      {
        delta = -1;
        tc.tables -= 1;
        _total_tables -= 1;
      }

      hit->second -= 1;
      if (hit->second == 0) tc.table_histogram.erase(hit);
      break;
    }
  }
  if (r > 0) {
    std::cerr << r << " " << count(dish) << " " << _a << " " << num_tables(dish) << std::endl;
    assert(false);
  }

  // remove the customer
  dcit->second -= 1;
  _total_customers -= 1;
  assert(dcit->second >= 0);
  if (dcit->second == 0) {
    erase(dcit);
    _dish_tables.erase(dtit);
    //std::cerr << "\tAFTER histogram: Empty\n";
  }
  else {
    //std::cerr << "\tAFTER histogram: " << _dish_tables[dish].table_histogram << " ";
    //std::cerr << "count: " << count(dish) << " ";
    //std::cerr << "tables: " << _dish_tables[dish].tables << "\n";
  }

  return delta;
}

template <typename Dish, typename Hash>
int 
PYPORG<Dish,Hash>::num_tables(Dish dish) const
{
  typename unordered_map<Dish, TableCounter, Hash>::const_iterator 
    dtit = _dish_tables.find(dish);

  //assert(dtit != _dish_tables.end());
  if (dtit == _dish_tables.end())
    return 0;

  return dtit->second.tables;
}

template <typename Dish, typename Hash>
int 
PYPORG<Dish,Hash>::num_tables() const
{
  return _total_tables;
}

template <typename Dish, typename Hash>
std::ostream&
PYPORG<Dish,Hash>::debug_info(std::ostream& os) const
{
  int hists = 0, tables = 0;
  for (typename unordered_map<Dish, TableCounter, Hash>::const_iterator 
       dtit = _dish_tables.begin(); dtit != _dish_tables.end(); ++dtit)
  {
    hists += dtit->second.table_histogram.size();
    tables += dtit->second.tables;

//    if (dtit->second.tables <= 0)
//      std::cerr << dtit->first << " " << count(dtit->first) << std::endl;
    assert(dtit->second.tables > 0);
    assert(!dtit->second.table_histogram.empty());

//    os << "Dish " << dtit->first << " has " << count(dtit->first) << " customers, and is sitting at " << dtit->second.tables << " tables.\n"; 
    for (std::map<int,int>::const_iterator 
         hit = dtit->second.table_histogram.begin();
         hit != dtit->second.table_histogram.end(); ++hit) {
//      os << "    " << hit->second << " tables with " << hit->first << " customers." << std::endl; 
      assert(hit->second > 0);
    }
  }

  os << "restaurant has " 
    << _total_customers << " customers; "
    << _total_tables << " tables; " 
    << tables << " tables'; " 
    << num_types() << " dishes; "
    << _dish_tables.size() << " dishes'; and "
    << hists << " histogram entries\n";

  return os;
}

template <typename Dish, typename Hash>
void 
PYPORG<Dish,Hash>::clear()
{
  this->unordered_map<Dish,int,Hash>::clear();
  _dish_tables.clear();
  _total_tables = _total_customers = 0;
}

// log_restaurant_prob returns the log probability of the PYPORG table configuration.
// Excludes Hierarchical P0 term which must be calculated separately.
template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_restaurant_prob() const {
  if (_total_customers < 1)
    return (double)0.0;

  double log_prob = 0.0;
  double lgamma1a = lgamma(1.0-_a);

  //std::cerr << "-------------------\n" << std::endl;
  for (typename DishTableType::const_iterator dish_it=_dish_tables.begin(); 
       dish_it != _dish_tables.end(); ++dish_it) {
    for (std::map<int, int>::const_iterator table_it=dish_it->second.table_histogram.begin(); 
         table_it !=dish_it->second.table_histogram.end(); ++table_it) {
      log_prob += (table_it->second * (lgamma(table_it->first - _a) - lgamma1a));
      //std::cerr << "|" << dish_it->first->parent << " --> " << dish_it->first->rhs << " " << table_it->first << " " << table_it->second << " " << log_prob;
    }
  }
  //std::cerr << std::endl;

  log_prob += (_a == (double)0.0 ? (_total_tables-1.0)*log(_b) : (_total_tables-1.0)*log(_a) + lgamma((_total_tables-1.0) + _b/_a) - lgamma(_b/_a));

  //std::cerr << "\t\t" << log_prob << std::endl;
  log_prob += (lgamma(1.0 + _b) - lgamma(_total_customers + _b));

  //std::cerr << _total_customers << " " << _total_tables << " " << log_prob << " " << log_prior() << std::endl;
  //std::cerr << _a << " " << _b << std::endl;
  if (!std::isfinite(log_prob)) {
    assert(false);
  }
  //return log_prob;
  if (log_prob > 0.0)
    std::cerr << log_prob << std::endl;
  return log_prob + log_prior();
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_prior() const {
  double prior = 0.0;
  if (_a_beta_a > 0.0 && _a_beta_b > 0.0 && _a > 0.0)
    prior += log_prior_a(_a, _a_beta_a, _a_beta_b);
  if (_b_gamma_s > 0.0 && _b_gamma_c > 0.0)
    prior += log_prior_b(_b, _b_gamma_c, _b_gamma_s);

  return prior;
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_prior_a(double a, double beta_a, double beta_b) {
  return lbetadist(a, beta_a, beta_b); 
}

template <typename Dish, typename Hash>
double 
PYPORG<Dish,Hash>::log_prior_b(double b, double gamma_c, double gamma_s) {
  return lgammadist(b, gamma_c, gamma_s); 
}

template <typename Dish, typename Hash>
long double PYPORG<Dish,Hash>::lbetadist(long double x, long double alpha, long double beta) {
  assert(x > 0);
  assert(x < 1);
  assert(alpha > 0);
  assert(beta > 0);
  return (alpha-1)*log(x)+(beta-1)*log(1-x)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta);
//boost::math::lgamma
}

template <typename Dish, typename Hash>
long double PYPORG<Dish,Hash>::lgammadist(long double x, long double alpha, long double beta) {
  assert(alpha > 0);
  assert(beta > 0);
  return (alpha-1)*log(x) - alpha*log(beta) - x/beta - lgamma(alpha);
}


template <typename Dish, typename Hash>
  template <typename Uniform01>
void 
PYPORG<Dish,Hash>::resample_prior(Uniform01& rnd) {
  for (int num_its=5; num_its >= 0; --num_its) {
    if (_b != 0.0) resample_prior_b(rnd);
    if (_a != 0.0) resample_prior_a(rnd);
  }
  if (_b != 0.0) resample_prior_b(rnd);
}

template <typename Dish, typename Hash>
  template <typename Uniform01>
void 
PYPORG<Dish,Hash>::resample_prior_b(Uniform01& rnd) {
  if (_total_tables == 0) 
    return;

  //int niterations = 10;   // number of resampling iterations
  int niterations = 5;   // number of resampling iterations
  //std::cerr << "\n## resample_prior_b(), initial a = " << _a << ", b = " << _b << std::endl;
  resample_b_type b_log_prob(_total_customers, _total_tables, _a, _b_gamma_c, _b_gamma_s);
  _b = slice_sampler1d(b_log_prob, _b, rnd, (double) 0.0, std::numeric_limits<double>::infinity(), 
  //_b = slice_sampler1d(b_log_prob, _b, mt_genrand_res53, (double) 0.0, std::numeric_limits<double>::infinity(), 
                       (double) 0.0, niterations, 100*niterations);
  //std::cerr << "\n## resample_prior_b(), final a = " << _a << ", b = " << _b << std::endl;
}

template <typename Dish, typename Hash>
  template <typename Uniform01>
void 
PYPORG<Dish,Hash>::resample_prior_a(Uniform01& rnd) {
  if (_total_tables == 0) 
    return;

  //int niterations = 10;
  int niterations = 5;
  //std::cerr << "\n## Initial a = " << _a << ", b = " << _b << std::endl;
  resample_a_type a_log_prob(_total_customers, _total_tables, _b, _a_beta_a, _a_beta_b, _dish_tables);
  _a = slice_sampler1d(a_log_prob, _a, rnd, std::numeric_limits<double>::min(), 
  //_a = slice_sampler1d(a_log_prob, _a, mt_genrand_res53, std::numeric_limits<double>::min(), 
                       (double) 1.0, (double) 0.0, niterations, 100*niterations);
}


#endif
