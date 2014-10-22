#ifndef _TextStats_h_
#define _TextStats_h_

#include <iostream>
#include <cstdlib>   /* for atof() */
#include <cmath>


#ifndef M_LN10
#define M_LN10  2.30258509299404568402
#endif

using namespace std;

/*
 * Types
 */
typedef float LogP;   /* A log-base-10 probability */
typedef double LogP2;   /* A log-base-10 probability, double-size */
typedef double Prob;    /* a straight probability */

/*
 * Constants
 */
extern const LogP LogP_Zero;    /* log(0) = -Infinity */
extern const LogP LogP_Inf;   /* log(Inf) = Infinity */
extern const LogP LogP_One;   /* log(1) = 0 */


inline Prob LogPtoProb(LogP2 prob)
{
  if (prob == LogP_Zero) {
    return 0;
  } else {
    return exp(prob * M_LN10);
  }
}
inline Prob LogPtoPPL(LogP prob)
{
  return exp(- prob * M_LN10);
}

class TextStats
{
public:
    TextStats() : prob(0.0), zeroProbs(0),
	numSentences(0), numWords(0), numOOVs(0) {};

    void reset() { prob = 0.0, zeroProbs = 0,
	numSentences = numWords = numOOVs = 0; };
    TextStats &increment(const TextStats &stats);

    LogP prob;
    unsigned zeroProbs;
    unsigned numSentences;
    unsigned numWords;
    unsigned numOOVs;
};

ostream &operator<<(ostream &, const TextStats &stats);

#endif /* _TextStats_h_ */

