#include "TextStats.h"

const LogP LogP_Zero = -HUGE_VAL;   /* log(0) */
const LogP LogP_Inf = HUGE_VAL;     /* log(Inf) */
const LogP LogP_One = 0.0;      /* log(1) */

/*
 * Increments from other source
 */
TextStats &
TextStats::increment(const TextStats &stats)
{
    numSentences += stats.numSentences;
    numWords += stats.numWords;
    numOOVs += stats.numOOVs;
    prob += stats.prob;
    zeroProbs += stats.zeroProbs;

    return *this;
}

/*
 * Format stats for stream output
 */
ostream &
operator<< (ostream &stream, const TextStats &stats)
{

    stream << stats.numSentences << " sentences, " 
           << stats.numWords << " words, "
	   << stats.numOOVs << " OOVs" << endl;
    if (stats.numWords + stats.numSentences > 0) {
	stream << stats.zeroProbs << " zeroprobs, "
	       << "logprob= " << stats.prob;

	int denom = stats.numWords - stats.numOOVs - stats.zeroProbs
							+ stats.numSentences;

	if (denom > 0) {
	    stream << " ppl= " << LogPtoPPL(stats.prob / denom);
	} else {
	    stream << " ppl= undefined";
	}

	denom -= stats.numSentences;

	if (denom > 0) {
	    stream << " ppl1= " << LogPtoPPL(stats.prob / denom);
	} else {
	    stream << " ppl1= undefined";
	}

	stream << endl;
    }
    return stream;
}

