#ifndef _RULE_HH
#define _RULE_HH

#include "text.hh"
#include "factorise.hh"
#include "alignment.hh"
#include "prior.hh"
#include "serialise.pb.hh"

#include <map>

#include <boost/functional/hash.hpp>

class Rule;
typedef boost::shared_ptr<Rule> RulePtr;

class Rule {
public:
  friend std::size_t 
    hash_value(Rule const& r)
    {
      std::size_t seed = 0;
      boost::hash_combine(seed, r.m_label);
      boost::hash_combine(seed, r.m_source);
      boost::hash_combine(seed, r.m_target);

      return seed;
    }
  enum Type { PRODUCTION, EMISSION, NULL_EMISSION };

public:
  Rule() { /*update_hash();*/ }

  Rule(const IndexedToken& label, const IndexedString& source, const IndexedString& target)
    : m_label(label), m_source(source), m_target(target) { /*update_hash();*/ }

  Rule(const Rule& r)
    : m_label(r.label()), m_source(r.source()), m_target(r.target()) { /*update_hash();*/ }

  Rule(const DecompositionTree& parse_tree, const DecompositionTree::iterator& node_it, 
      const Sentence& source_sentence, const Sentence& target_sentence) {
    int num_children = parse_tree.number_of_children(node_it);

    if (num_children == 0) {
      label(-1);
      source(source_sentence.span(node_it->source_span.start, node_it->source_span.end));
      target(target_sentence.span(node_it->target_span.start, node_it->target_span.end));
    }
    else {
      // sort the phrases on source order (assume they're already sorted on the target)
      std::map<Phrase, int, PhraseLt> source_sorted_phrases; 
      int counter=0;
      m_target.resize(num_children);
      for (DecompositionTree::sibling_iterator child_it=parse_tree.begin(node_it); 
          child_it != parse_tree.end(node_it); ++child_it, ++counter) {
        source_sorted_phrases.insert(std::make_pair(*child_it, counter));
        m_target.at(counter) = (parse_tree.number_of_children(child_it) == 0 ? -1 : -1);;
      }

      m_source.resize(num_children);
      counter=0;
      label(-1);
      for (std::map<Phrase, int, PhraseLt>::iterator child_it=source_sorted_phrases.begin(); 
          child_it != source_sorted_phrases.end(); ++child_it, ++counter) {
        m_source.at(child_it->second) = counter;;
      }
    }

    //update_hash();
  }

  Rule(const std::string& str) {
    // [a]->[c,d|||e,f]
    assert(str[0] == '['); assert(str[str.size()-1] == ']'); 
    std::string::size_type first_close = str.find(']');
    m_label = boost::lexical_cast<int>(str.substr(1, first_close - 1));

    assert(str[first_close+1] == '-'); 
    assert(str[first_close+2] == '>'); 
    assert(str[first_close+3] == '['); 

    m_source.clear();
    std::string::size_type comma = str.find(',');
    std::string::size_type start = first_close + 4;;
    std::string::size_type source_end= str.find('|');

    while (comma != std::string::npos && comma < source_end) {
      assert(comma > start);
      m_source.push_back(boost::lexical_cast<int>(str.substr(start, comma - start)));
      start = comma + 1;
      comma = str.find(',', start);
    }
    if (source_end > start)
      m_source.push_back(boost::lexical_cast<int>(str.substr(start, source_end - start)));

    assert(str[source_end+1] == '|'); 
    assert(str[source_end+2] == '|'); 

    m_target.clear();
    start = source_end + 3;
    comma = str.find(',', start);
    std::string::size_type target_end=str.find(']', start);

    while (comma != std::string::npos && comma < target_end) {
      assert(comma > start);
      m_target.push_back(boost::lexical_cast<int>(str.substr(start, comma - start)));
      start = comma + 1;
      comma = str.find(',', start);
    }
    if (target_end > start)
      m_target.push_back(boost::lexical_cast<int>(str.substr(start, target_end - start)));

    //update_hash();
  }

  Rule(const Rule& first, const Rule& second, bool monotone=true) {
    // copy the source/target in the correct order
    m_target = first.target();
    m_target.insert(m_target.end(), second.target().begin(), second.target().end());

    if (monotone) { // Monotone
      m_source = first.source();
      m_source.insert(m_source.end(), second.source().begin(), second.source().end());
    }
    else { // Reordered
      m_source = second.source();
      m_source.insert(m_source.end(), first.source().begin(), first.source().end());
    }

    m_label = -1; // FIXME
  }

  bool
    is_null() const { return m_source.empty() && m_target.empty(); }

  bool
    is_half_null() const { return m_source.empty() || m_target.empty(); }

  const IndexedToken&
    label() const { return m_label; }

  void
    label(const IndexedToken& l) { m_label = l; }

  const IndexedString&
    source() const { return m_source; }

  void
    source(const IndexedString& s) { m_source = s; }

  const IndexedString&
    target() const { return m_target; }

  void
    target(const IndexedString& t) { m_target = t; }

  // Replace child label
  void
    replace(int target_index, IndexedToken new_node) {
      m_target.at(target_index) = new_node;
    }

  // Replace adjacent pair of nodes with a single new one
  void
    replace(int target_first, int target_second, IndexedToken new_node) {
      assert(std::abs(target_first - target_second) == 1);
      assert(m_source.size() > 2 && m_target.size() > 2);
      assert(std::abs(m_source.at(target_first) - m_source.at(target_second)) == 1);
      assert (m_source.size() == m_target.size());

      int lowest_target_index = std::min(target_first, target_second);
      int highest_target_index = std::max(target_first, target_second);
      int lowest_mapping_index = std::min(m_source.at(target_first), m_source.at(target_second));

      // erase/replace the target pair
      m_target.erase(m_target.begin() + highest_target_index); 
      m_target.at(lowest_target_index) = new_node;
      m_source.erase(m_source.begin() + highest_target_index); 

      // remap source indexes
      for (IndexedString::iterator it=m_source.begin(); it != m_source.end(); ++it) {
        if (*it > lowest_mapping_index)
          --(*it);
      }

      assert (m_source.size() == m_target.size());
      //update_hash();
    }

  void
    replace(int target_index, int new_left, int new_right, bool reordered) {
      int source_index = m_source.at(target_index);  

      // update the target
      m_target.at(target_index) = new_right;
      m_target.insert(m_target.begin() + target_index, new_left);

      // update the source 
      m_source.at(target_index) = source_index + (reordered ? 0 : 1);
      m_source.insert(m_source.begin() + target_index, source_index + (reordered ? 1 : 0));

      // remap source indexes
      int counter=0;
      for (IndexedString::iterator it=m_source.begin(); it != m_source.end(); ++it, ++counter) {
        if (counter != target_index && counter != target_index+1 && *it > source_index)
          (*it)++;
      }

      //update_hash();
    }

  // return all splits of this rules source and target
  std::pair< std::pair<RulePtr, RulePtr>, std::pair<RulePtr, RulePtr> >
    split(int s, int t, bool include_monotone, bool include_reordered) const {
      assert(s <= (int)m_source.size() && s >= 0);
      assert(t <= (int)m_target.size() && t >= 0);
      IndexedString left_source_segment(m_source.begin(), m_source.begin() + s); 
      IndexedString left_target_segment(m_target.begin(), m_target.begin() + t);

      IndexedString right_source_segment, right_target_segment;
      if (s < (int)m_source.size()) 
        right_source_segment.insert(right_source_segment.begin(), m_source.begin() + s, m_source.end()); 
      if (t < (int)m_target.size()) 
        right_target_segment.insert(right_target_segment.begin(), m_target.begin() + t, m_target.end());
      
      std::pair<RulePtr,RulePtr> monotone_pair, reordered_pair;
      if (include_monotone)
        monotone_pair = std::make_pair(new Rule(m_label, left_source_segment, left_target_segment), 
                                       new Rule(m_label, right_source_segment, right_target_segment));
      if (include_reordered)
        reordered_pair = std::make_pair(new Rule(m_label, right_source_segment, left_target_segment), 
                                        new Rule(m_label, left_source_segment, right_target_segment));

      return std::make_pair(monotone_pair, reordered_pair);
    }
  
  Type
    type() const
    {
      bool found_terminal=false, found_non_terminal=false;
      for (IndexedString::const_iterator it=m_target.begin(); it != m_target.end(); ++it) {
        if (*it < 0) found_non_terminal = true;
        else found_terminal = true;
      }
      assert (found_terminal || found_non_terminal || m_target.empty());
      assert (!(found_terminal && found_non_terminal));
      return (found_non_terminal ? PRODUCTION : EMISSION);
    }

  size_t hash() const { return hash_value(*this); }
  size_t source_hash() const { return hash_value(m_source); }
  size_t target_hash() const { return hash_value(m_target); }

//  void update_hash() { 
//    m_source_hash = hash_value(m_source); 
//    m_target_hash = hash_value(m_target); 
//    m_hash = hash_value(*this); 
//  }

  void
    to_message(SerialisedRule& sr) const {
      sr.set_label(m_label);
      for (IndexedString::const_iterator it=m_source.begin(); it!=m_source.end(); ++it)
        sr.add_source(*it);
      for (IndexedString::const_iterator it=m_target.begin(); it!=m_target.end(); ++it)
        sr.add_target(*it);
    }

  void
    from_message(const SerialisedRule& sr) {
      m_label = sr.label();

      m_source.resize(sr.source_size());
      for (int i=0; i < sr.source_size(); ++i)
        m_source.at(i) = sr.source(i);

      m_target.resize(sr.target_size());
      for (int i=0; i < sr.target_size(); ++i)
        m_target.at(i) = sr.target(i);
    }

  std::ostream&
    serialise(std::ostream& out) const {
      SerialisedRule sr;
      to_message(sr);
      sr.SerializeToOstream(&out);

      return out;
    }

  std::istream&
    deserialise(std::istream& in) {
      SerialisedRule sr;
      sr.ParseFromIstream(&in);
      from_message(sr);

      return in;
    }

private:
//  size_t          m_hash, m_source_hash, m_target_hash;
  IndexedToken    m_label;  
  IndexedString   m_source;  
  IndexedString   m_target;  
};

struct RuleLt {
  bool operator()(const Rule& lhs, const Rule& rhs) const {
    if (lhs.label() < rhs.label()) 
      return true;
    if (lhs.label() == rhs.label()) {
      if (lhs.source() < rhs.source()) 
        return true;
      if (lhs.source() == rhs.source()) {
        if (lhs.target() < rhs.target()) 
          return true;
      }
    }
    return false;
  }
};

inline std::ostream& 
print_rule(std::ostream& out, const Rule& r) {
  out << "[" << r.label() << "]->[";

  int num_source_tokens = r.source().size();
  if (num_source_tokens > 0)
    out << r.source().at(0);
  for (int i=1; i<num_source_tokens; ++i)
    out << "," << r.source().at(i);

  out << "|||";

  int num_target_tokens = r.target().size();
  if (num_target_tokens > 0)
    out << r.target().at(0);
  for (int i=1; i<num_target_tokens; ++i)
    out << "," << r.target().at(i);

  out << "]";

  return out;
}

inline std::ostream&
operator<<(std::ostream& out, const Rule &r) {
  return print_rule(out, r);
}

inline bool 
operator<(const Rule& left, const Rule& right) { 
  return RuleLt()(left, right);
}

inline bool 
operator==(const Rule& left, const Rule& right) { 
  return left.label() == right.label() && left.source() == right.source() && left.target() == right.target();
}

typedef tree<Rule> RuleTree;

bool
read_tree(std::istream &in, RuleTree &result);

void
decomposition_tree_to_rule_tree(const DecompositionTree& decomposition_tree, 
    const Sentence& source_sentence, const Sentence& target_sentence, 
    RuleTree& rule_tree);

std::pair<IndexedString, IndexedString>
rule_tree_yield(const RuleTree& rule_tree);

std::pair<IndexedString, IndexedString>
_rule_tree_yield(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it);

Alignment
rule_tree_to_alignment(const RuleTree& rule_tree);

std::ostream&
extract_source_tree(const RuleTree &tree, std::ostream& out, VocabPtr vocab);

std::ostream&
extract_target_tree(const RuleTree &tree, std::ostream& out, VocabPtr vocab);

void
rule_tree_node_widths(const RuleTree& rule_tree, const RuleTree::iterator&
                      rule_tree_it, std::list< std::pair<int,int> >& widths);

struct YieldItem {
  YieldItem(RuleTree::iterator t_it, int ti)
    : target_it(t_it), target_index(ti) {}
  RuleTree::iterator target_it;
  int target_index;
};
//std::vector<int>
typedef std::vector<YieldItem> YieldSequence;

YieldSequence
rule_tree_yield_sequence(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it, int& target_index, bool compress=false);

void 
split_node_with_m1(RuleTree& rule_tree, RuleTree::iterator& node_it, Model1TablePtr forward_prior, Model1TablePtr backward, VocabPtr vocab);

double
log_prior(const IndexedString& s, const IndexedString& t, VocabPtr vocab, 
          Model1TablePtr forward_prior, Model1TablePtr backward_prior);

class RuleTreeTerminalsItem;
typedef boost::shared_ptr<RuleTreeTerminalsItem> RuleTreeTerminalsItemPtr;

struct RuleTreeTerminalsItem {
  RuleTreeTerminalsItemPtr next_source;
  RuleTreeTerminalsItemPtr next_target;
  RuleTreeTerminalsItemPtr previous_source;
  RuleTreeTerminalsItemPtr previous_target;
  RuleTree::iterator node;
  int source_index;
  int target_index;
};

struct RuleTreeTerminalsList {
  RuleTreeTerminalsList() {}
  RuleTreeTerminalsList(RuleTreeTerminalsItemPtr node)
    : source_head(node), target_head(node), source_tail(node), target_tail(node) {}

  RuleTreeTerminalsItemPtr source_head;
  RuleTreeTerminalsItemPtr target_head;
  RuleTreeTerminalsItemPtr source_tail;
  RuleTreeTerminalsItemPtr target_tail;
};

RuleTreeTerminalsList
terminals(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it);

#endif // _RULE_HH
