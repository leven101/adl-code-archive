#include "rule.hh"
#include "random.hh"

#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace boost;

  void
_decomposition_tree_to_rule_tree(const DecompositionTree::iterator& decomposition_tree_it, 
    const DecompositionTree& decomposition_tree, 
    const Sentence& source_sentence, const Sentence& target_sentence, 
    const RuleTree::iterator& rule_tree_it, RuleTree& rule_tree)
{
  if (decomposition_tree.number_of_children(decomposition_tree_it) == 0)
    return;

  for (DecompositionTree::sibling_iterator child_it=decomposition_tree.begin(decomposition_tree_it); 
      child_it != decomposition_tree.end(decomposition_tree_it); ++child_it) { 

    Rule new_child_rule_node(decomposition_tree, child_it, source_sentence, target_sentence);

    RuleTree::iterator new_child_rule_node_it;
    new_child_rule_node_it = rule_tree.append_child(rule_tree_it, new_child_rule_node);
    _decomposition_tree_to_rule_tree(child_it, decomposition_tree, source_sentence, target_sentence, new_child_rule_node_it, rule_tree);
  }

  return;
}

  void
decomposition_tree_to_rule_tree(const DecompositionTree& decomposition_tree, 
    const Sentence& source_sentence, const Sentence& target_sentence, 
    RuleTree& rule_tree)
{
  if (decomposition_tree.begin() == decomposition_tree.end())
    return;

  const DecompositionTree::iterator decomposition_tree_it=decomposition_tree.begin();
  Rule new_rule_node(decomposition_tree, decomposition_tree_it, source_sentence, target_sentence);
  RuleTree::iterator top_rule = rule_tree.set_head(new_rule_node);

  // recursive tree construction
  _decomposition_tree_to_rule_tree(decomposition_tree.begin(), decomposition_tree, source_sentence, target_sentence, top_rule, rule_tree);

  return;
}

  bool
read_tree(istream &in, RuleTree &result)
{
  string buf, token;
  if (!getline(in, buf)) 
    return false;

  TokenString tokens;
  istringstream ss(buf);

  TokenString surface_tokens;
  RuleTree::iterator current_parent = result.begin();

  while(ss >> token)
  {
    //cerr << "  Token: " << token<< endl;
    // new constituent
    if (token[0] == '(') {
      Rule new_rule(token.substr(1,token.size() - 1));
      //cerr << "New Phrase "; print_phrase(cerr, new_phrase) << endl;
      if (current_parent == result.end())
        current_parent = result.insert(current_parent, new_rule);
      else
        current_parent = result.append_child(current_parent, new_rule);
    }
    // leaf and constituent closing brackets
    else {
      // add the surface form
      string::size_type bracket = token.find(')');
      Rule new_leaf(token.substr(0,token.size() - (bracket == string::npos ? 0 : token.size() - bracket)));
      //cerr << "New Leaf "; print_phrase(cerr, new_leaf) << endl;
      if (current_parent == result.end())
        current_parent = result.insert(current_parent, new_leaf);
      else
        result.append_child(current_parent, new_leaf);

      // pop the constituents
      if (bracket != string::npos) {
        unsigned num_constituents = token.size() - bracket;
        for (unsigned i=0; i<num_constituents; i++) {
          assert(current_parent != result.end());
          //cerr << "Phrase: "; print_phrase(cerr, *current_parent);
          current_parent = result.parent(current_parent);
        }
      }
    }
    //cerr << "Tree "; kptree::print_tree_bracketed(result, std::cerr);
    //cerr << endl;
  }

  return true;
}

pair<IndexedString, IndexedString>
_rule_tree_yield(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it) {
  int number_of_children = rule_tree.number_of_children(rule_tree_it);

  // A terminal
  if (number_of_children == 0) {
    return make_pair(rule_tree_it->source(), rule_tree_it->target());
  }
  
  assert (number_of_children == static_cast<int>(rule_tree_it->source().size()));
  assert (number_of_children == static_cast<int>(rule_tree_it->target().size()));

  IndexedString target_yield;
  vector<IndexedString> reordered_source_yield(number_of_children);
  int counter=0;
  for (RuleTree::sibling_iterator child_it=rule_tree.begin(rule_tree_it); 
      child_it != rule_tree.end(rule_tree_it); ++child_it, ++counter) {
    // recurse to get the childs yield
    pair<IndexedString, IndexedString> child_yield = _rule_tree_yield(rule_tree, child_it);

    // place the source in the appropriate order
    reordered_source_yield.at(rule_tree_it->source().at(counter)) = child_yield.first; 

    // append to the target in linear order
    target_yield.insert(target_yield.end(), child_yield.second.begin(), child_yield.second.end());
  } 

  IndexedString source_yield;
  for (vector<IndexedString>::iterator source_it=reordered_source_yield.begin();
      source_it != reordered_source_yield.end(); ++source_it) {
    source_yield.insert(source_yield.end(), source_it->begin(), source_it->end());
  }

  return make_pair(source_yield, target_yield);
}

pair<IndexedString, IndexedString>
rule_tree_yield(const RuleTree& rule_tree) {
  return _rule_tree_yield(rule_tree, rule_tree.begin());
}

//vector<int>
YieldSequence
rule_tree_yield_sequence(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it, int& target_index, bool compress) {
  int number_of_children = rule_tree.number_of_children(rule_tree_it);

  // A terminal
  if (number_of_children == 0) {
    int number_of_indexes = (compress ? 1 : rule_tree_it->source().size());
    YieldSequence target_index_vector(number_of_indexes, YieldItem(rule_tree_it, target_index));
    //cerr << "Added: " << number_of_indexes << " " << target_index << endl;
    target_index++;
    return target_index_vector;
  }
  
  assert (number_of_children == static_cast<int>(rule_tree_it->source().size()));
  assert (number_of_children == static_cast<int>(rule_tree_it->target().size()));

  vector< YieldSequence > reordered_source_yield(number_of_children);
  int counter=0;
  for (RuleTree::sibling_iterator child_it=rule_tree.begin(rule_tree_it); 
      child_it != rule_tree.end(rule_tree_it); ++child_it, ++counter) {
    // recurse to get the childs yield
    YieldSequence child_yield = rule_tree_yield_sequence(rule_tree, child_it, target_index, compress);

    // place the source in the appropriate order
    reordered_source_yield.at(rule_tree_it->source().at(counter)) = child_yield; 
  } 

  YieldSequence source_yield;
  for (vector<YieldSequence>::iterator source_it=reordered_source_yield.begin();
      source_it != reordered_source_yield.end(); ++source_it) {
    source_yield.insert(source_yield.end(), source_it->begin(), source_it->end());
  }

  return source_yield;
}

// Collect a post-order sequence of the number of source/target words spanned by
// each tree node.
void
rule_tree_node_widths(const RuleTree& rule_tree, 
                      const RuleTree::iterator& rule_tree_it, 
                      list< pair<int,int> >& widths) {
  int number_of_children = rule_tree.number_of_children(rule_tree_it);
  int source_width=0, target_width=0;

  // A terminal
  if (number_of_children == 0) {
    source_width = rule_tree_it->source().size();
    target_width = rule_tree_it->target().size();
  }
  else {
    assert (number_of_children == static_cast<int>(rule_tree_it->source().size()));
    assert (number_of_children == static_cast<int>(rule_tree_it->target().size()));
    for (int c=number_of_children-1; c >= 0; --c) {
      RuleTree::iterator child_it = RuleTree::child(rule_tree_it, c);
      rule_tree_node_widths(rule_tree, child_it, widths);
      source_width += widths.front().first;
      target_width += widths.front().second;
    } 
  }
  //cerr << source_width << " " << target_width << " " << *rule_tree_it << endl;
  widths.push_front(make_pair(source_width, target_width));
} // rule_tree_node_widths

Alignment
rule_tree_to_alignment(const RuleTree& rule_tree) {
  int target_index=0;
  //vector<int> source_yield_sequence = rule_tree_yield_sequence(rule_tree, rule_tree.begin(), target_index);
  YieldSequence source_yield_sequence = rule_tree_yield_sequence(rule_tree, rule_tree.begin(), target_index);
  int source_index=0;
  multimap<int, int> target_leaf_to_source_index;
  for (YieldSequence::iterator source_seq_it=source_yield_sequence.begin(); 
      source_seq_it != source_yield_sequence.end(); ++source_seq_it, ++source_index) {
    target_leaf_to_source_index.insert(make_pair(source_seq_it->target_index, source_index)); 
  }

  Alignment alignment_result;
  int current_target_token_index=0;
  int current_target_segment_index=0;
  for (RuleTree::leaf_iterator leaf_it=rule_tree.begin_leaf(); leaf_it != rule_tree.end_leaf(); ++leaf_it, ++current_target_segment_index) {
    // for each source position aligned to this target rule
    multimap<int, int>::iterator source_it, source_end;
    tie(source_it, source_end) = target_leaf_to_source_index.equal_range(current_target_segment_index);

    for (; source_it != source_end; ++source_it) {
      // for each target token in this rule
      for (size_t token_index=0; token_index < leaf_it->target().size(); ++token_index) {
        alignment_result.insert(AlignmentPoint(source_it->second, current_target_token_index + token_index));
      }
    }
    current_target_token_index += leaf_it->target().size();
  }

  return alignment_result;
}

void 
split_node_with_m1(RuleTree& rule_tree, RuleTree::iterator& node_it, Model1TablePtr forward_prior, Model1TablePtr backward_prior, VocabPtr vocab) {
  bool debug=false;
  if (debug) 
    cerr << "Split Node with M1" << endl;

  int num_children = rule_tree.number_of_children(node_it);
  const IndexedString& rule_source = node_it->source();
  const IndexedString& rule_target = node_it->target();

  if (num_children != 0 || rule_source.size() <= 1 || rule_target.size() <= 1) {
    cerr << "split_node_with_m1(it="; print_rule(cerr, *node_it) << "): can't split rule." << endl;
    return;
  }

  boost::tuple< pair<RulePtr, RulePtr>, double, bool > max_split_option;
  double max_split_prob= -numeric_limits<double>::infinity();
  bool uninitialised=true;

  // create the source and target indexes for the proposed parent rules (both monotone and reordering)
  IndexedString new_source(2), new_target(2);
  new_target.at(0) = -1; new_target.at(1) = -1;
  new_source.at(0) = 0; new_source.at(1) = 1;
  Rule proposed_monotone_rule(-1, new_source, new_target);
  new_source.at(0) = 1; new_source.at(1) = 0;
  Rule proposed_reordered_rule(-1, new_source, new_target);

  for (size_t source_split=0; source_split <= rule_source.size(); ++source_split) {
    for (size_t target_split=0; target_split <= rule_target.size(); ++target_split) {
      //        if (debug) 
      //          cerr << "Split at (" << source_split << "," << target_split << ") ";

      // create the proposed child rules
      pair<RulePtr, RulePtr> monotone_rules, reordered_rules;
      boost::tie(monotone_rules, reordered_rules) = node_it->split(source_split, target_split, true, true);

    
      if (debug) {
        cerr << "Monotone proposed: ";
        print_rule(cerr, proposed_monotone_rule) << " -> ";
        print_rule(cerr, *(monotone_rules.first));
        cerr << "--";
        print_rule(cerr, *(monotone_rules.second));
        cerr << " Reordered child proposed: ";
        print_rule(cerr, proposed_reordered_rule) << " -> ";
        print_rule(cerr, *(reordered_rules.first));
        cerr << "--";
        print_rule(cerr, *(reordered_rules.second));
        cerr << endl;
      }

      // probability of generating the monotone configuration
      double proposed_monotone_log_prob=0;
      if (!(monotone_rules.first->is_half_null() || monotone_rules.second->is_half_null())) {
        proposed_monotone_log_prob = log_prior(monotone_rules.first->source(), monotone_rules.first->target(), vocab, forward_prior, backward_prior);
        proposed_monotone_log_prob += log_prior(monotone_rules.first->target(), monotone_rules.first->source(), vocab, forward_prior, backward_prior);

        proposed_monotone_log_prob += log_prior(monotone_rules.second->source(), monotone_rules.second->target(), vocab, forward_prior, backward_prior);
        proposed_monotone_log_prob += log_prior(monotone_rules.second->target(), monotone_rules.second->source(), vocab, forward_prior, backward_prior);

        proposed_monotone_log_prob *= 0.5;

        if (uninitialised || proposed_monotone_log_prob > max_split_prob) {
          max_split_option = boost::make_tuple(monotone_rules, proposed_monotone_log_prob, false);
          max_split_prob = proposed_monotone_log_prob;
          uninitialised = false;
        }
      }

      // probability of generating the reordering configuration
      double proposed_reordered_log_prob=0;
      if (!(reordered_rules.first->is_half_null() || reordered_rules.second->is_half_null())) {
        proposed_reordered_log_prob = log_prior(reordered_rules.first->source(), reordered_rules.first->target(), vocab, forward_prior, backward_prior);
        proposed_reordered_log_prob += log_prior(reordered_rules.first->target(), reordered_rules.first->source(), vocab, forward_prior, backward_prior);

        proposed_reordered_log_prob += log_prior(reordered_rules.second->source(), reordered_rules.second->target(), vocab, forward_prior, backward_prior);
        proposed_reordered_log_prob += log_prior(reordered_rules.second->target(), reordered_rules.second->source(), vocab, forward_prior, backward_prior);

        proposed_reordered_log_prob *= 0.5;

        if (uninitialised || proposed_reordered_log_prob > max_split_prob) {
          max_split_option = boost::make_tuple(reordered_rules, proposed_reordered_log_prob, true);
          max_split_prob = proposed_reordered_log_prob;
          uninitialised = false;
        }
      }

      if (debug) 
        cerr << "  Prop. Mono. LogProb=" << proposed_monotone_log_prob << " Pro. Reord. LogProb=" << proposed_reordered_log_prob << endl;
    }
  }

  if (debug) 
    cerr << "  Best LogProb==" << max_split_prob << " Best Pair==" << *(max_split_option.get<0>().first) << " ||| " << *(max_split_option.get<0>().second) << endl;

  assert (!uninitialised);

  if (max_split_option.get<2>()) {
    const pair<RulePtr, RulePtr>& reordered_rules = max_split_option.get<0>();
    RuleTree::sibling_iterator new_right_it = rule_tree.append_child(node_it, *(reordered_rules.second));
    rule_tree.insert(new_right_it, *(reordered_rules.first));
    *node_it = proposed_reordered_rule;
  }
  else {
    const pair<RulePtr, RulePtr>& monotone_rules = max_split_option.get<0>();
    RuleTree::sibling_iterator new_right_it = rule_tree.append_child(node_it, *(monotone_rules.second));
    rule_tree.insert(new_right_it, *(monotone_rules.first));
    *node_it = proposed_monotone_rule;
  }
}

RuleTreeTerminalsList
_terminals(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it, int& target_index) {
  int number_of_children = rule_tree.number_of_children(rule_tree_it);

  // A terminal
  if (number_of_children == 0) {
    RuleTreeTerminalsItemPtr terminal_node(new RuleTreeTerminalsItem);
    terminal_node->node = rule_tree_it;
    terminal_node->target_index = target_index++;
    return RuleTreeTerminalsList(terminal_node);
  }
  
  assert (number_of_children == static_cast<int>(rule_tree_it->source().size()));
  assert (number_of_children == static_cast<int>(rule_tree_it->target().size()));

  RuleTreeTerminalsList previous_child, result;
  vector< RuleTreeTerminalsList > reordered_source_yield(number_of_children);

  int counter=0;
  for (RuleTree::sibling_iterator child_it=rule_tree.begin(rule_tree_it); 
      child_it != rule_tree.end(rule_tree_it); ++child_it, ++counter) {
    // recurse to get the childs yield
    RuleTreeTerminalsList child_yield = _terminals(rule_tree, child_it, target_index);
    child_yield.target_head->previous_target = previous_child.target_tail;

    if (result.target_head.get() == 0)
      result.target_head = child_yield.target_head;

    // place the source in the appropriate order
    reordered_source_yield.at(rule_tree_it->source().at(counter)) = child_yield; 
  } 
  result.target_tail = previous_child.target_tail;

  previous_child = RuleTreeTerminalsList();
  for (vector< RuleTreeTerminalsList>::iterator source_it=reordered_source_yield.begin();
      source_it != reordered_source_yield.end(); ++source_it) {
    if (result.source_head.get() == 0)
      result.source_head = source_it->source_head;
    
    source_it->source_head->previous_source = previous_child.source_tail;
  }
  result.source_tail = previous_child.source_tail;

  return result;
}

RuleTreeTerminalsList
terminals(const RuleTree& rule_tree, const RuleTree::iterator& rule_tree_it) {
  int num_leaves=0;
  RuleTreeTerminalsList leaves = _terminals(rule_tree, rule_tree_it, num_leaves); 

  int source_counter=0;
  for (RuleTreeTerminalsItemPtr leaf=leaves.source_head; leaf.get() != 0; leaf = leaf->next_source, ++source_counter)
    leaf->source_index = source_counter;

  assert(source_counter == num_leaves);
  return leaves;
}

double
log_prior(const IndexedString& s, const IndexedString& t, VocabPtr vocab, 
          Model1TablePtr forward_prior, Model1TablePtr backward_prior) {
  static Poisson poisson;

  // nulls have a simple language model prior on their non-null side
  if (s.empty() || t.empty()) {
    int non_null_length = std::max(s.size(), t.size());
    double lm_log_prob = poisson.log_prob(non_null_length, 1) + log(0.0001)*(non_null_length);

    return (log(0.5) + lm_log_prob);
  }

  double lm_log_prob_target = poisson.log_prob(t.size(), 1);// + log(0.00001)*(t.size());
  double lm_log_prob_source = poisson.log_prob(s.size(), 1);// + log(0.00001)*(s.size());
  double log_prob_target_length = poisson.log_prob(t.size(), s.size());// + log(0.00001)*(t.size());
  double log_prob_source_length = poisson.log_prob(s.size(), t.size());// + log(0.00001)*(s.size());

  if (backward_prior.get() == 0 || forward_prior.get() == 0) {
    return 0.5 * (log_prob_source_length + lm_log_prob_target + log_prob_target_length + lm_log_prob_source);
  }

  // Model 1 prior
  double log_prob_backward = backward_prior->log_prob(s, t, vocab) + lm_log_prob_target + log_prob_source_length;
  double log_prob_forward = forward_prior->log_prob(t, s, vocab) + lm_log_prob_source + log_prob_target_length;

  return 0.5 * (log_prob_backward + log_prob_forward);
}

ostream&
extract_subtree(const RuleTree &t, RuleTree::iterator iRoot, ostream& str, VocabPtr vocab, bool source)
{
  if(t.empty()) 
    return str;

  int number_of_children = t.number_of_children(iRoot);
  if (number_of_children == 0) {
    const TokenString leaves = vocab->lookup(source ? iRoot->source() : iRoot->target());
    if (leaves.size() == 1)
      str << " (" << iRoot->label() << "- " << map_brakets(leaves.at(0)) << ")";
    else if (leaves.size() > 1) {
      str << " (" << iRoot->label() << "-";
      for (TokenString::const_iterator it=leaves.begin(); it != leaves.end(); ++it)
        str << " (-L- " << map_brakets(*it) << ")";
      str << ")";
    }
  }
  else {
    assert(number_of_children == static_cast<int>(iRoot->source().size()));
    assert(number_of_children == static_cast<int>(iRoot->target().size()));

    // parent
    str << " (" << iRoot->label() << "-";

    // get the ordering
    vector<RuleTree::sibling_iterator> ordered_children(number_of_children);

    int siblingNum;
    RuleTree::sibling_iterator iChildren;
    for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren, ++siblingNum) {
      if (source)
        ordered_children.at(iRoot->source().at(siblingNum)) = iChildren;
      else
        ordered_children.at(siblingNum) = iChildren;
    }

    // print the children in source order
    for (int i=0; i < static_cast<int>(ordered_children.size()); ++i)
      extract_subtree(t, ordered_children.at(i), str, vocab, source);

    str << ")";
  }

  return str;
}

ostream&
extract_tree(const RuleTree &tree, ostream& str, VocabPtr vocab, bool source)
{
  int headCount = tree.number_of_siblings(tree.begin());
  int headNum = 0;
  for(RuleTree::sibling_iterator iRoots = tree.begin(); iRoots != tree.end(); ++iRoots, ++headNum) {
    extract_subtree(tree, iRoots, str, vocab, source);
    if (headNum != headCount) {
      str << std::endl;
    }
  }

  return str;
}

ostream&
extract_source_tree(const RuleTree &tree, ostream& str, VocabPtr vocab) {
  return extract_tree(tree, str, vocab, true);
}

ostream&
extract_target_tree(const RuleTree &tree, ostream& str, VocabPtr vocab) {
  return extract_tree(tree, str, vocab, false);
}
