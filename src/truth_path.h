/*
  By: Ke @ PSU
  Last edited: 04/23/25
*/

#ifndef __TRUTH_PATH_H__
#define __TRUTH_PATH_H__

#include <unordered_set>
#include "splice_graph.h"

class truth_path
{
public:
    double weight;
    std::unordered_set<edge_descriptor> edges;

    truth_path(double w);
    void add_edge(edge_descriptor e);
    void remove_edge(edge_descriptor e);
    
    bool contains_edge(edge_descriptor e) const;
    // return true if edges contains each edge in s
    bool contains_subseq(const std::vector<edge_descriptor>& s) const;
};

#endif
