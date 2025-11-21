#include "truth_path.h"

truth_path::truth_path(double w) : weight(w)
{}

void truth_path::add_edge(edge_descriptor e)
{
    auto result = edges.insert(e);
    assert(result.second);
}

void truth_path::remove_edge(edge_descriptor e)
{
    auto result = edges.find(e);
    if(result != edges.end()) edges.erase(result);
}

bool truth_path::contains_edge(edge_descriptor e) const
{
    return edges.find(e) != edges.end();
}

bool truth_path::contains_subseq(const std::vector<edge_descriptor>& s) const
{
    for(edge_descriptor e : s)
    {
	if(!contains_edge(e)) return false;
    }
    return true;
}
