/*
  Part of Scallop Transcript Assembler
  (c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
  See LICENSE for licensing.
*/

#ifndef __VERTEX_BASE_H__
#define __VERTEX_BASE_H__

#include <set>
#include "edge_base.h"

using namespace std;

class vertex_base
{
public:
    vertex_base();
    virtual ~vertex_base();

protected:
    struct edge_pointer_cmp
    {
	bool operator()(const edge_base* a, const edge_base* b) const
	{
	    if((a->source()) == (b->source()))
	    {
		if((a->target()) == (b->target()))
		{
		    return a < b;
		}
		else
		{
		    return (a->target()) < (b->target());
		}
	    }
	    else
	    {
		return (a->source()) < (b->source());
	    }
	}
    };
    set<edge_base*, edge_pointer_cmp> si;		// in_edges
    set<edge_base*, edge_pointer_cmp> so;		// out_edges

public:
    virtual int add_in_edge(edge_base *e);
    virtual int add_out_edge(edge_base *e);
    virtual int remove_in_edge(edge_base *e);
    virtual int remove_out_edge(edge_base *e);
    virtual int degree() const;
    virtual int in_degree() const;
    virtual int out_degree() const;
    virtual PEEI in_edges() const;
    virtual PEEI out_edges() const;
    virtual int print() const;
};

#endif
