/*
  Part of Catfish
  (c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
  (c) 2025 by Ke Chen @ Penn State.
  See LICENSE for licensing.
*/

#include "catfish.h"
#include "subsetsum.h"
#include "config.h"
#include "nested_graph.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <fstream>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

catfish::catfish()
{}

catfish::catfish(const string &s, const splice_graph &g)
    : name(s), gr(g)
{
    round = 0;
    has_truth = false;
    //assert(gr.check_fully_connected() == true);
    try
    {
	gurobi_env = new GRBEnv();
	gurobi_base_model_built = false;
    }
    catch(GRBException& e)
    {
	std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	std::cerr << e.getMessage() << std::endl;
    }
}

catfish::catfish(const string &name, const splice_graph &gr,
		 const std::vector<std::pair<double, std::vector<int>>>& truth)
    :catfish(name, gr)
{
    has_truth = true;
    build_truth(truth);
}

catfish::~catfish()
{
    delete gurobi_env;
}

int catfish::clear()
{
    name = "";
    gr.clear();
    e2i.clear();
    i2e.clear();
    mev.clear();
    round = -1;
    paths.clear();
    return 0;
}

int catfish::save(catfish &sc)
{
    sc.clear();

    sc.name = name;
    MEE x2y, y2x;
    sc.gr.copy(gr, x2y, y2x);

    for(MEI::iterator it = e2i.begin(); it != e2i.end(); it++)
    {
	edge_descriptor e = it->first;
	int k = it->second;
	assert(k >= 0 && k < i2e.size());
	if(i2e[k] == null_edge)
	{
	    assert(x2y.find(e) == x2y.end());
	    sc.e2i.insert(PEI(e, k));
	}
	else
	{
	    assert(x2y.find(e) != x2y.end());
	    sc.e2i.insert(PEI(x2y[e], k));
	}
    }

    for(int i = 0; i < i2e.size(); i++)
    {
	if(i2e[i] == null_edge)
	{
	    sc.i2e.push_back(null_edge);
	}
	else
	{
	    assert(x2y.find(i2e[i]) != x2y.end());
	    sc.i2e.push_back(x2y[i2e[i]]);
	}
    }

    for(MEV::iterator it = mev.begin(); it != mev.end(); it++)
    {
	edge_descriptor e = it->first;
	vector<int> v = it->second;
	if(x2y.find(e) == x2y.end())
	{
	    if(e2i.find(e) == e2i.end()) continue;
	    int k = e2i[e];
	    assert(i2e[k] == null_edge);
	    assert(sc.i2e[k] == null_edge);
	    sc.mev.insert(PEV(e, v));
	}
	else
	{
	    int k = e2i[e];
	    assert(i2e[k] == e);
	    assert(sc.e2i[x2y[e]] == k);
	    assert(sc.i2e[k] == x2y[e]);
	    sc.mev.insert(PEV(x2y[e], v));
	}
    }

    sc.i2backtrack = i2backtrack; // auto deep copy
    sc.round = round;
    sc.paths = paths;

    return 0;
}

int catfish::load(catfish &sc)
{
    clear();

    name = sc.name;

    MEE x2y, y2x;
    gr.copy(sc.gr, x2y, y2x);

    for(MEI::iterator it = sc.e2i.begin(); it != sc.e2i.end(); it++)
    {
	edge_descriptor e = it->first;
	int k = it->second;
	assert(k >= 0 && k < sc.i2e.size());
	if(sc.i2e[k] == null_edge)
	{
	    assert(x2y.find(e) == x2y.end());
	    e2i.insert(PEI(e, k));
	}
	else
	{
	    assert(x2y.find(e) != x2y.end());
	    e2i.insert(PEI(x2y[e], k));
	}
    }

    for(int i = 0; i < sc.i2e.size(); i++)
    {
	if(sc.i2e[i] == null_edge)
	{
	    i2e.push_back(null_edge);
	}
	else
	{
	    assert(x2y.find(sc.i2e[i]) != x2y.end());
	    i2e.push_back(x2y[sc.i2e[i]]);
	}
    }

    for(MEV::iterator it = sc.mev.begin(); it != sc.mev.end(); it++)
    {
	edge_descriptor e = it->first;
	vector<int> v = it->second;
	if(x2y.find(e) == x2y.end())
	{
	    if(sc.e2i.find(e) == sc.e2i.end()) continue;
	    int k = sc.e2i[e];
	    assert(sc.i2e[k] == null_edge);
	    assert(i2e[k] == null_edge);
	    mev.insert(PEV(e, v));
	}
	else
	{
	    int k = sc.e2i[e];
	    assert(sc.i2e[k] == e);
	    assert(e2i[x2y[e]] == k);
	    assert(i2e[k] == x2y[e]);
	    mev.insert(PEV(x2y[e], v));
	}
    }

    round = sc.round;
    paths = sc.paths;
    i2backtrack = sc.i2backtrack;

    return 0;
}

int catfish::assemble()
{
    int c = classify();
    //if(c == TRIVIAL) return 0;

    if(algo == "core") return assemble1();
    if(algo == "full") return assemble2();
    if(algo == "greedy") return greedy();
    return 0;
}

int catfish::classify()
{
    assert(gr.num_vertices() >= 2);
    if(gr.num_vertices() == 2) return TRIVIAL;

    string s;	

    long p0 = gr.compute_num_paths();
    long p1 = gr.num_edges() - gr.num_vertices() + 2;
    for(int i = 0; i < gr.num_vertices(); i++) 
    {
	if(gr.degree(i) == 0) p1++;
    }

    printf("vertices = %lu, edges = %lu, p0 = %ld, p1 = %ld\n", gr.num_vertices(), gr.num_edges(), p0, p1);

    assert(p0 >= p1);

    bool b = (p0 == p1) ? true : false;

    printf("\nprocess %s %s\n", name.c_str(), b ? "TRIVIAL" : "NORMAL");

    if(p0 == p1) return TRIVIAL;
    else return NORMAL;
}

int catfish::assemble0()
{
    if(output_tex_files == true) gr.draw(name + "." + tostring(round) + ".tex");
    round += 1;

    //gr.round_weights();
    //remove_empty_edges();

    //if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

    init_super_edges();
    reconstruct_splice_graph();
    gr.get_edge_indices(i2e, e2i);

    i2backtrack = std::vector<std::vector<edge_descriptor>>(i2e.size());
    for(int i = 0; i < i2e.size(); ++i)
    {
	edge_descriptor e = i2e[i];
	i2backtrack[i].emplace_back(e);
	original_e2st.emplace(e, std::make_pair(e->source(), e->target()));
    }
								  

    //print();
    //collect_existing_st_paths();
    //printf("%s catfish0 solution %lu paths\n", name.c_str(), paths.size());

    return 0;
}

int catfish::assemble1()
{
    assemble0();

    int f = iterate();

    collect_existing_st_paths();
    //print();

    printf("%s core solution %lu paths, iteration = %d\n\n", name.c_str(), paths.size(), f);

    return 0;
}

int catfish::assemble2()
{
    assemble0();

    int f = iterate();

    collect_existing_st_paths();
    size_t paths_before_greedy = paths.size();
    printf("%s before greedy, decomposed %zu paths\n", name.c_str(), paths_before_greedy);
    print();

    greedy_decompose();
    assert(gr.num_edges() == 0);

    //print();

    printf("%s full solution %lu paths include %zu greedy paths\n\n",
	   name.c_str(), paths.size(), paths.size() - paths_before_greedy);

    return 0;
}

int catfish::greedy()
{
    assemble0();

    greedy_decompose();
    assert(gr.num_edges() == 0);

    //print();

    printf("%s greedy solution %lu paths\n\n", name.c_str(), paths.size());

    return 0;
}

int catfish::iterate()
{   
    print();
    while(true)
    {
	bool b = false;

	b = decompose_trivial_vertices();
	if(b)
	{
	    print();
	    continue;
	}

	b = decompose_adjacent_equal_edges();
	if(b)
	{
	    print();
	    continue;
	}

	b = decompose_with_equations();
	if(b)
	{
	    print();
	    continue;
	}

	if(gurobi_base_model_built)
	{
	    try
	    {
		// last equation test may have caused the model to be infeasible
		if(gurobi_model->get(GRB_IntAttr_Status) != GRB_OPTIMAL)
		{
		    gurobi_is_feasible();
		}
		b = decompose_one_path_from_LP();
		if(b)
		{
		    print();
		    continue;
		}
	    }
	    catch(GRBException& e)
	    {
		std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
		std::cerr << e.getMessage() << std::endl;
		throw;
	    }
	}
	
	break;
    }

    gurobi_clear_model();
    return 0;
}

bool catfish::decompose_with_equations(int level)
{
    bool b = false;
    if(level == 0)
    {
	if(b == false) b = identify_equations0();
    }
    else if(level == 1)
    {
	if(b == false) b = identify_equations1();
	if(b == false) b = identify_equations2();
    }
    else if(level == 2)
    {
	if(b == false) b = identify_equations3();
    }

    if(eqns.size() == 0) return false;

    sort(eqns.begin(), eqns.end(), equation_cmp2);

    printf("equations = %lu\n", eqns.size());

    for(int i = 0; i < eqns.size(); i++) 
    {
	printf(" ");
        // eqns[i].print(i);
	print_equation(eqns[i], i);
    }

    if(fast_mode == false && (eqns[0].f < 2 || eqns[0].d > 0))
    {
	sort(eqns.begin(), eqns.end(), equation_cmp1);
	resolve_vertex_with_equations(eqns);
	sort(eqns.begin(), eqns.end(), equation_cmp2);
    }


    if(eqns[0].f == 3) return true;

    if(eqns[0].f == 2 && eqns[0].d == 0)
    {
	resolve_equation(eqns[0]);
	return true;
    }

    return false;
}

int catfish::init_super_edges()
{
    mev.clear();
    edge_iterator it1, it2;
    for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
    {
	vector<int> v;
	int s = (*it1)->source();
	v.push_back(s);
	mev.insert(PEV(*it1, v));
    }
    return 0;
}

int catfish::reconstruct_splice_graph()
{
    while(true)
    {
	bool flag = false;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
	    bool b = init_trivial_vertex(i);
	    if(b == true) flag = true;
	}
	if(flag == false) break;
    }
    return 0;
}

bool catfish::init_trivial_vertex(int x)
{
    int id = gr.in_degree(x);
    int od = gr.out_degree(x);

    if(id <= 0 || od <= 0) return false;
    if(id >= 2 && od >= 2) return false;
    //if(id <= 1 && od <= 1) return false;

    edge_iterator it1, it2;
    edge_iterator ot1, ot2;
    for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
    {

	for(tie(ot1, ot2) = gr.out_edges(x); ot1 != ot2; ot1++)
	{
	    int s = (*it1)->source();
	    int t = (*ot1)->target();

	    double w1 = gr.get_edge_weight(*it1);
	    double w2 = gr.get_edge_weight(*ot1);

	    double w = w1 < w2 ? w1 : w2;
	    
	    edge_descriptor p = gr.add_edge(s, t);
	    gr.set_edge_weight(p, w);

	    // adjust grount truth paths
	    if(has_truth)
	    {
		std::vector<edge_descriptor> both(2);
		both[0] = *it1;
		both[1] = *ot1;
		for(truth_path& tp : truth)
		{
		    if(tp.contains_subseq(both))
		    {
			for(edge_descriptor e : both)
			{
			    tp.remove_edge(e);
			}
			tp.add_edge(p);
		    }
		}
	    }

	    assert(mev.find(*it1) != mev.end());
	    assert(mev.find(*ot1) != mev.end());

	    vector<int> v1 = mev[*it1];
	    vector<int> v2 = mev[*ot1];
	    v1.insert(v1.end(), v2.begin(), v2.end());

	    if(mev.find(p) != mev.end()) mev[p] = v1;
	    else mev.insert(PEV(p, v1));
	}
    }
    gr.clear_vertex(x);
    return true;
}

int catfish::split_merge_path(const VE &p, double wx, vector<int> &vv)
{
    vector<int> v;
    for(int i = 0; i < p.size(); i++)
    {
	assert(p[i] != null_edge);
	assert(e2i.find(p[i]) != e2i.end());
	v.push_back(e2i[p[i]]);
    }
    return split_merge_path(v, wx, vv);
}

int catfish::split_merge_path(const vector<int> &p, double ww, vector<int> &vv)
{
    vv.clear();
    if(p.size() == 0) return -1;
    int ee = p[0];
    int x = split_edge(p[0], ww);
    vv.push_back(x);
    for(int i = 1; i < p.size(); i++)
    {
	x = split_edge(p[i], ww);
	vv.push_back(x);
	ee = merge_adjacent_equal_edges(ee, p[i]);
    }
    return ee;
}

int catfish::merge_adjacent_equal_edges(int x, int y)
{
    if(i2e[x] == null_edge) return -1;
    if(i2e[y] == null_edge) return -1;

    edge_descriptor xx = i2e[x];
    edge_descriptor yy = i2e[y];

    int xs = (xx)->source();
    int xt = (xx)->target();
    int ys = (yy)->source();
    int yt = (yy)->target();

    if(xt != ys && yt != xs) return -1;
    if(yt == xs) return merge_adjacent_equal_edges(y, x);
	
    assert(xt == ys);

    edge_descriptor p = gr.add_edge(xs, yt);

    int n = i2e.size();
    i2e.push_back(p);
    assert(e2i.find(p) == e2i.end());
    e2i.insert(PEI(p, n));

    double wx0 = gr.get_edge_weight(xx);
    double wy0 = gr.get_edge_weight(yy);

    assert(fabs(wx0 - wy0) <= SMIN);

    gr.set_edge_weight(p, wx0);

    vector<int> v = mev[xx];
    v.insert(v.end(), mev[yy].begin(), mev[yy].end());

    if(mev.find(p) != mev.end()) mev[p] = v;
    else mev.insert(PEV(p, v));

    std::vector<edge_descriptor> vb = i2backtrack[x];
    vb.insert(vb.end(), i2backtrack[y].begin(), i2backtrack[y].end());
    i2backtrack.push_back(vb);

    assert(i2e[n] == p);
    assert(e2i.find(p) != e2i.end());
    assert(e2i[p] == n);
    assert(e2i[i2e[n]] == n);

    e2i.erase(xx);
    e2i.erase(yy);
    i2e[x] = null_edge;
    i2e[y] = null_edge;
    gr.remove_edge(xx);
    gr.remove_edge(yy);

    return n;
}

int catfish::merge_adjacent_edges(int x, int y)
{
    if(i2e[x] == null_edge) return -1;
    if(i2e[y] == null_edge) return -1;

    if(has_truth && equation_comply_with_truth(
	   equation(std::vector<int>{x}, std::vector<int>{y})))
    {
	printf("    merge is OK with ground truth\n");
    }

    edge_descriptor xx = i2e[x];
    edge_descriptor yy = i2e[y];

    double wx = gr.get_edge_weight(xx);
    double wy = gr.get_edge_weight(yy);
    double ww = (wx <= wy) ? wx : wy;

    split_edge(x, ww);
    split_edge(y, ww);

    return merge_adjacent_equal_edges(x, y);
}

int catfish::split_edge(int ei, double w)
{
    assert(i2e[ei] != null_edge);
    edge_descriptor ee = i2e[ei];

    double ww = gr.get_edge_weight(ee);

    if(fabs(ww - w) <= SMIN) return ei;

    assert(ww >= w + SMIN);

    int s = ee->source();
    int t = ee->target();

    edge_descriptor p2 = gr.add_edge(s, t);

    gr.set_edge_weight(ee, w);
    gr.set_edge_weight(p2, ww - w);

    if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
    else mev.insert(PEV(p2, mev[ee]));

    i2backtrack.push_back(i2backtrack[ei]);

    int n = i2e.size();
    i2e.push_back(p2);
    e2i.insert(PEI(p2, n));

    return n;
}

bool catfish::verify_equation_mergable(equation &eqn)
{
    if(eqn.s.size() == 0) return false; 
    if(eqn.t.size() == 0) return false;

    SE xx, yy;
    for(int i = 0; i < eqn.s.size(); i++)
    {
	edge_descriptor x = i2e[eqn.s[i]];
	SE se;
	gr.bfs(x->target(), se);
	xx.insert(se.begin(), se.end());
	se.clear();
	gr.bfs_reverse(x->source(), se);
	xx.insert(se.begin(), se.end());
    }
    for(int i = 0; i < eqn.t.size(); i++)
    {
	edge_descriptor x = i2e[eqn.t[i]];
	SE se;
	gr.bfs(x->target(), se);
	yy.insert(se.begin(), se.end());
	se.clear();
	gr.bfs_reverse(x->source(), se);
	yy.insert(se.begin(), se.end());
    }

    for(int i = 0; i < eqn.s.size(); i++)
    {
	edge_descriptor x = i2e[eqn.s[i]];
	if(yy.find(x) == yy.end()) return false;
    }

    for(int i = 0; i < eqn.t.size(); i++)
    {
	edge_descriptor x = i2e[eqn.t[i]];
	if(xx.find(x) == xx.end()) return false;
    }
    return true;
}

bool catfish::verify_equation_nontrivial(equation &eqn)
{
    set<edge_descriptor> fbs;
    set<edge_descriptor> fbt;
    vector<int> vsx;
    vector<int> vsy;
    vector<int> vtx;
    vector<int> vty;
    for(int i = 0; i < eqn.s.size(); i++)
    {
	edge_descriptor &e = i2e[eqn.s[i]];
	assert(e != null_edge);
	fbs.insert(e);
	vsx.push_back(e->source());
	vsy.push_back(e->target());
    }
    for(int i = 0; i < eqn.t.size(); i++)
    {
	edge_descriptor &e = i2e[eqn.t[i]];
	assert(e != null_edge);
	fbt.insert(e);
	vtx.push_back(e->source());
	vty.push_back(e->target());
    }

    int n = gr.num_vertices() - 1;
    bool b1 = gr.bfs(vsy, n, fbt);
    bool b2 = gr.bfs_reverse(vtx, 0, fbs);
    if(b1 == false && b2 == false) return false;

    b1 = gr.bfs(vty, n, fbs);
    b2 = gr.bfs_reverse(vsx, 0, fbt);
    if(b1 == false && b2 == false) return false;

    return true;
}

// adjacent equal edge
bool catfish::identify_equations0()
{
    for(int i = 0; i < gr.num_vertices(); i++)
    {
	if(gr.degree(i) == 0) continue;
	if(gr.in_degree(i) == 1) continue;
	if(gr.out_degree(i) == 1) continue;

	MI m1;
	MI m2;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
	{
	    int w = (int)(gr.get_edge_weight(*it1));
	    int e = e2i[*it1];
	    if(m1.find(w) != m1.end()) continue;
	    m1.insert(PI(w, e));
	}

	for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
	{
	    int w = (int)(gr.get_edge_weight(*it1));
	    int e = e2i[*it1];
	    if(m2.find(w) != m2.end()) continue;
	    m2.insert(PI(w, e));
	}

	for(MI::iterator it = m1.begin(); it != m1.end(); it++)
	{
	    int w = it->first;
	    if(m2.find(w) == m2.end()) continue;
	    vector<int> s;
	    vector<int> t;
	    s.push_back(it->second);
	    t.push_back(m2[w]);

	    equation eqn(s, t, 0);
	    eqn.f = 2;
	    eqn.a = 1;
	    eqn.d = 0;
	    eqn.w = w;

	    eqns.push_back(eqn);
	    return true;
	}
    }
    return false;
}

// distant equal edges
bool catfish::identify_equations1()
{
    typedef map<int, vector<int> > MIV;
    typedef pair<int, vector<int> > PIV;

    MIV m;
    edge_iterator it1, it2;
    for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
    {
	int w = (int)(gr.get_edge_weight(*it1));
	int e = e2i[*it1];
	if(m.find(w) != m.end()) 
	{
	    m[w].push_back(e);
	}
	else
	{
	    vector<int> v;
	    v.push_back(e);
	    m.insert(PIV(w, v));
	}
    }

    nested_graph nt(gr);
    for(MIV::iterator it = m.begin(); it != m.end(); it++)
    {
	vector<int> &v = it->second;
	if(v.size() <= 1) continue;
	for(int i = 0; i < v.size(); i++)
	{
	    for(int j = i + 1; j < v.size(); j++)
	    {
		bool b = check_adjacent_mergable(v[i], v[j], nt);
		int l = check_distant_mergable(v[i], v[j], it->first);

		vector<int> s;
		vector<int> t;
		s.push_back(v[i]);
		t.push_back(v[j]);

		equation eqn(s, t, 0);

		if(b == true) eqn.f = 2; // can be resolved by reverse operations
		else if(l >= 0) eqn.f = 1; // can be connected by a path with enough weight
		else eqn.f = 0; // not connected

		if(b == true)
		{
		    eqn.a = 1;
		    eqn.d = 0;
		}
		else
		{
		    eqn.a = 0;
		    eqn.d = 1;
		}

		eqn.w = it->first;

		if(b == true)
		{
		    eqns.push_back(eqn);
		    return true;
		}

		if(fast_mode == true)
		{
		    resolve_vertex_with_equation(eqn);
		    if(eqn.f == 3)
		    {
			eqns.push_back(eqn);
			return true;
		    }
		}

		// retain distant equal edges that can be connected by a path with enough weight
		if(l >= 0)
		{
		    eqns.push_back(eqn);
		}
	    }
	}
    }
    return false;
}

bool catfish::identify_equations2()
{
    vector<PI> p;
    for(int i = 0; i < i2e.size(); i++)
    {
	edge_descriptor &e = i2e[i];
	if(e == null_edge) continue;
	if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
	int w = (int)(gr.get_edge_weight(e));
	if(w <= 0) continue;
	p.push_back(PI(w, i));
    }

    if(p.size() == 0) return false;

    subsetsum sss1(p, p);
    sss1.solve();

    nested_graph nt(gr);
    for(int i = 0; i < sss1.eqns.size(); i++)
    {
	equation &eqn = sss1.eqns[i];
	assert(eqn.s.size() == 1);
	if(eqn.t.size() == 1) continue;
	if(verify_equation_mergable(eqn) == false) continue;
	if(verify_equation_nontrivial(eqn) == false) continue;

	if(fast_mode == true)
	{
	    for(int k = 0; k < eqn.t.size(); k++)
	    {
		bool b = check_adjacent_mergable(eqn.s[0], eqn.t[k], nt);
		if(b == false) eqn.d++;
		else eqn.a++;
	    }
	    resolve_vertex_with_equation(eqn);
	    eqns.push_back(eqn);
	    if(eqn.f == 3) return true;
	}
	else
	{
	    catfish sc;
	    save(sc);
	    resolve_equation(eqn);
	    load(sc);
	    eqns.push_back(eqn);
	    if(eqn.f == 2 && eqn.d == 0) return true;
	}
    }
    return false;
}

bool catfish::identify_equations3()
{
    vector<PI> p;
    for(int i = 0; i < i2e.size(); i++)
    {
	edge_descriptor &e = i2e[i];
	if(e == null_edge) continue;
	if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
	int w = (int)(gr.get_edge_weight(e));
	if(w <= 0) continue;
	p.push_back(PI(w, i));
    }

    if(p.size() == 0) return false;

    // get all pairs of edges (x, y) as the target of subsetsum
    // produces equations of type x+y=a+b+...
    vector<PI> ppi;
    vector<PI> ppw;
    for(int i = 0; i < p.size(); i++)
    {
	for(int j = i + 1; j < p.size(); j++)
	{
	    int w = p[i].first + p[j].first;
	    ppw.push_back(PI(w, ppw.size()));
	    ppi.push_back(PI(p[i].second, p[j].second));
	}
    }

    subsetsum sss1(p, ppw);
    sss1.solve();

    nested_graph nt(gr);
    for(int i = 0; i < sss1.eqns.size(); i++)
    {
	equation &eqn = sss1.eqns[i];
	assert(eqn.s.size() == 1);
	if(eqn.t.size() == 1) continue;

	int si = eqn.s[0];
	eqn.w = ppw[si].first;
	eqn.s.clear();
	eqn.s.push_back(ppi[si].first);
	eqn.s.push_back(ppi[si].second);

	set<int> ss(eqn.t.begin(), eqn.t.end());
	if(ss.find(eqn.s[0]) != ss.end()) continue;
	if(ss.find(eqn.s[1]) != ss.end()) continue;

	if(verify_equation_mergable(eqn) == false) continue;
	if(verify_equation_nontrivial(eqn) == false) continue;

	if(fast_mode == true)
	{
	    for(int j = 0; j < eqn.s.size(); j++)
	    {
		for(int k = 0; k < eqn.t.size(); k++)
		{
		    bool b = check_adjacent_mergable(eqn.s[j], eqn.t[k], nt);
		    if(b == false) eqn.d++;
		    else eqn.a++;
		}
	    }
	    resolve_vertex_with_equation(eqn);
	    eqns.push_back(eqn);
	    if(eqn.f == 3) return true;
	}
	else
	{
	    catfish sc;
	    save(sc);
	    resolve_equation(eqn);
	    load(sc);
	    eqns.push_back(eqn);
	    if(eqn.f == 2 && eqn.d == 0) return true;
	}
    }

    return false;
}

bool catfish::resolve_vertex_with_equations(vector<equation> &eqns)
{
    for(int i = 0; i < eqns.size(); i++)
    {
	equation &eqn = eqns[i];
	resolve_vertex_with_equation(eqn);
	if(eqn.f == 3) return true;
    }
    return false;
}

bool catfish::resolve_vertex_with_equation(equation &eqn)
{
    resolve_vertex_with_equation1(eqn);
    resolve_vertex_with_equation2(eqn);
    if(eqn.f == 3) return true;

    vector<int> t = eqn.s;
    eqn.s = eqn.t;
    eqn.t = t;
    resolve_vertex_with_equation1(eqn);
    resolve_vertex_with_equation2(eqn);
    if(eqn.f == 3) return true;
    return false;
}

bool catfish::resolve_vertex_with_equation1(equation &eqn)
{
    if(eqn.s.size() == 0) return false;	
	
    SE sex;
    edge_descriptor ex = i2e[eqn.s[0]];
    sex.insert(ex);
    int xt = ex->target();
    for(int i = 1; i < eqn.s.size(); i++)
    {
	edge_descriptor ex = i2e[eqn.s[i]];
	sex.insert(ex);
	if(ex->target() != xt) return false;
    }

    for(int i = 0; i < eqn.t.size(); i++)
    {
	edge_descriptor ey = i2e[eqn.t[i]];
	int ys = ey->source();
	if(gr.check_path(xt, ys) == false) return false;
    }

    assert(gr.in_degree(xt) >= 1);
    if(gr.in_degree(xt) != 1 + eqn.s.size()) return false;

    edge_descriptor ex2 = null_edge;
    edge_iterator it1, it2;
    for(tie(it1, it2) = gr.in_edges(xt); it1 != it2; it1++)
    {
	edge_descriptor e = (*it1);
	if(sex.find(e) != sex.end()) continue;
	ex2 = e;
    }
    assert(ex2 != null_edge);
    double w2 = gr.get_edge_weight(ex2);

    SE sey;
    for(int i = 0; i < eqn.t.size(); i++)
    {
	edge_descriptor ey = i2e[eqn.t[i]];
	sey.insert(ey);
	int ys = ey->source();
	SE se;
	gr.bfs_reverse(ys, se);
	sey.insert(se.begin(), se.end());
    }

    VE ve;
    double w1 = 0;
    for(tie(it1, it2) = gr.out_edges(xt); it1 != it2; it1++)
    {
	edge_descriptor e = (*it1);
	if(sey.find(e) != sey.end()) continue;
	ve.push_back(e);
	w1 += gr.get_edge_weight(e);
    }

    if(ve.size() <= 0) return false;
    if(ve.size() >= 1 && w1 > w2 + SMIN) return false;

    int ei = e2i[ex2];
    for(int i = 0; i < ve.size(); i++)
    {
	edge_descriptor e = ve[i];
	double ww = gr.get_edge_weight(e);
	int ee = split_edge(ei, ww);

	printf("resolve vertex with equation: (%d, %d)\n", ei, e2i[e]);

	merge_adjacent_equal_edges(ei, e2i[e]);
	ei = ee;
    }

    eqn.f = 3;
    return true;
}

bool catfish::resolve_vertex_with_equation2(equation &eqn)
{
    if(eqn.t.size() == 0) return false;	

    SE sey;
    edge_descriptor ey = i2e[eqn.t[0]];
    sey.insert(ey);
    int ys = ey->source();
    for(int i = 1; i < eqn.t.size(); i++)
    {
	edge_descriptor ey = i2e[eqn.t[i]];
	sey.insert(ey);
	if(ey->source() != ys) return false;
    }

    for(int i = 0; i < eqn.s.size(); i++)
    {
	edge_descriptor ex = i2e[eqn.s[i]];
	int xt = ex->target();
	if(gr.check_path(xt, ys) == false) return false;
    }

    assert(gr.out_degree(ys) >= 1);
    if(gr.out_degree(ys) != 1 + eqn.t.size()) return false;

    edge_iterator it1, it2;
    edge_descriptor ey2 = null_edge;
    for(tie(it1, it2) = gr.out_edges(ys); it1 != it2; it1++)
    {
	edge_descriptor e = (*it1);
	if(sey.find(e) != sey.end()) continue;
	ey2 = e;
    }
    assert(ey2 != null_edge);
    double w2 = gr.get_edge_weight(ey2);

    SE sex;
    for(int i = 0; i < eqn.s.size(); i++)
    {
	edge_descriptor ex = i2e[eqn.s[i]];
	sex.insert(ex);
	int xt = ex->target();
	SE se;	
	gr.bfs(xt, se);
	sex.insert(se.begin(), se.end());
    }
    VE ve;
    double w1 = 0;
    for(tie(it1, it2) = gr.in_edges(ys); it1 != it2; it1++)
    {
	edge_descriptor e = (*it1);
	if(sex.find(e) != sex.end()) continue;
	ve.push_back(e);
	w1 += gr.get_edge_weight(e);
    }

    if(ve.size() <= 0) return false;
    if(ve.size() >= 1 && w1 > w2 + SMIN) return false;

    int ei = e2i[ey2];
    for(int i = 0; i < ve.size(); i++)
    {
	edge_descriptor e = ve[i];
	double ww = gr.get_edge_weight(e);
	int ee = split_edge(ei, ww);

	printf("resolve vertex with equation: (%d, %d)\n", ei, e2i[e]);

	merge_adjacent_equal_edges(ei, e2i[e]);
	ei = ee;
    }

    eqn.f = 3;
    return true;
}

int catfish::resolve_equation(equation &eqn)
{
    vector<int> s = eqn.s;
    vector<int> t = eqn.t;
    eqn.a = eqn.d = 0;
    eqn.f = resolve_equation(s, t, eqn.a, eqn.d);
    return 0;
}

int catfish::resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md)
{
    if(s.size() == 0 && t.size() == 0) return 2;

    assert(s.size() >= 1);
    assert(t.size() >= 1);

    for(int i = 0; i < s.size(); i++)
    {
	for(int j = 0; j < t.size(); j++)
	{
	    int x = s[i];
	    int y = t[j];
	    assert(i2e[x] != null_edge);
	    assert(i2e[y] != null_edge);

	    vector<PI> p;
	    if(check_adjacent_mergable(x, y, p) == false) continue;

	    build_adjacent_edges(p);

	    double wx = gr.get_edge_weight(i2e[x]);
	    double wy = gr.get_edge_weight(i2e[y]);
	    double ww = (wx <= wy) ? wx : wy;

	    vector<int> v;
	    v.push_back(x);
	    v.push_back(y);

	    vector<int> vv;
	    split_merge_path(v, ww, vv);

	    ma++;

	    //// printf("merge (adjacent) edge pair (%d, %d)\n", x, y);

	    if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
	    if(i2e[vv[1]] == null_edge) assert(vv[1] == y);

	    if(vv[0] == x) s.erase(s.begin() + i);
	    else s[i] = vv[0];

	    if(vv[1] == y) t.erase(t.begin() + j);
	    else t[j] = vv[1];

	    int f = resolve_equation(s, t, ma, md);
	    if(f == 2) return 2;
	    else return 1;
	}
    }

    set<int> ss;
    for(int i = 0; i < s.size(); i++)
    {
	assert(ss.find(s[i]) == ss.end());
	ss.insert(s[i]);
    }
    for(int i = 0; i < t.size(); i++)
    {
	assert(ss.find(t[i]) == ss.end());
	ss.insert(t[i]);
    }

    for(int i = 0; i < s.size(); i++)
    {
	for(int j = 0; j < t.size(); j++)
	{
	    int x = s[i];
	    int y = t[j];
	    assert(i2e[x] != null_edge);
	    assert(i2e[y] != null_edge);

	    double wx = gr.get_edge_weight(i2e[x]);
	    double wy = gr.get_edge_weight(i2e[y]);
	    double ww = (wx <= wy) ? wx : wy;

	    VE p;
	    int l = check_distant_mergable(x, y, ww, p);
	    if(l < 0) continue;

	    assert(l >= 2);
	    assert(p.size() >= 2);
	    assert(i2e[x] == p[0]);
	    assert(i2e[y] == p[p.size() - 1]);

	    bool c = false;
	    for(int k = 1; k < p.size() - 1; k++)
	    {
		int e = e2i[p[k]];
		if(ss.find(e) != ss.end()) c = true;
		if(c == true) break;
	    }

	    if(c == true) continue;

	    vector<int> vv;
	    split_merge_path(p, ww, vv);

	    //// printf("connect (distant) edge pair (%d, %d)\n", x, y);

	    md++;

	    int n = vv.size() - 1;

	    if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
	    if(i2e[vv[n]] == null_edge) assert(vv[n] == y);

	    if(vv[0] == x) s.erase(s.begin() + i);
	    else s[i] = vv[0];

	    if(vv[n] == y) t.erase(t.begin() + j);
	    else t[j] = vv[n];

	    int f = resolve_equation(s, t, ma, md);
	    if(f == 2) return 2;
	    else return 1;
	}
    }

    return 0;
}

bool catfish::check_adjacent_mergable(int ex, int ey, nested_graph &nt)
{
    assert(i2e[ex] != null_edge);
    assert(i2e[ey] != null_edge);

    int xs = i2e[ex]->source();
    int xt = i2e[ex]->target();
    int ys = i2e[ey]->source();
    int yt = i2e[ey]->target();

    if(xt == ys) return true;
    if(yt == xs) return true;

    bool b = false;
    vector<PI> xp, yp;
    if(gr.check_path(i2e[ex], i2e[ey])) b = nt.link(xs, xt, ys, yt, xp, yp);
    else if(gr.check_path(i2e[ey], i2e[ex])) b = nt.link(ys, yt, xs, xt, yp, xp);
    else return false;

    return b;
}

bool catfish::check_adjacent_mergable(int ex, int ey, vector<PI> &p)
{
    assert(i2e[ex] != null_edge);
    assert(i2e[ey] != null_edge);

    int xs = i2e[ex]->source();
    int xt = i2e[ex]->target();
    int ys = i2e[ey]->source();
    int yt = i2e[ey]->target();

    if(xt == ys) return true;
    if(yt == xs) return true;

    vector<PI> xp, yp;
    bool b = false;

    nested_graph nt(gr);

    if(gr.check_path(i2e[ex], i2e[ey])) b = nt.link(xs, xt, ys, yt, xp, yp);
    else if(gr.check_path(i2e[ey], i2e[ex])) b = nt.link(ys, yt, xs, xt, yp, xp);
    else return false;
	
    if(b == false) return false;

    p = xp;
    p.insert(p.end(), yp.begin(), yp.end());

    return true;
}

int catfish::check_distant_mergable(int x, int y, double w)
{
    VE p;
    return check_distant_mergable(x, y, w, p);
}

int catfish::check_distant_mergable(int x, int y, double w, VE &p)
{
    p.clear();

    assert(i2e[x] != null_edge);
    assert(i2e[y] != null_edge);

    edge_descriptor xx = i2e[x];
    edge_descriptor yy = i2e[y];

    int xs = (xx)->source();
    int xt = (xx)->target();
    int ys = (yy)->source();
    int yt = (yy)->target();

    if(gr.check_path(yt, xs) == true)
    {
	int r = check_distant_mergable(y, x, w, p);
	reverse(p);
	return r;
    }

    if(gr.check_path(xt, ys) == false) return -1;

    int l = gr.compute_shortest_path_w(xt, ys, w, p);
    if(l < 0) return -1;

    p.insert(p.begin(), xx);
    p.insert(p.end(), yy);

    assert(p.size() == l + 2);
    return l + 2;
}

int catfish::build_adjacent_edges(const vector<PI> &p)
{
    for(int i = 0; i < p.size(); i++)
    {
	int x = p[i].first;
	int y = p[i].second;
	if(y == -1)
	{
	    int l = gr.compute_in_partner(x);
	    int r = gr.compute_out_partner(x);
	    gr.exchange(l, x, r);
	}
	else
	{
	    gr.rotate(x, y);
	}
    }
    return 0;
}

bool catfish::decompose_trivial_vertices()
{
    bool flag = false;
    for(int i = 1; i < gr.num_vertices() - 1; i++)
    {
	if(gr.degree(i) == 0) continue;
	if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;

	printf("decompose trivial vertex %d\n", i);

	equation eqn(0);
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
	{
	    int e = e2i[*it1];
	    eqn.s.push_back(e);
	}
	for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
	{
	    int e = e2i[*it1];
	    eqn.t.push_back(e);
	}

	printf(" ");
	// eqn.print(i);
	print_equation(eqn, i);

	resolve_equation(eqn);

	flag = true;
    }
    return flag;
}

int catfish::greedy_decompose()
{
    while(true)
    {
	VE v;
	vector<int> vv;
	double w = gr.compute_maximum_path_w(v);
	if(w <= 0.0) break;
	int e = split_merge_path(v, w, vv);
	collect_path(e);
    }
    return 0;
}

int catfish::collect_existing_st_paths()
{
    for(int i = 0; i < i2e.size(); i++)
    {
	if(i2e[i] == null_edge) continue;
	if(i2e[i]->source() != 0) continue;
	if(i2e[i]->target() != gr.num_vertices() - 1) continue;
	collect_path(i);
    }
    return 0;
}

int catfish::collect_path(int e)
{
    assert(i2e[e] != null_edge);
    assert(i2e[e]->source() == 0);
    assert(i2e[e]->target() == gr.num_vertices() - 1);

    assert(mev.find(i2e[e]) != mev.end());
    vector<int> v = mev[i2e[e]];
    sort(v.begin(), v.end());

    assert(v[0] == 0);
    assert(v[v.size() - 1] < gr.num_vertices() - 1);
    v.push_back(gr.num_vertices() - 1);

    path p;
    p.abd = gr.get_edge_weight(i2e[e]);
    p.v = v;
    paths.push_back(p);

    gr.remove_edge(i2e[e]);
    e2i.erase(i2e[e]);
    i2e[e] = null_edge;

    return 0;
}

int catfish::remove_empty_edges()
{
    for(int i = 0; i < i2e.size(); i++)
    {
	if(i2e[i] == null_edge) continue;
	double w = gr.get_edge_weight(i2e[i]);
	if(w >= 1) continue;
	assert(w <= 0);
	e2i.erase(i2e[i]);
	gr.remove_edge(i2e[i]);
	i2e[i] = null_edge;
    }
    return 0;
}

int catfish::print()
{
    int n = 0;
    for(int i = 0; i < gr.num_vertices(); i++) 
    {
	if(gr.degree(i) >= 1) n++;
    }

    int p1 = gr.compute_num_paths();
    int p2 = gr.compute_decomp_paths();
    printf("statistics: %lu edges, %d vertices, total %d paths, %d required\n", gr.num_edges(), n, p1, p2);
    printf("finish round %d\n\n", round);

    if(output_tex_files == true)
    {
	draw_splice_graph(name + "." + tostring(round) + ".tex");
	nested_graph nt(gr);
	nt.draw(name + "." + tostring(round) + ".nt.tex");
    }

    round++;

    return 0;
}

int catfish::draw_splice_graph(const string &file) 
{
    MIS mis;
    char buf[10240];
    for(int i = 0; i < gr.num_vertices(); i++)
    {
	double w = gr.get_vertex_weight(i);
	sprintf(buf, "%d:%.0lf", i, w);
	mis.insert(PIS(i, buf));
    }

    MES mes;
    for(int i = 0; i < i2e.size(); i++)
    {
	if(i2e[i] == null_edge) continue;
	double w = gr.get_edge_weight(i2e[i]);
	sprintf(buf, "%d:%.0lf", i, w);
	//sprintf(buf, "%d", i);
	mes.insert(PES(i2e[i], buf));
    }
    gr.draw(file, mis, mes, 4.5);
    return 0;
}

int catfish::write(const string &file)
{
    ofstream fout(file.c_str());
    write(fout);
    fout.close();
    return 0;
}

int catfish::write(ofstream &fout)
{
    for(int i = 0; i < paths.size(); i++)
    {
	//fout <<"path " << i + 1 << ", weight = " << (int)(paths[i].abd) << ", vertices = ";
	fout << (int)(paths[i].abd) << " ";
	for(int k = 0; k < paths[i].v.size(); k++)
	{
	    fout << paths[i].v[k] << " ";
	}
	fout << endl;
    }
    return 0;
}

void catfish::gurobi_setup_base_model()
{
    gurobi_get_edge_weights(); // prepare i2w
    int ne = i2w.size();
    gurobi_vars = std::vector<GRBVar*>(ne);
    for(int i = 0; i < ne; ++i)
    {
	edge_descriptor ei = i2e[i];
	if(ei == null_edge) continue;
	
	// lb(0.0), ub, obj(0.0), type(continuous), names(default), count
	gurobi_vars[i] = gurobi_model->addVars(NULL, i2w.data(), NULL,
					       NULL, NULL, ne);

	// LP constraint 1: flow saturate edge i
	gurobi_vars[i][i].set(GRB_DoubleAttr_LB, i2w[i]);
    }

    // LP constraint 2: symmetry
    for(int i = 1; i < ne; ++i)
    {
	if(i2e[i] == null_edge) continue;
	for(int j = 0; j < i; ++j)
	{
	    if(i2e[j] == null_edge) continue;
	    gurobi_model->addConstr(gurobi_vars[j][i],
				    GRB_EQUAL,
				    gurobi_vars[i][j]);
	}
    }
    
    // LP constraint 3:
    // flow conservation constraints, exclude s(0) and t(num_vertices-1)
    edge_iterator it, itEd;
    for(int i = 1; i < gr.num_vertices() - 1; ++i)
    {
	std::vector<GRBLinExpr> in(ne);
	for(tie(it, itEd) = gr.in_edges(i); it != itEd; ++it)
	{
	    int idx = e2i[*it];
	    for(int j = 0; j < ne; ++j)
	    {
		if(i2e[j] != null_edge)
		{
		    in[j] += gurobi_vars[j][idx];
		}
	    }
	}
	std::vector<GRBLinExpr> out(ne);
	for(tie(it, itEd) = gr.out_edges(i); it != itEd; ++it)
	{
	    int idx = e2i[*it];
	    for(int j = 0; j < ne; ++j)
	    {
		if(i2e[j] != null_edge)
		{
		    out[j] += gurobi_vars[j][idx];
		}
	    }
	}

	for(int j = 0; j < ne; ++j)
	{
	    gurobi_model->addConstr(in[j], GRB_EQUAL, out[j]);
	}
    }

    // LP constraint 4: the i-th flow must all go through edge i
    std::vector<GRBLinExpr> source(ne);
    for(tie(it, itEd) = gr.out_edges(0); it != itEd; ++it)
    {
	int idx = e2i[*it];
	for(int i = 0; i < ne; ++i)
	{
	    if(i2e[i] != null_edge)
	    {
		source[i] += gurobi_vars[i][idx];
	    }
	}
    }
    for(int i = 0; i < ne; ++i)
    {
	if(i2e[i] != null_edge)
	{
	    gurobi_model->addConstr(source[i], GRB_EQUAL, i2w[i]);
	}
    }
}

bool catfish::gurobi_add_equation(equation& eqn)
{
    std::vector<GRBConstr> new_constraints;
    int ne = i2w.size();
    
    for(int i = 0; i < ne; ++i)
    {
	if(i2e[i] != null_edge)
	{
	    GRBLinExpr lhs;
	    for(int x : eqn.s)
	    {
		lhs += gurobi_vars[x][i];
	    }

	    // LP constraint 1:
	    // combined flow of the equation does not exceed the weight of each edge
	    new_constraints.push_back(gurobi_model->addConstr(lhs,
							      GRB_LESS_EQUAL,
							      i2w[i]));

	    GRBLinExpr rhs;
	    for(int x : eqn.t)
	    {
		rhs += gurobi_vars[x][i];
	    }

	    // LP constraint 2:
	    // for each edge, combined flow on lhs equals that on rhs
	    new_constraints.push_back(gurobi_model->addConstr(lhs,
							      GRB_EQUAL,
							      rhs));
	}
    }

    if(gurobi_is_feasible())
    {
	return true;
    }
    else
    {
	for(int i = 0; i < new_constraints.size(); ++i)
	{
	    gurobi_model->remove(new_constraints[i]);
	}
	// gurobi_model->update();
	return false;
    }
}

bool catfish::gurobi_is_feasible()
{
    gurobi_model->optimize();
    int result = gurobi_model->get(GRB_IntAttr_Status);

    if(result == GRB_INF_OR_UNBD)
    {
	gurobi_model->set(GRB_IntParam_Presolve, 0);
	gurobi_model->optimize();
	result = gurobi_model->get(GRB_IntAttr_Status);
    }

    if(result == GRB_OPTIMAL)
    {
	return true;
    }
    else if(result == GRB_INFEASIBLE)
    {
	return false;
    }
    else
    {
	std::cerr << "Gurobi stopped with status: " << result << std::endl;
	return false;
    }
}

void catfish::gurobi_get_edge_weights()
{
    int ne = i2e.size();
    i2w = std::vector<double>(ne);
    for(int i = 0; i < ne; ++i)
    {
	i2w[i] = (i2e[i] == null_edge) ? 0.0 : gr.get_edge_weight(i2e[i]);
    }
}

bool catfish::decompose_with_equations_induced_merge()
{
    //eqns filled by decompose_with_equations(X)
    if(eqns.empty()) return false;

    if(!gurobi_base_model_built)
    {
	try
	{
	    gurobi_model = new GRBModel(*gurobi_env);
	    gurobi_model->set(GRB_IntParam_OutputFlag, 0); // no logging output
	    gurobi_model->set(GRB_IntParam_Threads, GRB_THREAD_LIMIT); // limit number of threads gurobi can use
	    gurobi_setup_base_model();
	}
	catch(GRBException& e)
	{
	    std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	    std::cerr << e.getMessage() << std::endl;
	    throw;
	}
	gurobi_base_model_built = true;
    }

    printf("Check feasibility of equations\n");
    std::vector<equation> feasible_eqns;
    try
    {
	for(int i = 0; i < eqns.size(); ++i)
	{
	    bool feasible = gurobi_add_equation(eqns[i]);

	    printf("%10s ", feasible ? "feasible" : "infeasible");
	    // eqns[i].print(i);
	    print_equation(eqns[i], i);
	    if(has_truth && equation_comply_with_truth(eqns[i])) printf("  OK with ground truth\n");
		
	    if(feasible)
	    {
		// print solution
		/*
		int ne = gurobi_vars.size();
		for(int x : eqns[i].s)
		{
		    edge_descriptor ex = i2e[x];
		    printf("    flow through edge %d (%d, %d):\n",
			   x, ex->source(), ex->target());
		    for(int i = 0; i < ne; ++i)
		    {
			if(i2e[i] == null_edge) continue;
		        double sol = gurobi_vars[x][i].get(GRB_DoubleAttr_X);
			if(sol > 0)
			{
			    edge_descriptor ei = i2e[i];
			    printf("        %d (%d, %d) [%.1f]\n",
				   i, ei->source(), ei->target(), sol);
			}
		    }
		}
		for(int x : eqns[i].t)
		{
		    edge_descriptor ex = i2e[x];
		    printf("    flow through edge %d (%d, %d):\n",
			   x, ex->source(), ex->target());
		    for(int i = 0; i < ne; ++i)
		    {
			if(i2e[i] == null_edge) continue;
		        double sol = gurobi_vars[x][i].get(GRB_DoubleAttr_X);
			if(sol > 0)
			{
			    edge_descriptor ei = i2e[i];
			    printf("        %d (%d, %d) [%.1f]\n",
				   i, ei->source(), ei->target(), sol);
			}
		    }
		}
		*/

		feasible_eqns.push_back(std::move(eqns[i]));
	    }
	}

	eqns.clear();
	
	if(feasible_eqns.empty()) return false;
	else
	{
	    // if the last tested equation makes the model infeasible,
	    // re-run the model to get a feasible solution to be used
	    // for finding mergeable edges
	    int result = gurobi_model->get(GRB_IntAttr_Status);
	    if(result != GRB_OPTIMAL)
	    {
		gurobi_is_feasible();
	    }
	    return gurobi_find_mergeable_edges();
	}
    }
    catch(GRBException& e)
    {
	std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	std::cerr << e.getMessage() << std::endl;
	throw;
    }
}

void catfish::gurobi_clear_model()
{
    if(gurobi_base_model_built)
    {
	for(int i = 0; i < gurobi_vars.size(); ++i)
	{
	    delete[] gurobi_vars[i];
	}
	delete gurobi_model;
	gurobi_base_model_built = false;
    }
}

bool catfish::gurobi_find_mergeable_edges()
{
    if(output_tex_files) gurobi_print_current_solution();
    printf("Use feasible model to find mergeable edges\n");
    bool any = false;
    int ne = gurobi_vars.size();
    
    // save feasible solution
    std::vector<std::vector<double> > sol(ne, std::vector<double>(ne));
    for(int i = 0; i < ne; ++i)
    {
	if(i2e[i] == null_edge) continue;
	for(int j = 0; j < ne; ++j)
	{
	    if(i2e[j] == null_edge) continue;
	    sol[i][j] = gurobi_vars[i][j].get(GRB_DoubleAttr_X);
	}
    }

    // mark if an edge has been involved in a merge
    std::vector<bool> merged(ne, false);

    // for each flow that saturates one edge e, check if there are adjacent edges
    // that must have flow w(e), if found, they can be merged
    for(int i = 0; i < ne; ++i)
    {
	if(merged[i]) continue;
	edge_descriptor e = i2e[i];
	if(e == null_edge) continue;
	
	int se = e->source();
	int te = e->target();
	edge_iterator it, itEd;
	for(tie(it, itEd) = gr.in_edges(se); it != itEd; ++it)
	{
	    int x = e2i[*it];
	    if(x >= ne || merged[x]) continue;
	    if(fabs(sol[i][x] - i2w[i]) < gurobi_slack)
	    {
		GRBConstr tmp = gurobi_model->addConstr(gurobi_vars[i][x],
							GRB_LESS_EQUAL,
							i2w[i] - gurobi_slack);
		bool still_feasible = gurobi_is_feasible();	
		gurobi_model->remove(tmp);

		if(!still_feasible)
		{
		    // all flow through i must go through x
		    printf("  All flow through %2d (%2d, %2d) must use %2d (%2d, %2d), merged\n",
			   i, se, te, x, (*it)->source(), se);
		    print_edge_with_backtrack(i);
		    print_edge_with_backtrack(x);

		    merged[i] = true;
		    merged[x] = true;
		    merge_adjacent_edges(x, i);
		    any = true;
		    break;
		}
	    }
	}

	if(merged[i]) continue;
	for(tie(it, itEd) = gr.out_edges(te); it != itEd; ++it)
	{
	    int x = e2i[*it];
	    if(x >= ne || merged[x]) continue;
	    if(fabs(sol[i][x] - i2w[i]) < gurobi_slack)
	    {
		GRBConstr tmp = gurobi_model->addConstr(gurobi_vars[i][x],
							GRB_LESS_EQUAL,
							i2w[i] - gurobi_slack);
		bool still_feasible = gurobi_is_feasible();
		gurobi_model->remove(tmp);

		if(!still_feasible)
		{
		    // all flow through i must go through x
		    printf("  All flow through %2d (%2d, %2d) must use %2d (%2d, %2d), merged\n",
			   i, se, te, x, te, (*it)->target());
		    print_edge_with_backtrack(i);
		    print_edge_with_backtrack(x);

		    merged[i] = true;
		    merged[x] = true;
		    merge_adjacent_edges(i, x);
		    any = true;
		    break;
		}
	    }
	}
    }

    // for each flow that saturates one edge e, find if any edge must use w(e),
    // if found, check if they can be made adjacent by reverse operations
    // only do at most one such merge per round, so the nested graph does not
    // change, build it here and use it to test if edges can be made adjacent.
    nested_graph nt(gr);
    for(int i = 0; i < ne; ++i)
    {
	if(merged[i]) continue;
	edge_descriptor e = i2e[i];
	if(e == null_edge) continue;
	
	int se = e->source();
	int te = e->target();
	for(int j = 0; j < ne; ++j)
	{
	    if(j == i) continue;
	    if(fabs(sol[i][j] - i2w[i]) < gurobi_slack && !merged[j])
	    {
		GRBConstr tmp = gurobi_model->addConstr(gurobi_vars[i][j],
							GRB_LESS_EQUAL,
							i2w[i] - gurobi_slack);
		bool still_feasible = gurobi_is_feasible();
		gurobi_model->remove(tmp);

		if(!still_feasible)
		{   
		    edge_descriptor ej = i2e[j];
		    // all flow through i must go through x
		    printf("  All flow through %2d (%2d, %2d) must use %2d (%2d, %2d), ",
			   i, se, te, j, ej->source(), ej->target());

		    bool result = check_adjacent_mergable(i, j, nt);
		    if(result)
		    {
			printf("can be made adjacent, merged\n");
		    }
		    else
		    {
			printf("cannot be made adjacent\n");
		    }
		    print_edge_with_backtrack(i);
		    print_edge_with_backtrack(j);
		    if(has_truth && equation_comply_with_truth(
			   equation(std::vector<int>{i}, std::vector<int>{j})))
		    {
			printf("    OK with ground truth\n");
		    }

		    if(result)
		    {
			// copied from catfish::resolve_equation
		        vector<PI> p;
			check_adjacent_mergable(i, j, p);
			build_adjacent_edges(p);
			vector<int> v;
			v.push_back(i);
			v.push_back(j);
			vector<int> vv;
			double ww = i2w[i];
			split_merge_path(v, ww, vv);

		        return true;
		    }
		}

	    }
	}
    }
    
    
    return any;
}


void catfish::gurobi_print_current_solution()
{
    printf("Current LP solution:\n");
    for(int i = 0; i < i2e.size(); ++i)
    {
	edge_descriptor ei = i2e[i];
	if(ei == null_edge) continue;
	printf("  Flow through edge %d (%d, %d):\n",
	       i, ei->source(), ei->target());

	for(int j = 0; j < i2e.size(); ++j)
	{
	    if(i2e[j] == null_edge) continue;
	    double sol = gurobi_vars[i][j].get(GRB_DoubleAttr_X);
	    if(fabs(sol) >= SMIN)
	    {
		edge_descriptor ej = i2e[j];
		printf("    %2d (%2d, %2d) [%6.1f/%6.1f]\n",
		       j, ej->source(), ej->target(), sol,
		       gr.get_edge_weight(ej));
	    }
	}
    }
}

void catfish::print_equation(const equation& eqn, int idx)
{
    eqn.print(idx);

    for(int x : eqn.s)
    {
	print_edge_with_backtrack(x);
    }
    for(int x : eqn.t)
    {
	print_edge_with_backtrack(x);
    }
}

void catfish::build_truth(const std::vector<std::pair<double, std::vector<int>>>& truth_paths)
{
    truth.reserve(truth_paths.size());
    for(const auto& p : truth_paths)
    {
	truth.emplace_back(p.first);
	int a = p.second[0];
	assert(a == 0);
	for(int i = 1; i < p.second.size(); ++i)
	{
	    int b = p.second[i];
	    PEB result = gr.edge(a, b);
	    assert(result.second);
	    truth.back().add_edge(result.first);
	    a = b;
	}
	assert(a == gr.num_vertices() - 1);
    }
}

bool catfish::equation_comply_with_truth(const equation& eqn)
{   
    // one entry for each edge in the lhs of the equation
    // record ground truth path ids that uses this edge
    // but not any other in the lhs
    std::vector<std::vector<int>> lhs_independent_paths(eqn.s.size());
    for(int i = 0; i < truth.size(); ++i)
    {
	int x = -1;
	for(int j = 0; j < eqn.s.size(); ++j)
	{
	    if(truth[i].contains_subseq(i2backtrack[eqn.s[j]]))
	    {
		if(x >= 0)
		{
		    x = -1;
		    break; // path contains multiple edges in lhs, discard
		}
		else
		{
		    x = j;
		}
	    }
	}
	if(x >= 0)
	{
	    lhs_independent_paths[x].push_back(i);
	}
    }

    std::map<int, std::pair<int, int>> p2e;
    // find true paths that go through one edge in lhs and one edge in rhs
    for(int i = 0; i < lhs_independent_paths.size(); ++i)
    {
	for(int j : lhs_independent_paths[i])
	{
	    int x = -1;
	    for(int k = 0; k < eqn.t.size(); ++k)
	    {
		if(truth[j].contains_subseq(i2backtrack[eqn.t[k]]))
		{
		    if(x >= 0)
		    {
			x = -1;
			break; // path contains multiple edges in rhs
		    }
		    else
		    {
			x = k;
		    }
		}
	    }
	    if(x >= 0) // path j contains a single edge i in lhs and a single edge x in rhs
	    {
		auto result = p2e.emplace(j, std::pair<int, int>{i, x});
		assert(result.second); // insertion should always be successful
	    }
	}
    }

    int eqn_w_l = 0;
    for(int x : eqn.s)
    {
	eqn_w_l += (int)gr.get_edge_weight(i2e[x]);
    }
    int eqn_w_r = 0;
    for(int x : eqn.t)
    {
	eqn_w_r += (int)gr.get_edge_weight(i2e[x]);
    }
    // since we may pass in two edges to be merged as an equation, they
    // may not have the same weight
    int eqn_w = eqn_w_l < eqn_w_r ? eqn_w_l : eqn_w_r; 
    
    int max_flow = get_max_flow(eqn.s, eqn.t, p2e);
    if(eqn_w != max_flow)
    {
	printf("  do not comply with ground truth, truth supports %d, equation weight %d\n",
	       max_flow, eqn_w);
	return false;
    }

    return true;
}

void catfish::print_edge_with_backtrack(int i)
{
    edge_descriptor e = i2e[i];
    if(e == null_edge) printf("    edge %2d (changed):", i);
    else printf("    edge %2d (%2d, %2d) [%.1f]:", i, e->source(), e->target(), gr.get_edge_weight(e));
    for(edge_descriptor x : i2backtrack[i])
    {
	auto search = original_e2st.find(x);
	assert(search != original_e2st.end());
	printf(" (%2d, %2d)", search->second.first, search->second.second);
    }
    printf("\n");
}


//adjusted from @sehe's answer here https://stackoverflow.com/questions/77601238/using-boostpush-relabel-max-flow
int catfish::get_max_flow(const std::vector<int>& lhs,
			  const std::vector<int>& rhs,
			  const std::map<int, std::pair<int, int>>& p2e)
{
    using namespace boost;
    typedef adjacency_list_traits<vecS, vecS, directedS> boost_traits;
    typedef boost_traits::vertex_descriptor boost_v;
    typedef boost_traits::edge_descriptor boost_e;
    struct boost_v_props {};
    struct boost_e_props {
	int edge_capacity = 0;
	boost_e reverse_edge  = {};

	boost_e_props(long c = 0, boost_e rev_e = {}) :edge_capacity(c), reverse_edge(rev_e){};
    };
    typedef adjacency_list<vecS, vecS, directedS, boost_v_props, boost_e_props> boost_graph;

    int n = lhs.size() + rhs.size() + 2;
    boost_graph g(n);

    auto ecap = get(&boost_e_props::edge_capacity, g);
    auto reve = get(&boost_e_props::reverse_edge, g);

    auto insert = [&](boost_v s, boost_v t, boost_e_props p) {
        auto fwd  = add_edge(s, t, p, g).first;
        auto bck  = add_edge(t, s, {}, g).first;
        reve[fwd] = bck;
        reve[bck] = fwd;
    };

    boost_v source = n - 2, sink = n - 1;
    
    for(int i = 0; i < lhs.size(); ++i)
    {
        insert(source, i, {(int)gr.get_edge_weight(i2e[lhs[i]])});
    }
    for(int i = 0; i < rhs.size(); ++i)
    {
	insert(i + lhs.size(), sink, {(int)gr.get_edge_weight(i2e[rhs[i]])});
    }

    for(auto it = p2e.begin(); it != p2e.end(); ++it)
    {
	insert(it->second.first, it->second.second + lhs.size(),
	       {(int)truth[it->first].weight});
    }

    std::map<boost_e, int> residuals;
    return push_relabel_max_flow(g, source, sink,
				 capacity_map(ecap)
				 .reverse_edge_map(reve)
				 .residual_capacity_map(make_assoc_property_map(residuals)));
}

bool catfish::decompose_adjacent_equal_edges()
{
    bool any = false;
    for(int i = 1; i < gr.num_vertices() - 1; ++i)
    {
	if(gr.degree(i) == 0) continue;
	//otherwise should have been handled by decompose_trivial_vertices
	assert(gr.in_degree(i) > 1);
	assert(gr.out_degree(i) > 1);

	std::unordered_map<int, std::vector<int>> w2in;
	edge_iterator it, itEd;
	for(tie(it, itEd) = gr.in_edges(i); it != itEd; ++it)
	{
	    int w = (int) gr.get_edge_weight(*it);
	    int ei = e2i[*it];
	    auto search = w2in.find(w);
	    if(search == w2in.end())
	    {
		w2in.emplace(w, std::vector<int>{ei});
	    }
	    else
	    {
		search->second.push_back(ei);
	    }
	}

	std::vector<PI> pairs;

	for(tie(it, itEd) = gr.out_edges(i); it != itEd; ++it)
	{
	    int w = (int) gr.get_edge_weight(*it);
	    auto search = w2in.find(w);
	    if(search != w2in.end())
	    {
		int eo = e2i[*it];
		int ei = search->second.back();
		if(search->second.size() == 1)
		{
		    w2in.erase(search);
		}
		else{
		    search->second.pop_back();
		}

		pairs.emplace_back(ei, eo);
	    }
	}

	for(std::pair<int, int>& p : pairs)
	{
	    edge_descriptor ei = i2e[p.first];
	    edge_descriptor eo = i2e[p.second];
	    printf("merge adjacent equal edges %d (%d, %d) and %d (%d, %d)\n",
		   p.first, ei->source(), ei->target(),
		   p.second, eo->source(), eo->target());
	    
	    print_edge_with_backtrack(p.first);
	    print_edge_with_backtrack(p.second);
	    if(has_truth && equation_comply_with_truth(
		   equation(std::vector<int>{p.first}, std::vector<int>{p.second})))
	    {
		printf("    OK with ground truth\n");
	    }

	    merge_adjacent_equal_edges(p.first, p.second);
	    any = true;
	}
	
    }
    return any;
}

void catfish::find_equations_lhs1(std::list<equation>& results)
{
    std::vector<PI> p;
    // find all equal edges
    std::unordered_map<int, std::vector<int>> w2i;
    for(int i = 0; i < i2e.size(); ++i)
    {
        edge_descriptor e = i2e[i];
        if(e == null_edge) continue;
        if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
        int w = (int)(gr.get_edge_weight(e));
        if(w <= 0) continue;
        p.push_back(PI(w, i));
        auto search = w2i.find(w);
        if(search == w2i.end())
        {
            w2i.emplace(w, std::vector<int>{i});
        }
        else
        {
            for(int x : search->second)
            {
                equation eqn(std::vector<int>{x}, std::vector<int>{i});
                if(verify_equation_mergable(eqn) == false) continue;
                results.push_back(eqn);
            }
            search->second.push_back(i);
        }
    }

    // find all A = B + C
    for(auto it1 = w2i.begin(); it1 != w2i.end(); ++it1)
    {
	int w1 = it1->first;
	for(auto it2 = it1; it2 != w2i.end(); ++it2)
	{
	    int w = it2->first + w1;
	    auto search = w2i.find(w);
	    if(search != w2i.end())
	    {
		for(int a : search->second)
		{
		    if(it1 == it2)
		    {
			const std::vector<int>& bc = it1->second;
			for(int bi = 0; bi < bc.size(); ++bi)
			{
			    int b = bc[bi];
			    for(int ci = bi + 1; ci < bc.size(); ++ ci)
			    {
				int c = bc[ci];
				equation eqn(std::vector<int>{a}, std::vector<int>{b, c});
				if(verify_equation_mergable(eqn) == false) continue;
				results .push_back(eqn);
			    }
			}
		    }
		    else
		    {
			for(int b : it1->second)
			{
			    for(int c : it2->second)
			    {
				equation eqn(std::vector<int>{a}, std::vector<int>{b, c});
				if(verify_equation_mergable(eqn) == false) continue;
				results .push_back(eqn);
			    }
			}
		    }
		}
	    }
	}
    }

    if(p.size() == 0) return;

    subsetsum sss(p, p);
    sss.solve();

    for(equation& eqn : sss.eqns)
    {
	assert(eqn.s.size() == 1);
	if(eqn.t.size() <= 2) continue;
	if(verify_equation_mergable(eqn) == false) continue;
	if(verify_equation_nontrivial(eqn) == false) continue;
	results.push_back(eqn);
    }
}

void catfish::find_equations_lhs2(std::list<equation>& results)
{
    std::vector<PI> p;
    for(int i = 0; i < i2e.size(); ++i)
    {
	edge_descriptor e = i2e[i];
	if(e == null_edge) continue;
	if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
	int w = (int)(gr.get_edge_weight(e));
	if(w <= 0) continue;
	p.push_back(PI(w, i));
    }

    if(p.size() == 0) return;

    // get all pairs of edges (x, y) as the target of subsetsum
    // produces equations of type x+y=a+b+...
    std::vector<PI> ppi;
    std::vector<PI> ppw;
    for(int i = 0; i < p.size(); i++)
    {
	for(int j = i + 1; j < p.size(); j++)
	{
	    int w = p[i].first + p[j].first;
	    ppw.push_back(PI(w, ppw.size()));
	    ppi.push_back(PI(p[i].second, p[j].second));
	}
    }
    
    subsetsum sss(p, ppw);
    sss.solve();

    for(equation& eqn : sss.eqns)
    {
	assert(eqn.s.size() == 1);
	if(eqn.t.size() == 1) continue; // covered in find_equations_lhs1

	int si = eqn.s[0];
	eqn.w = ppw[si].first;
	eqn.s.clear();
	eqn.s.push_back(ppi[si].first);
	eqn.s.push_back(ppi[si].second);

	std::set<int> ss(eqn.t.begin(), eqn.t.end());
	if(ss.find(eqn.s[0]) != ss.end()) continue;
	if(ss.find(eqn.s[1]) != ss.end()) continue;
	
	if(verify_equation_mergable(eqn) == false) continue;
	if(verify_equation_nontrivial(eqn) == false) continue;
	results.push_back(eqn);
    }
}


void catfish::gurobi_filter_equations(std::list<equation>& eqns)
{
    if(eqns.empty()) return;

    if(!gurobi_base_model_built)
    {
	try
	{
	    gurobi_model = new GRBModel(*gurobi_env);
	    gurobi_model->set(GRB_IntParam_OutputFlag, 0); // no logging output
	    gurobi_model->set(GRB_IntParam_Threads, GRB_THREAD_LIMIT); // limit number of threads gurobi can use
	    gurobi_setup_base_model();
	}
	catch(GRBException& e)
	{
	    std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	    std::cerr << e.getMessage() << std::endl;
	    throw;
	}
	gurobi_base_model_built = true;
    }

    printf("check feasibility of equations\n");
    try
    {
	int idx = 0;
	for(auto it = eqns.begin(); it != eqns.end();)
	{
	    bool feasible = gurobi_add_equation(*it);
	    printf("%10s ", feasible ? "feasible" : "infeasible");
	    //(*it).print(idx++);
	    print_equation(*it, idx++);
	    if(has_truth && equation_comply_with_truth(*it)) printf("  OK with ground truth\n");
	    
	    if(feasible)
	    {
		++it;
	    }
	    else
	    {
		it = eqns.erase(it);
	    }
	}
    }
    catch(GRBException& e)
    {
	std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	std::cerr << e.getMessage() << std::endl;
	throw;
    }
}

bool catfish::decompose_with_equations()
{
    std::list<equation> eqns;
    find_equations_lhs1(eqns); // A = B + ...

    gurobi_clear_model();

    int size = 2;
    for(; size <= 3; ++size)
    {
	printf("\nequations of size %d\n", size);
	if(filter_and_try_resolve_equations(size, eqns)) return true;
    }

    find_equations_lhs2(eqns); // A + B = C + D + ...

    while(!eqns.empty())
    {
	printf("\nequations of size %d\n", size);
	if(filter_and_try_resolve_equations(size, eqns)) return true;
	++size;
    }

    return false;
}


bool catfish::filter_and_try_resolve_equations(int size, std::list<equation>& eqns)
{
    // std::vector<equation> feasible_eqns;
    std::list<equation> cur_eqns;
    for(auto it = eqns.begin(); it != eqns.end();)
    {
	if(size == (*it).s.size() + (*it).t.size())
	{
	    cur_eqns.push_back(*it);
	    it = eqns.erase(it);
	}
	else
	{
	    ++it;
	}
    }
    gurobi_filter_equations(cur_eqns);
    // feasible_eqns.(feasible_eqns.end(), cur_eqns.begin(), cur_eqns.end());

    printf("feasible equations = %lu\n", cur_eqns.size());
    int idx = 0; 
    for(equation& eqn : cur_eqns)
    {
	printf("try to resolve ");
	eqn.print(idx++);
	    
	catfish sc;
	save(sc);
	resolve_equation(eqn);
	load(sc);
	if(eqn.f == 2 && eqn.d == 0)
	{
	    printf("  succeed\n");
	    resolve_equation(eqn);
	    return true;
	}
	else
	{
	    printf("  failed, f = %1d, d = %d\n", eqn.f, eqn.d);
	}
    }

    if(cur_eqns.empty()) return false;
    else
    {
	try
	{
	    // the last test in gurobi_filter_equation may have made the model
	    // infeasible, re-run the model to get a feasible solution to be
	    // used for finding mergeable edges
	    int result = gurobi_model->get(GRB_IntAttr_Status);
	    if(result != GRB_OPTIMAL)
	    {
		gurobi_is_feasible();
	    }
	    return gurobi_find_mergeable_edges();
	}
	catch(GRBException& e)
	{
	    std::cerr << "Gurobi error code = " << e.getErrorCode() << std::endl;
	    std::cerr << e.getMessage() << std::endl;
	    throw;
	}
    }
}

bool catfish::decompose_one_path_from_LP()
{
    assert(gurobi_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL);
    if(output_tex_files) gurobi_print_current_solution();
    printf("Use feasible model to extract one path:\n");
    int ne = gurobi_vars.size();

    assert(i2w.size() == ne);
    std::vector<int> edgeIdxByWeight(ne);
    for(int i = 0; i < ne; ++i)
    {
	edgeIdxByWeight[i] = i;
    }
    std::sort(edgeIdxByWeight.begin(), edgeIdxByWeight.end(),
	      [this](int i1, int i2) { return this->i2w[i1] > this->i2w[i2]; });

    for(int cur : edgeIdxByWeight)
    {
	if(i2e[cur] == null_edge) break; // i2w[cur] should be 0, and all later edges in edgeIdxByWeight are null

	double w = i2w[cur];
	std::vector<int> p;
	bool isPath = true;
	for(int j = 0; j < ne; ++j)
	{
	    if(i2e[j] == null_edge) continue;
	    double sol = gurobi_vars[cur][j].get(GRB_DoubleAttr_X);
	    if(fabs(sol - w) <= SMIN)
	    {
		p.push_back(j);
	    }
	    else if(fabs(sol) <= SMIN) continue;
	    else
	    {
		isPath = false;
		break;
	    }
	}
	if(isPath && p.size() > 1)
	{
	    // sorting by end point does not work as reverse operations may have
	    // disrupted the order of vertices, so we may have a path 0->6->4->19
	    //std::sort(p.begin(), p.end(),
	    //	      [this](int i1, int i2) { return this->i2e[i1]->source() < this->i2e[i2]->source(); });

	    // map from source to index of the edge (in i2e) for all edges in p
	    int nv = gr.num_vertices();
	    std::vector<int> s2i(nv, -1); 
	    for(int x : p)
	    {
		int s = i2e[x]->source();
		assert(s2i[s] == -1);
		s2i[s] = x;
	    }

	    std::vector<int> pp;
	    printf("  split merge path saturating edge %d of weight %.1f:\n", cur, w);

	    int cur_t = 0;
	    int pre_e = -1;
	    while(cur_t != nv - 1)
	    {
		int cur_e = s2i[cur_t];
		assert(cur_e != -1);
		cur_t = i2e[cur_e]->target();
		pp.push_back(cur_e);
		// printf(" %d", cur_e);
		print_edge_with_backtrack(cur_e);
                if(pre_e >= 0 && has_truth && equation_comply_with_truth(
                       equation(std::vector<int>{pre_e}, std::vector<int>{cur_e})))
                {
                    printf("    merge %d and %d OK with ground truth\n", pre_e, cur_e); 
                }
		pre_e = cur_e;
	    }
	    std::vector<int> vv;
	    split_merge_path(pp, w, vv);
	    return true;
	}
    }

    printf("  no path found\n");
    return false;
}
