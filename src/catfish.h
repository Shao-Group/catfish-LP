/*
  Part of Catfish
  (c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
  (c) 2025 by Ke Chen @ Penn State.
  See LICENSE for licensing.
*/

#ifndef __SCALLOP1_H__
#define __SCALLOP1_H__

#include "path.h"
#include "equation.h"
#include "splice_graph.h"
#include "nested_graph.h"
#include <unordered_set>
#include <list>
#include "gurobi_c++.h"
#include "truth_path.h"

#define GRB_THREAD_LIMIT 5

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for perfectly estimated splice graph
class catfish
{
public:
    catfish();
    catfish(const string &name, const splice_graph &gr);
    catfish(const string &name, const splice_graph &gr,
	    const std::vector<std::pair<double, std::vector<int>>>& truth);
    virtual ~catfish();

public:
    string name;			// name for this gene
    splice_graph gr;		// splice graph
    /*
      nested_graph nt;		// nested graph
      int gr_version;			// version of splice graph
      int nt_version;			// version of nested graph
    */

    MEI e2i;				// edge map, from edge to index
    VE i2e;					// edge map, from index to edge
    MEV mev;				// super edges
    int round;				// round in iteration

    vector<path> paths;		// predicted transcripts

private:
    // store set of equations for each round,
    // fill by identify_equationsX(), reset by iterate()
    std::vector<equation> eqns;
    void print_edge_with_backtrack(int i);

    // gurobi environment, documentation suggests one env per program (thread)
    // set in the parameterized constructor but not the default constructor
    GRBEnv* gurobi_env = nullptr;
    // a bit ugly to have these class variables, but it makes things easier to
    // reuse the model within the same round for different equantions found
    // by decompose_with_equations(X).
    bool gurobi_base_model_built = false;
    GRBModel* gurobi_model = nullptr;
    std::vector<double> i2w;
    std::vector<GRBVar*> gurobi_vars;
    // used by gurobi_find_mergeable_edges to test if the model remains
    // feasible when requiring an edge to take a value less than its
    // current solution (namely, <= current_solution - gurobi_slack).
    const double gurobi_slack = 0.1;

    /*
     * Variables for backtrack each edge current in gr to a list of edges
     * in the original graph. This helps determine if equations agree
     * with the ground truth and if merging edges conform to the ground truth.
     *
     * Same size and indexing as i2e.
     */
    std::vector<std::vector<edge_descriptor>> i2backtrack;
    // maps original edges in the graph (contents of i2backtrack) to
    // their source and target, does not change after initialization
    // (in assemble0) and does not participate in save/load
    std::map<edge_descriptor, std::pair<int, int>> original_e2st;

    /*
     * Variables and functions for checking if equations/merge operations
     * comply with the provided ground truth paths.
     * The variables do not change (after reconstruct_splice_graph()),
     * and do not participate in save and load.
     */
    bool has_truth;
    std::vector<truth_path> truth;
    // map a path represented as a list of vertices to a set of edge_descriptors
    // according to gr
    void build_truth(const std::vector<std::pair<double, std::vector<int>>>& truth_paths);
    // returns if the given equation complies with the ground truth paths,
    // print violation info if it does not.
    bool equation_comply_with_truth(const equation& eqn);
    // use boost::push_relabel_max_flow to determine if the equation complies
    // with the ground truth paths, lhs and rhs contains edge indices (in i2e)
    // of the two sides of the equation. p2e maps the index of a ground truth
    // path (in truth) to the indices of the two edges (one in lhs, and one
    // in rhs) that it passes through. Note that the edge indices in p2e are
    // indices in lhs and rhs, namely, 3->{0,1} means truth[3] passes through
    // edges lhs[0] and rhs[1].
    int get_max_flow(const std::vector<int>& lhs,
		     const std::vector<int>& rhs,
		     const std::map<int, std::pair<int, int>>& p2e);
    

public:
    int assemble();
    int write(const string &file);
    int write(ofstream &fout);

private:
    int clear();
    int save(catfish &sc);
    int load(catfish &sc);

    int assemble0();
    int assemble1();
    int assemble2();
    int greedy();

    // iteratively decompose
    int iterate();
    bool decompose_trivial_vertices();
    bool decompose_with_equations(int level);

    // trivial, or hard
    int classify();

    // simplify the splice graph and init all data structures
    int init_super_edges();
    int reconstruct_splice_graph();
    bool init_trivial_vertex(int x);

    // remove empty edges, not use now
    int remove_empty_edges();

    // identify and handle equations
    
    // adjacent equal edges, can always resolve, return true if found (f=2)
    bool identify_equations0();

    // distant equal edges, return true if
    // 1. found one that can be resolved by reverse operations (f=2)
    // 2. in fast mode and found one that can resolve_vertex_with_equation (f=3)
    // otherwise return false, in this case, eqns contains all distant equal
    // edges that can be connected by a path with enough weight (f=1)
    //
    // If used for inferring independent edges, do not need to check with LP.
    // (current implementation still checks).
    bool identify_equations1();

    // x = a + b + ..., return true if
    // 1. found one fully resolvable equation (f=2 && d=0)
    // 2. in fast mode and found one that can resolve_vertex_with_equation (f=3)
    // otherwise return false, in this case, eqns contains all such equations
    bool identify_equations2();

    // x + y = a + b + ..., return true if
    // 1. found one fully resolvable equation (f=2 && d=0)
    // 2. in fast mode and found one that can resolve_vertex_with_equation (f=3)
    // otherwise return false, in this case, eqns contains all such equations
    bool identify_equations3();
    bool verify_equation_nontrivial(equation &eqn);
    bool verify_equation_mergable(equation &eqn);

    // resolve equation
    int resolve_equation(equation &eqn);
    int resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md);

    // use equation to decompose trivial vertex, new feature
    bool resolve_vertex_with_equations(vector<equation> &eqns);
    bool resolve_vertex_with_equation(equation &eqn);
    bool resolve_vertex_with_equation1(equation &eqn);
    bool resolve_vertex_with_equation2(equation &eqn);

    // split, and merge
    int split_edge(int exi, double w);
    int merge_adjacent_equal_edges(int x, int y);
    int split_merge_path(const vector<int> &p, double ww, vector<int> &vv);
    int split_merge_path(const VE &p, double ww, vector<int> &vv);
    int merge_adjacent_edges(int x, int y);

    // check, and make two edges adjacent
    bool check_adjacent_mergable(int ex, int ey, vector<PI> &p);
    bool check_adjacent_mergable(int ex, int ey, nested_graph &nt);
    int check_distant_mergable(int x, int y, double w, VE &p);
    int check_distant_mergable(int x, int y, double w);
    int build_adjacent_edges(const vector<PI> &p);

    // decompose the graph with greedy algorithm
    int greedy_decompose();

    // collect existing s-t path e
    int collect_path(int e);
    int collect_existing_st_paths();

    // test, print and draw
    int draw_splice_graph(const string &file);
    int print();

    /**********************************************************************
     *** functions for checking if equations are feasible using LP
     **********************************************************************
     * Given a splice graph and a set of equations, the equations are said to be
     * feasible if there exists a flow decompostion of the graph that respects 
     * all the equations. We check a necessary condition using LP as follows:
     *
     * Base model:
     * Let m be the number of edges in the splice graph, the base model contains
     * m^2 continuous variables. For each edge e, a set of m variables v_e(.) are
     * defined, one for each edge in the graph; together, they represent a valid
     * sub-st-flow that saturates the edge e. By symmetry, v_a(b) should equal
     * v_b(a) for all pairs of edges (a, b).
     * If e has parallel edges e', the resulting flow v_e(.) may not all go
     * through e, so we additionally set v_e(e')==0.
     *
     * Additional constraints model an equation:
     * Given an equation a+b=x+y+z (actually, w(a)+w(b)=w(x)+w(y)+w(z)), we add
     * 2m constraints, two for each edge e which are 
     * v_a(e)+v_b(e)==v_x(e)+v_y(e)+v_z(e)<=w(e). Note for example that 
     * v_a(a)+v_b(a)<=w(a) (together with the base requirement that v_a(a)==w(a))
     * implies that v_b(a)==0, i.e., the flows through a and those through b are
     * independent. 
     *
     * Notes:
     * 1. we could explicitly use the same variable for v_a(b) and v_b(a) which
     *    gives m(m+1)/2 variables instead of m^2, the current implementation
     *    leaves this to the solver.
     * 2. for a single equation, the additional constraints hold if and only if
     *    the equation is feasible; on the other hand, this only models a necessary
     *    condition for a set of equations, for example, two equations may both
     *    (independently) require w(e) flow on edge e so they cannot be feasible
     *    at the same time but it does not violate any constraints described above.
     */

    // Entry point of checking equations and build an LP model with a set
    // of feasible equations. Then call gurobi_find_mergeable_edges().
    // Return true if graph modified.
    bool decompose_with_equations_induced_merge();

    // Assume gurobi_model has been initialized.
    // All variables are stored in gurobi_vars, one entry per edge
    // each points to a list of variables with the same indices as i2w
    void gurobi_setup_base_model();
    
    // Assume gurobi_setup_base_model() has been called.
    // Add new constraints according to the given equation
    // return true if the model remains feasible
    // otherwise the new constraints are removed and return false
    //
    // Note that gurobi uses lazy approach, updates are queued but
    // not applied until optimize() or update() is called.
    bool gurobi_add_equation(equation& eqn);
    
    // Optimize the model, return true if gurobi_model is feasible
    bool gurobi_is_feasible();

    // Prepare i2w used by above functions
    void gurobi_get_edge_weights();

    // Clear flags and free gurobi variables
    void gurobi_clear_model();

    // Assume an LP model with a set of feasible equations has been solved.
    // For each edge e and the LP solution v_e(.) that saturates e, check
    // if any edge x adjacent to e has v_e(x)==w(e); if found, set a
    // temporary constraint that v_e(x)<w(e) and resolve; if the model
    // becomes infeasible, we know e should be merged with x (assuming
    // the feasible equations in the model are all correct).
    bool gurobi_find_mergeable_edges();

    // Assume gurobi_is_feasible has been called after the last change
    // of the model, otherwise the solution may be from an earilier
    // version of the model.
    void gurobi_print_current_solution();

    // Add each equation in eqns one by one to the current model and
    // check for feasibility. Modifies eqns to only keep feasible equations.
    void gurobi_filter_equations(std::list<equation>& eqns);
    
    /**********************************************************************
     *** end of functions for checking if equations are feasible using LP
     *********************************************************************/

    // call eqn.print(i) and also show backtrack info for each edge
    void print_equation(const equation& eqn, int idx);

    /**********************************************************************
     *** functions for the new pipeline: classify equations by size, check
     *** by LP, try resolve, try LP merge if none can be resolved
     *********************************************************************/

    // iterate through each vertex and merge adjacent equal edges, return
    // true if any merge happened
    bool decompose_adjacent_equal_edges();

    // similar to identify_equations2/3, but without the heuristics/try
    // to resolve part
    void find_equations_lhs1(std::list<equation>& results);
    void find_equations_lhs2(std::list<equation>& results);

    // entry point for the new pipeline, called in iterate(), return true
    // if any modification to the graph happened
    bool decompose_with_equations();
    // handles equations of the given size, all such equations will be
    // removed from eqns, return true if any modification to the graph
    // according to the selected equations happened
    bool filter_and_try_resolve_equations(int size, std::list<equation>& eqns);

    // if all strategies failed, try to see if the latest LP solution for
    // any edge is a single path, if found, extract the heaviest path and
    // return true
    bool decompose_one_path_from_LP();
    
    /**********************************************************************
     *** end of functions for the new pipeline
     *********************************************************************/
};

#endif
