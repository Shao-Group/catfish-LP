bin_PROGRAMS=catfish

GTFDIR = ../lib/gtf
UTILDIR = ../lib/util
GRAPHDIR = ../lib/graph

catfish_CPPFLAGS = -std=c++11 -I$(GTFDIR) -I$(GRAPHDIR) -I$(UTILDIR) -I$(GUROBI_HOME)/include
catfish_LDFLAGS = -L$(GTFDIR) -L$(GRAPHDIR) -L$(UTILDIR) -L$(GUROBI_HOME)/lib
catfish_LDADD = -lgtf -lgraph -lutil -lgurobi_c++ -lgurobi110

catfish_SOURCES = splice_graph.h splice_graph.cc \
				  sgraph_compare.h sgraph_compare.cc \
				  vertex_info.h vertex_info.cc \
				  edge_info.h edge_info.cc \
				  nested_graph.h nested_graph.cc \
				  interval_map.h interval_map.cc \
				  subsetsum.h subsetsum.cc \
				  config.h config.cc \
				  path.h path.cc \
				  equation.h equation.cc \
				  gtf.h gtf.cc \
				  truth_path.h truth_path.cc \
				  catfish.h catfish.cc \
				  main.cc
