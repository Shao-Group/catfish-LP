/*
Part of Catfish
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "config.h"
#include "catfish.h"
#include "gtf.h"
#include "genome.h"
#include "nested_graph.h"
#include "sgraph_compare.h"
#include "boost/optional.hpp"

using namespace std;

int assemble(const string &file);
int transform(const string &file);
vector<vector<pair<double, vector<int>>>> parse_ground_truth_file(const string& filename);


int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		print_copyright();
		print_help();
		return 0;
	}

	srand(time(0));
	parse_arguments(argc, argv);
	//print_parameters();

	if(input_file == "" && output_file != "" && algo == "simulate")
	{
		splice_graph sg;
		sg.simulate(simulation_num_vertices, simulation_num_edges, simulation_max_edge_weight);
		sg.write(output_file);
		return 0;
	}

	if(input_file != "" && output_file != "" && algo == "transform")
	{
		string s = input_file.substr(input_file.size() - 3, 3);
		if(s != "gtf") return 0;
		transform(input_file);
		return 0;
	}

	if(input_file != "")
	{
		assemble(input_file);
	}
	return 0;
}

int assemble(const string &file)
{
	ifstream fin(file);
	if(fin.fail()) return 0;

	ofstream fout;
	if(output_file != "") fout.open(output_file);

	char buf[10240];
	string header;
	stringstream sstr;

	vector<vector<pair<double, vector<int>>>> truth_paths;
	if(truth_file != "")
	{
	    truth_paths = parse_ground_truth_file(truth_file);
	}

	int ct = 0;
	while(fin.getline(buf, 10240, '\n'))
	{
		string s(buf);
		if(s.size() == 0) continue;
		if(s[0] == '#')
		{
			if(header != "")
			{
				splice_graph sg;
				sg.build(sstr);

				optional<catfish> sc;
				
				if(truth_file != "")
				{
				    if(ct >= truth_paths.size())
				    {
					cerr << "truth file only has " << truth_paths.size()
					     << " sets of paths, but " << header << " is graph #"
					     << ct+1 << " in the input file." << endl;
					sc.emplace(header.substr(1), sg);
				    }
				    else
				    {
					sc.emplace(header.substr(1), sg, truth_paths[ct]);
					++ct;
				    }
				}
				else
				{
				    sc.emplace(header.substr(1), sg);
				}
				
				sc->assemble();
				if(output_file != "") fout << header.c_str() << " paths = " << sc->paths.size() << endl;
				if(output_file != "") sc->write(fout);
			}

			sstr.clear();
			header = s;
		}
		else
		{
			sstr << s << "\n";
		}
	}

	if(header != "")
	{
		splice_graph sg;
		sg.build(sstr);

		optional<catfish> sc;
		if(truth_file != "")
		{
		    if(ct >= truth_paths.size())
		    {
			cerr << "truth file only has " << truth_paths.size()
			     << " sets of paths, but " << header << " is graph #"
			     << ct+1 << " in the input file." << endl;
			sc.emplace(header.substr(1), sg);
		    }
		    else
		    {
			sc.emplace(header.substr(1), sg, truth_paths[ct]);
		    }
		}
		else
		{
		    sc.emplace(header.substr(1), sg);
		}

		sc->assemble();
		if(output_file != "") fout << header.c_str() << " paths = " << sc->paths.size() << endl;
		if(output_file != "") sc->write(fout);
	}

	fin.close();
	if(output_file != "") fout.close();

	return 0;
}

int transform(const string &file)
{
	genome g(file);

	ofstream fout1(output_file + ".graph");
	ofstream fout2(output_file + ".truth");

	for(int i = 0; i < g.genes.size(); i++)
	{
		gtf gg(g.genes[i]);

		string name = gg.get_gene_id();

		fout1 << "# graph number = " << i << " name = " << name.c_str() << endl;
		fout2 << "# graph number = " << i << " name = " << name.c_str() << endl;

		splice_graph gr;
		gg.build_splice_graph(gr);

		gr.write(fout1);
		gg.write_transcript_paths(fout2);
	}

	fout1.close();
	fout2.close();
	return 0;
}

vector<vector<pair<double, vector<int>>>> parse_ground_truth_file(const string& filename)
{
    ifstream fin(filename);

    vector<vector<pair<double, vector<int>>>> result;
    vector<pair<double, vector<int>>> current;
    string line;

    while(getline(fin, line))
    {
        if(line.empty()) continue;

        if(line[0] == '#')
	{
            if(!current.empty())
	    {
                result.push_back(current);
                current.clear();
            }
        }
	else
	{
            stringstream ss(line);
            double w;
            ss >> w;

            vector<int> path;
            int v;
            while(ss >> v)
	    {
                path.push_back(v);
            }

            current.emplace_back(w, path);
        }
    }

    if (!current.empty()) {
        result.push_back(current);
    }

    fin.close();
    return result;
}
