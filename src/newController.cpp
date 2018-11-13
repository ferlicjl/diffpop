/*
 * =====================================================================================
 *
 *       Filename:  newController.cpp
 *
 *    Description:  Contains simulating and input file parsing functions
 *
 *        Version:  1.0
 *        Created:  10/1/2017 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (s), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */

// Cell Population classes
#include "FixedPop.h"
#include "DiffTree.h"
#include "GrowingPop.h"
#include "DiffTriangle.h"
#include "NodeList.h"

// Random Distributions and Helper functions
#include "constantRVFunctions.h"
#include "helpers.h"

// Includes
#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <iterator>
#include <chrono>
#include <cstring>
#include <cmath>
#include <limits>
#include <gsl/gsl_randist.h>

// Rcpp headers
#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

// GSL random number generators
gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
double seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();


// Helper Methods to read in data from text files
template
< typename T
  , template<typename ELEM, typename ALLOC=std::allocator<ELEM> > class Container
  >
std::ostream& operator<< (std::ostream& o, const Container<T>& container)
{
    typename Container<T>::const_iterator beg = container.begin();

    o << "["; // 1

    while(beg != container.end())
    {
        o << " " << *beg++; // 2
    }

    o << " ]"; // 3

    return o;
}

// trim from left
static inline std::string &ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s)
{
    return ltrim(rtrim(s));
}


// Read in a vector matrix
std::vector<std::vector<double> > fileToVectorMatrix(std::string name)
{
    std::vector<std::vector<double> > result;
    std::ifstream input (name);
    std::string lineData;

    while(getline(input, lineData))
    {
        double d;
        std::vector<double> row;
        std::stringstream lineStream(lineData);

        while (lineStream >> d)
            row.push_back(d);

        result.push_back(row);
    }

    return result;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T> >& v)
{
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

// Read in a vector
std::vector<double> fileToVector(std::string name)
{
    std::vector<std::vector<double> > result;
    std::ifstream input (name);
    std::string lineData;

    while(getline(input, lineData))
    {
        double d;
        std::vector<double> row;
        std::stringstream lineStream(lineData);

        while (lineStream >> d)
            row.push_back(d);

        result.push_back(row);
    }

    return flatten(result);
}

// Read in a text file, output a vector of strings, one string per line, trimmed
std::vector<std::string> inputStringVector(std::string in)
{
    std::string s;
    std::vector<std::string> a;
    std::ifstream infile(in);
    if(infile.is_open())
    {
        while(getline(infile, s))
        {
            // allow comments
            std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
            s.erase(end_pos, s.end());
            a.push_back(trim(s));
        }
    }

    return(a);
}

std::vector<std::vector<std::string>> funct(std::string filename){
	std::ifstream csv(filename);
	std::string line;
	std::vector <std::vector<std::string>> items;

	if (csv.is_open()) {
			for (std::string row_line; std::getline(csv, row_line);)
			{
				items.emplace_back();
				std::istringstream row_stream(row_line);
				for(std::string column; std::getline(row_stream, column, ',');)
					items.back().push_back(column);
			}
	}
	else {
		std::cout << "Unable to open file";
	}

	return items;
}

NodeList readFile(std::string indir, std::string fname){
	NodeList n;
	std::vector<std::vector<std::string>> l = funct(indir + fname);
	//std::cout << "Inserting " << l.size() << " rows" << std::endl;
	for(size_t i = 0; i < l.size(); i++){
		long int b = stol(l[i][0]);
		std::string m = l[i][1];
		double f = stod(l[i][2]);
		int c = stoi(l[i][3]);
		
		//std::cout << "Inserting: " << std::to_string(b) << " " << m << " "<< std::to_string(f) << " "<< std::to_string(c) << std::endl;
		n.insert(b, m, f, c);
	}
	return n;
}

//' readFile
//'
//' Runs simulation of a tree containing Fixed Population Size objects
//' [FixedPops & DiffTriangles].
//'
//'@param indir input directory for tree files
//'@param filename output directory
//'@return Nodelist nodelist represented by file
//'
//' @export
// [[Rcpp::export]]
void readFileR(std::string indir, std::string fname){
	NodeList n = readFile(indir, fname);
	n.traversePrint();
	
	return;
}

// Reads input files generated by the R code, passes this information onto the Differentiation Tree
void readInputFiles(std::string indir, std::map<std::string, GrowingPop*> &m, DiffTree &tree){

	//std::cout << "in readIntputFiles" << std::endl;
	std::vector <std::vector<std::string>> nodes = funct(indir + "nodes.csv");

	//std::ifstream  src(indir + "bfs.txt", std::ios::binary);
    //std::ofstream  dst(tree.opath + ".bfs",   std::ios::binary);

    //dst << src.rdbuf();
	
	long int offset = 0;

    // Create Populations
    for(size_t i = 0; i < nodes.size(); i++)
    {
        double l = stod(nodes[i][4]);
        if (nodes[i][1] == "CellPopulation")
        {
            long int s = stol(nodes[i][3]);
            m.insert(std::make_pair(nodes[i][0], new GrowingPop(nodes[i][0], s, offset, l)));
			offset += s;
        }
        else if (nodes[i][1] == "FixedPopCell")
        {
            long int s = stol(nodes[i][3]);
            m.insert(std::make_pair(nodes[i][0], new FixedPop(nodes[i][0], s, offset, l)));
			offset += s;
        }
        else if(nodes[i][1] == "DiffTriangle")
        {
            int h = stoi(nodes[i][5]);
            int f = stoi(nodes[i][6]);
            m.insert(std::make_pair(nodes[i][0], new DiffTriangle(nodes[i][0], h, f) ));
        }

		// Check if this population is the root
        if (nodes[i][2] == "TRUE")
        {
            tree.setRoot(m[nodes[i][0]]);
        }
    }

	std::vector <std::vector<std::string>> edges = funct(indir + "edges.csv");

    // Iterate over all transitions
    for(size_t i = 0; i < edges.size(); i++)
    {
        double r = stod(edges[i][3]);
        if (edges[i][2] == "alpha")
        {
            m[edges[i][0]]->setAlpha(r);
        }
        else if (edges[i][2] == "delta" || edges[i][2] == "death")
        {
            m[edges[i][0]]->setDeath(r);
        }
        else if (edges[i][2] == "gamma1")
        {
            m[edges[i][0]]->setDiff1(m[edges[i][1]], r);
        }
		else if (edges[i][2] == "gamma2")
        {
            m[edges[i][0]]->setDiff2(m[edges[i][1]], r);
        }
		else if (edges[i][2] == "beta")
        {
            m[edges[i][0]]->setAsymdiff(m[edges[i][1]], r);
        }
        else if(edges[i][2] == "dediff" || edges[i][2] == "zeta")
        {
            m[edges[i][0]]->addDediff(m[edges[i][1]], r);
        }
        else if(edges[i][2] == "mutation" || edges[i][2] == "mu")
        {
            m[edges[i][0]]->setMu(r);
        }
    }

	// Fitness Distribution
	std::vector<std::string> fp;
    fp = inputStringVector(indir + "fitnessdist.txt");
	if(fp.size() > 0){
		FitnessParameters f;
		f.fitness_distribution = fp[0];
		f.is_randfitness = true;
		f.alpha_fitness = stod(fp[1]);
		f.beta_fitness = stod(fp[2]);
		f.pass_prob = stod(fp[3]);

		if(fp[4] != "NA"){
			f.upper_fitness = stod(fp[4]);
		} else {
			f.upper_fitness = std::numeric_limits<double>::max();
		}

		if(fp[5] != "NA"){
			f.lower_fitness = stod(fp[5]);
		} else {
			f.lower_fitness = std::numeric_limits<double>::lowest();
		}
		tree.setFitnessDist(f);
	}
	
	// Read in BFS order of populations
	tree.bfs = inputStringVector(indir + "bfs.txt");
}


//' simulateFixedTreeCodeNew
//'
//' Runs simulation of a tree containing Fixed Population Size objects
//' [FixedPops & DiffTriangles].
//'
//'@param nObs number of observations for simulation time
//'@param traverseFrequency how often to output entire system population to filesystem
//'@param indir input directory for tree files
//'@param outdir output directory
//'@param seed numeric seed for random number generator
//'
//' @export
// [[Rcpp::export]]
int simulateFixedTreeCodeNew(int nObs = 10, int traverseFrequency = -1, std::string indir = "./", std::string outdir = "./", SEXP seed = R_NilValue)
{
	// Set the seed for the GSL random number generator
	 double seedcpp;
	if(Rf_isNull(seed)){
		seedcpp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	} else {
		seedcpp = Rf_asReal(seed);
	}
    gsl_rng_set(rng, seedcpp);

	// Create a DiffTree and give it the output path
    DiffTree tree(outdir);

	// Create a map to tie string population names to their actual object instances
    std::map<std::string, GrowingPop*> m;

	// Read input files
	readInputFiles(indir, m, tree);
	std::cout << "Read input files..." << std::endl;

	// Give the tree the map
	tree.setMap(m);

	// Initialize the tree, instantiates the NodeLIst data structures
	tree.initializeTree();
	std::cout << "Initialized the tree..." << std::endl;

	// Output the initial parameters given by the user
	std::ofstream p_file;
    p_file.open(tree.opath+"_params.csv");

	p_file << "nObs," << nObs << std::endl;
	p_file << "census," << traverseFrequency << std::endl;
	p_file << "indir," << indir << std::endl;
	p_file << "outdir," << outdir << std::endl;
	p_file << "seed," << seed << std::endl;
	p_file << "fp.distribution," << tree.fp.fitness_distribution << std::endl;
	p_file << "fp.alpha_fitness," << tree.fp.alpha_fitness << std::endl;
	p_file << "fp.beta_fitness," << tree.fp.beta_fitness << std::endl;
	p_file << "fp.pass_prob," << tree.fp.pass_prob << std::endl;
	p_file << "fp.lower_fitness," << tree.fp.lower_fitness << std::endl;
	p_file << "fp.upper_fitness," << tree.fp.upper_fitness << std::endl;

    for(auto const& x: tree.bfs)
    {
        p_file << "N_" << m[x]->name << "," << m[x]->cells.total << std::endl;
        p_file << "alpha_" << m[x]->name << ","<< m[x]->alpha << std::endl;
        p_file << "beta_" << m[x]->name << ","<< std::accumulate(m[x]->beta.begin(), m[x]->beta.end(), 0.0)<< std::endl;
        p_file << "gamma1_" << m[x]->name << ","<< std::accumulate(m[x]->gamma1.begin(), m[x]->gamma1.end(), 0.0) << std::endl;
		p_file << "gamma2_" << m[x]->name << ","<< std::accumulate(m[x]->gamma2.begin(), m[x]->gamma2.end(), 0.0) << std::endl;
		p_file << "zeta2_" << m[x]->name << ","<< std::accumulate(m[x]->zeta.begin(), m[x]->zeta.end(), 0.0) << std::endl;
        p_file << "delta_" << m[x]->name << ","<< m[x]->death << std::endl;
		p_file << "mu_" << m[x]->name << ","<< m[x]->mu << std::endl;
    }
    p_file << std::endl;
	p_file.close();

	// Calculate and adjust net proliferation for each population
	tree.calcDelta();

	// Re-output parameters after having adjusted net-proliferation
    p_file.open(tree.opath+"_params.csv", std::fstream::in | std::fstream::out | std::fstream::app);

	p_file << "# Rates after adjusting net proliferation for FixedPops" << std::endl;

    for(auto const& x: tree.bfs)
    {
        p_file << "N_" << m[x]->name << "_adj," << m[x]->cells.total << std::endl;
        p_file << "alpha_" << m[x]->name << "_adj,"<< m[x]->alpha << std::endl;
        p_file << "beta_" << m[x]->name << "_adj,"<< std::accumulate(m[x]->beta.begin(), m[x]->beta.end(), 0.0)<< std::endl;
        p_file << "gamma1_" << m[x]->name << "_adj,"<< std::accumulate(m[x]->gamma1.begin(), m[x]->gamma1.end(), 0.0) << std::endl;
		p_file << "gamma2_" << m[x]->name << "_adj,"<< std::accumulate(m[x]->gamma2.begin(), m[x]->gamma2.end(), 0.0) << std::endl;
		p_file << "zeta2_" << m[x]->name << "_adj,"<< std::accumulate(m[x]->zeta.begin(), m[x]->zeta.end(), 0.0) << std::endl;
        p_file << "delta_" << m[x]->name << "_adj,"<< m[x]->death << std::endl;
		p_file << "mu_" << m[x]->name << "_adj,"<< m[x]->mu << std::endl;
    }
    p_file << std::endl;
	p_file.close();

	std::cout << "Created parameters file..." << std::endl;

	// Print out the tree structure
    tree.print();

	// Simulate the hierarchy
    tree.time_steps(nObs, traverseFrequency);
	
	std::cout << "Completed simulation." << std::endl;

    // Ensure NodeLists get deleted
    for ( auto current = m.begin(); current != m.end(); ++ current )
    {
        current->second->cells.deleteList();
    }
    m.clear();
    return 0;
}

//' simulateTreeCodeNew
//'
//' Runs simulation of a tree structure using Gillespie algorithm
//'
//'@param nObs number of observations for simulation time
//'@param traverseFrequency how often to output entire system population to filesystem
//'@param indir input directory for tree files
//'@param outdir output directory
//'@param seed numeric seed for random number generator
//'
//' @export
// [[Rcpp::export]]
int simulateTreeCodeNew(int nObs = 10, int traverseFrequency = -1, std::string indir = "./", std::string outdir = "./", SEXP seed = R_NilValue)
{
    int t=time(NULL);
    srand(t);

    double seedcpp;
	if(Rf_isNull(seed)){
		seedcpp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	} else {
		seedcpp = Rf_asReal(seed);
	}
    //seed = 1;
    gsl_rng_set(rng, seedcpp);


    DiffTree tree(outdir);

    std::map<std::string, GrowingPop*> m;

	readInputFiles(indir, m, tree);
	std::cout << "Read input files..." << std::endl;
	tree.setMap(m);

	tree.initializeTree();
	std::cout << "Initialized the tree..." << std::endl;

	// Parameters File
	std::ofstream p_file;
    p_file.open(tree.opath+"_params.csv");

	p_file << "nObs," << nObs << std::endl;
	p_file << "census," << traverseFrequency << std::endl;
	p_file << "indir," << indir << std::endl;
	p_file << "outdir," << outdir << std::endl;
	p_file << "seed," << seed << std::endl;
	p_file << "fp.distribution," << tree.fp.fitness_distribution << std::endl;
	p_file << "fp.alpha_fitness," << tree.fp.alpha_fitness << std::endl;
	p_file << "fp.beta_fitness," << tree.fp.beta_fitness << std::endl;
	p_file << "fp.pass_prob," << tree.fp.pass_prob << std::endl;
	p_file << "fp.lower_fitness," << tree.fp.lower_fitness << std::endl;
	p_file << "fp.upper_fitness," << tree.fp.upper_fitness << std::endl;	

    for(auto const& x: tree.bfs)
    {
        p_file << "N_" << m[x]->name << "," << m[x]->cells.total << std::endl;
        p_file << "alpha_" << m[x]->name << ","<< m[x]->alpha << std::endl;
        p_file << "beta_" << m[x]->name << ","<< std::accumulate(m[x]->beta.begin(), m[x]->beta.end(), 0.0)<< std::endl;
        p_file << "gamma1_" << m[x]->name << ","<< std::accumulate(m[x]->gamma1.begin(), m[x]->gamma1.end(), 0.0) << std::endl;
		p_file << "gamma2_" << m[x]->name << ","<< std::accumulate(m[x]->gamma2.begin(), m[x]->gamma2.end(), 0.0) << std::endl;
		p_file << "zeta2_" << m[x]->name << ","<< std::accumulate(m[x]->zeta.begin(), m[x]->zeta.end(), 0.0) << std::endl;
        p_file << "delta_" << m[x]->name << ","<< m[x]->death << std::endl;
		p_file << "mu_" << m[x]->name << ","<< m[x]->mu << std::endl;
    }
    p_file << std::endl;
	p_file.close();

	std::cout << "Created parameters file..." << std::endl;
    tree.print();

    tree.simulate(nObs, traverseFrequency);

	std::cout << "Completed simulation." << std::endl;

    // Make sure NodeLists get deleted
    for ( auto current = m.begin(); current != m.end(); ++ current )
    {
		current->second->cells.deleteList();
	}
    m.clear();
    return 0;
}


