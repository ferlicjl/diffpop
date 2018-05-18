/*
 * =====================================================================================
 *
 *       Filename:  DiffTriangle.cpp
 *
 *    Description:  Class for differentiation triangles
 *
 *        Version:  1.0
 *        Created:  03/14/2017 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * ============================================= ========================================
 */

#include "DiffTriangle.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <numeric>
#include <algorithm>
#include <assert.h>     /* assert */
#include <map>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iterator>
#include <deque>
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;

// All Constructors simply call CellType constructors
DiffTriangle::DiffTriangle() : GrowingPop()
{
    num_levels = 0;
    mfactor = 0;
    updated = false;
}

DiffTriangle::DiffTriangle(std::string cname) : GrowingPop(cname)
{

}

DiffTriangle::DiffTriangle(std::string cname, int nlevels, int mfac) : GrowingPop(cname)
{
    num_levels = nlevels;
    mfactor = mfac;
    numCells = (std::pow(2, nlevels) - 1) * mfac;
    updated = false;

	/*
    //tri_ind = std::vector<std::vector<long int>>(mfac);
    long int ind = 0;
    for(int i = 0; i < mfac; i++)
    {
        for(int j = 1; j <= nlevels; j++)
        {
			// Insert Barcode, M, F, triangle, level, count
            cells.insert(-1, "", 1, i, j, std::pow(2, j-1));
        }
    }
	*/
	ptrs = cells.makeTriangle(nlevels, mfac);
}

/*
DiffTriangle::DiffTriangle(std::string cname, int nlevels, int mfac, long int offset) : GrowingPop(cname){
	// Calculate/set number of cells, levels, and mfactor
	numCells = (std::pow(2, nlevels) - 1) * mfac;
	num_levels = nlevels;
	mfactor = mfac;

	// Push back a unique code for each cell
	for(long int i = 0 + offset; i < numCells + offset; i++){
		codes.push_back(i);
	}

	// Set up the indices vector structure
	// Used to add cells to triangles
	indices = std::vector<std::vector<long int>>(nlevels);

	// Use k as index number
	int k = 0;

	// Iterate over each level i
	for(int i = 0; i < nlevels; i++){
		// Calculate j, the number of cells in the ith level
		for(int j = 0; j < std::pow(2, i) * mfac; j++){
			// Push back index k to the ith index vector
			indices[i].push_back(k);
			k++;
		}
	}
}
*/

/*
DiffTriangle::DiffTriangle(std::string cname, int nlevels, int mfac, long int offset, double eff) : GrowingPop(cname){
	// Calculate/set number of cells, levels, and mfactor
	numCells = (std::pow(2, nlevels) - 1) * mfac;
	num_levels = nlevels;
	mfactor = mfac;
	efficiency = eff;

	// Push back a unique code for each cell
	for(long int i = 0 + offset; i < numCells + offset; i++){
		//codes.push_back(i);
		int labelled = gsl_ran_bernoulli(rng, efficiency);
		if(labelled > 0){
			codes.push_back(true);
		} else {
			codes.push_back(false);
		}
	}

	// Set up the indices vector structure
	// Used to add cells to triangles
	indices = std::vector<std::vector<long int>>(nlevels);

	// Use k as index number
	int k = 0;

	// Iterate over each level i
	for(int i = 0; i < nlevels; i++){
		// Calculate j, the number of cells in the ith level
		for(int j = 0; j < std::pow(2, i) * mfac; j++){
			// Push back index k to the ith index vector
			indices[i].push_back(k);
			k++;
		}
	}
}
*/

// Destructor - no dynamic memory
DiffTriangle::~DiffTriangle()
{
}

/*
std::vector<DiffTriangle*> DiffTriangle::getChildren(){
	return triChildren;
}
*/

// Add a child GrowingPop to current GrowingPop
// initiates a 0 differentiation rate
// sets the child's parent to this GrowingPop
void DiffTriangle::addTriChild(DiffTriangle* child)
{
    child->setParent(this);
    triChildren.push_back(child);
    gamma1.push_back(0);
    numChildren++;
}

// Add a child GrowingPop with given differentiation rate
// sets the child's parent to this GrowingPop
void DiffTriangle::addTriChild(DiffTriangle* child, double diff_rate)
{
    child->setParent(this);
    triChildren.push_back(child);
    gamma1.push_back(diff_rate);
    numChildren++;
}

void DiffTriangle::initialize(long int offset)
{
    //offset = 0;
    //std::cout << "DiffTriangle Initialize: " << offset << std::endl;
    /*
    for(long int i = 0; i < numCells; i++){
    	if(codes[i] >= 0){
    		codes[i] += offset;
    	}
    }
    */
    return;
}


/*
// Enacts the events that would occur in one-unit of time according to rates determined by the cell type
// Importantly, events are ordered such that cell population size remains constant
void DiffTriangle::time_step(bool verbose){
	// Create vector of 3 probabilities:
	// 1. split(s), 2. die, 3. differentiate
	double div_rates = alpha + beta + std::accumulate(gamma1.begin(), gamma1.end(), 0.0);

	// Select number of events from these rates
	int n_events = codes.size() * div_rates;
	if(verbose)
		std::cout << n_events << std::endl;

	// Set up a vector to be able to choose an event type
	std::vector<double> init_rates;
	init_rates.push_back(alpha);
	init_rates.push_back(beta);

	std::vector<double> down_rates;
	down_rates.push_back(death);
	down_rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));

	// Set up simple string to print out event type name
	std::vector<std::string> init_names;
	init_names.push_back("Mitosis");
	init_names.push_back("Asym. Differentiation");

	std::vector<std::string> down_names;
		down_names.push_back("Apoptosis");
	for(int j = 0; j < numChildren; j++){
		down_names.push_back(children[j]->name);
	}

	// Iterate over the number of events
	for(int i = 0; i < n_events; i++){
		if(verbose)
			std::cout << "Step " << (i+1) << ": ";
		int choice = choose(init_rates);
		if(verbose)
			std::cout << init_names[choice] << ' ';

		// Gained a cell, need to remove by initiating a down event
		if(choice == 0){
			mitosis();
			int choice2 = choose(down_rates);
			// Cell death, or a differentiation event
			if(choice2 == 0)
				apoptosis();
			else{
				diff1();
				i++;
			}
		}
		// Asym. Diff doesn't change cell population size
		else
			asym_diff();
	}
	if(verbose)
		std::cout << std::endl;
}

*/

/*
// Enacts the events that would occur in one-unit of time according to rates determined by the cell type
// Importantly, events are ordered such that cell population size remains constant
void DiffTriangle::single_event(bool verbose){
	// Create vector of 3 probabilities:
	// 1. split(s), 2. die, 3. differentiate
	double div_rates = alpha + beta + std::accumulate(gamma1.begin(), gamma1.end(), 0.0);

	// Select number of events from these rates
	int n_events = codes.size() * div_rates;
	n_events = 1;
	if(verbose)
		std::cout << n_events << std::endl;

	// Set up a vector to be able to choose an event type
	std::vector<double> init_rates;
	init_rates.push_back(alpha);
	init_rates.push_back(beta);

	std::vector<double> down_rates;
	down_rates.push_back(death);
	down_rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));

	// Set up simple string to print out event type name
	std::vector<std::string> init_names;
	init_names.push_back("Mitosis");
	init_names.push_back("Asym. Differentiation");

	std::vector<std::string> down_names;
		down_names.push_back("Apoptosis");
	for(int j = 0; j < numChildren; j++){
		down_names.push_back(children[j]->name);
	}

	// Iterate over the number of events
	for(int i = 0; i < n_events; i++){
		if(verbose)
			std::cout << "Step " << (i+1) << ": ";
		int choice = choose(init_rates);
		if(verbose)
			std::cout << init_names[choice] << ' ';

		// Gained a cell, need to remove by initiating a down event
		if(choice == 0){
			mitosis();
			int choice2 = choose(down_rates);
			// Cell death, or a differentiation event
			if(choice2 == 0)
				apoptosis();
			else{
				diff1();
				i++;
			}
		}
		// Asym. Diff doesn't change cell population size
		else
			asym_diff();
	}
	if(verbose)
		std::cout << std::endl;
}

*/

void DiffTriangle::add_cell(long int code, std::string m, double f)
{
	addcellEvents += 1;
	/*
    updated = false;
    addcellEvents += 1;

    // Set up data structures to hold our codes to add and our placement indices
    std::deque<long int> toAdd;
    std::deque<long int> pIndices;

    std::vector<double> probs(mfactor, (double)1 / mfactor);
    int choice = choose(probs);


    //std::cout << "Choice " << name << " " << mfactor << ": " << choice << std::endl;
    //triangles[choice].add_cell(code);

    // Get cells from last row we need to pass on or die
    int lastLevelFirst = std::pow(2, num_levels - 1) - 1;
    int lastLevelLast = std::pow(2, num_levels) - 2;
    std::vector<long int> leftovers(codes.begin() + tri_ind[choice][lastLevelFirst], codes.begin() + tri_ind[choice][lastLevelLast]);

    // Get the index of the last cell in pyramid
    int lastIndex = std::pow(2, num_levels - 1) - 2;

    // Iterate from the bottom of the pyramid, setting a child node value to its parent's value
    for(int i = lastIndex; i > 0; i--)
    {
        codes[tri_ind[choice][i]] = codes[tri_ind[choice][(i-1)/2]];
    }

    // Set the first cell to our new value
    codes[tri_ind[choice][0]] = code;
	*/
	int tri = gsl_rng_uniform_int(rng, mfactor);
	
	/*
	OLD WAY!!!
	//std::cout << "Send to triangle " << tri << std::endl;
	NodeList matured = removeLastAndDouble(tri);
	//std::cout << "Here's the matured cells!" << std::endl;
	//matured.traversePrint();
	//doubleCells(tri);
	//std::cout << "I've doubled my cells!" << std::endl;
	cells.insert(code, m, f, tri, 1, 1);
	//std::cout << "I've added my new cell!" << std::endl;
	*/
	
	
	//NEW WAY
	
	//Save off information from last level
	NodeList matured;
	Node* t = ptrs[tri][num_levels - 1];
	matured.insert(t->barcode, t->mutation, t->fitness, tri, num_levels, t->count);
	cells.totalFitness -= (t->count * t->fitness);
	// Move all of the information up a level
	for(int j = num_levels-1; j >= 0; j--){
		if(j >= 1){
			Node* y = ptrs[tri][j];
			Node* x = ptrs[tri][j-1];
			y->barcode = x->barcode;
			y->mutation = x->mutation;
			y->fitness = x->fitness;
			cells.totalFitness -= (x->count * x->fitness);
			cells.totalFitness += (y->count * y->fitness);
		} else {
			Node* y = ptrs[tri][j];
			y->barcode = code;
			y->mutation = m;
			y->fitness = f;
			cells.totalFitness += f;
		}
	}
	
}

NodeList DiffTriangle::getRemoveLastLevel(int tri){
	NodeList ret;
	Node* t;
	std::vector<Node*> deleteList;
	for(t = cells.root; t != NULL; t = t->next){
		if(t->triangle == tri && t->level == num_levels){
			//t->printNode();
			ret.insert(t->barcode, t->mutation, t->fitness, tri, num_levels, t->count);
			deleteList.push_back(t);
		}
	}
	//std::cout << "made it here" << std::endl;
	for(size_t i = 0; i < deleteList.size(); i++){
		deleteList[i]->removeFromList();
	}
	return ret;
}

void DiffTriangle::doubleCells(int tri){
	Node* t;
	for(t = cells.root; t != NULL; t = t->next){
		if(t->triangle == tri){
			t->level = t->level + 1;
			cells.insert(t->barcode, t->mutation, t->fitness, tri, t->level, t->count);
		}
	}
}

NodeList DiffTriangle::removeLastAndDouble(int tri){
	NodeList ret;
	Node* t;
	std::vector<Node*> deleteList;
	for(t = cells.root; t != NULL; t = t->next){
		if(t->triangle == tri){
			if(t->level == num_levels){
				//t->printNode();
				ret.insert(t->barcode, t->mutation, t->fitness, tri, num_levels, t->count);
				deleteList.push_back(t);
			} else {
				t->level = t->level + 1;
				cells.insert(t->barcode, t->mutation, t->fitness, tri, t->level, t->count);
			}
		}
	}
	//std::cout << "made it here" << std::endl;
	for(size_t i = 0; i < deleteList.size(); i++){
		deleteList[i]->removeFromList();
	}
	return ret;
}

void DiffTriangle::printCodes()
{
    cells.traversePrint();
}

double DiffTriangle::getSumRates()
{
	return 0;
}

std::vector<double> DiffTriangle::getRates()
{
    std::vector<double> ret;
	ret.push_back(0.0);
    return ret;
}



// ORIGINAL
/*
std::vector<long int> randomIndices(std::vector<long int> ind, int size){
	//std::random_device rd;
    //std::mt19937 g(rd());
	//std::cout << "In randomIndices(): " << std::endl;
	//std::copy(ind.begin(), ind.end(), std::ostream_iterator<long int>(std::cout, " "));
	//std::cout << std::endl;
	//std::vector<long int> r;
	//for(int i = 0; i < size; i++)
	//	r.push_back(ind[i]);
	std::vector<long int> ret(ind);
	//std::cout << "HERE" << std::endl;;
	//std::cout << "HERE TOO" << std::endl;
	//for (std::vector<long int>::const_iterator i = ret.begin(); i != ret.end(); ++i)
	//	std::cout << *i << ' ';
	//std::cout << std::endl;

    std::random_shuffle(ret.begin(), ret.end());
	std::vector<long int> l(ret.begin(), ret.begin() + size);
	//std::copy(r.begin(), r.end(), std::ostream_iterator<long int>(std::cout, " "));
	//std::cout << std::endl;
	return l;
}
*/



void DiffTriangle::time_step(bool verbose)
{
    return;
}

double DiffTriangle::sdi()
{
    return(GrowingPop::sdi());
}

double DiffTriangle::coded()
{
    return(GrowingPop::coded());
}

void DiffTriangle::writeToFile(std::ofstream& of, int time){
	//
	cells.writeToFile2(of, time);
}

