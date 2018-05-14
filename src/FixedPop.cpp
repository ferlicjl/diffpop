/*
 * =====================================================================================
 *
 *       Filename:  FixedPop.cpp
 *
 *    Description:  Class for a fixed-size population
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

#include "FixedPop.h"
#include "DiffTree.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <numeric>
#include <algorithm>
#include <assert.h>     /* assert */
#include <map>
#include <math.h>
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;

// All Constructors simply call CellType constructors
FixedPop::FixedPop() : GrowingPop()
{
}

FixedPop::FixedPop(std::string cname) : GrowingPop(cname)
{

}

FixedPop::FixedPop(std::string cname, long int ncells) : GrowingPop(cname, ncells)
{
}

FixedPop::FixedPop(std::string cname, long int ncells, long int offset) : GrowingPop(cname, ncells, offset)
{

}

FixedPop::FixedPop(std::string cname, long int ncells, long int offset, double eff) : GrowingPop(cname, ncells, offset, eff)
{

}

// Destructor - no dynamic memory
FixedPop::~FixedPop()
{
}

void FixedPop::up_Event(bool verbose){
	std::vector<double> up_rates;
    up_rates.push_back(alpha);

	std::vector<std::string> up_names;
    up_names.push_back("Mitosis");

	int choice = choose(up_rates);
	if(verbose)
		std::cout << up_names[choice] << ' ';
	if(choice == 0){
		mitosis();
	} else {
		mitosis();
	}
}

void FixedPop::down_Event(bool verbose){
	verbose = false;
	std::vector<double> down_rates;
    down_rates.push_back(death);
    down_rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));
	down_rates.push_back(std::accumulate(gamma2.begin(), gamma2.end(), 0.0));
	down_rates.push_back(std::accumulate(zeta.begin(), zeta.end(), 0.0));

	std::vector<std::string> down_names;
    down_names.push_back("Apoptosis");
    down_names.push_back("Diff1");
	down_names.push_back("Diff2");
	down_names.push_back("Dedifferentiation");

	int choice = choose(down_rates);
	if(verbose)
		std::cout << down_names[choice] << ' ';
	if(choice == 0){
		apoptosis();
	} else if (choice == 1){
		diff1();
	} else if (choice == 2) {
		diff2();
	} else if (choice == 3) {
		dediff();
	} else {
		apoptosis();
	}
}

// Enacts the events that would occur in one-unit of time according to rates determined by the cell type
// Importantly, events are ordered such that cell population size remains constant
void FixedPop::time_step(bool verbose)
{
    //std::cout << "in time_step for " << name << std::endl;
	//verbose = true;
    // Select number of events from these rates
    //int n_events = getNumDownEvents() - tree->dediffTo(name);
	int n_events = getNumDownEvents();
    if(verbose)
        std::cout << "Down events: " <<  n_events << std::endl;

    // Set up a vector to be able to choose an event type
    std::vector<double> init_rates = getRates();

	// TRYING SOMETHING
	//std::cout << "HELLO" << std::endl;

    // Set up simple string to print out event type name
    std::vector<std::string> init_names;
    init_names.push_back("Mitosis");
    init_names.push_back("Asym. Differentiation");
    init_names.push_back("Diff1");
	init_names.push_back("Diff2");
	init_names.push_back("Death");
	init_names.push_back("Dedifferentiation");
	//init_names.push_back("Mutation");

	init_rates[0] = 0.0;
	//init_rates[1] = 0.0;

    // Iterate over the number of events
	if(verbose)
		std::cout << "I already made " << upstreamEvents << " events from upstream populations" << std::endl;
    for(int i = upstreamEvents; i < n_events; i++)
    {
        if(verbose)
            std::cout << "Step " << (i+1) << ": ";
        int choice = choose(init_rates);
        if(verbose)
            std::cout << init_names[choice] << ' ';

        // Gained a cell, need to remove by initiating a down event
        if(choice == 0)
        {
            mitosis();
            down_Event();
			i++;
        }
        //Asym. Diff doesn't change cell population size
        else if(choice == 1)
        {
            asym_diff();
        }
        // Differentiation causes us to need a mitosis event
        else if (choice == 2)
        {
            diff1();
			up_Event();
        }
		// Differentiation causes us to need a mitosis event
        else if (choice == 3)
        {
            diff2();
			up_Event();
        }
		// Death
		else if (choice == 4) {
			apoptosis();
			up_Event();
		}
		// Dedifferentiation
		else if (choice == 5) {
			dediff();
			up_Event();
		}
    }
    if(verbose)
        std::cout << std::endl;
    //std::cout << "end of time_step() for " << name << std::endl;
}

// Enacts the events that would occur in one-unit of time according to rates determined by the cell type
// Importantly, events are ordered such that cell population size remains constant
void FixedPop::single_event(bool verbose)
{
    int n_events = 1;
    if(verbose)
        std::cout << n_events << std::endl;

    // Set up a vector to be able to choose an event type
    std::vector<double> init_rates;
    init_rates.push_back(alpha);
    init_rates.push_back(std::accumulate(beta.begin(), beta.end(), 0.0));

    std::vector<double> down_rates;
    down_rates.push_back(death);
    down_rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));

    // Set up simple string to print out event type name
    std::vector<std::string> init_names;
    init_names.push_back("Mitosis");
    init_names.push_back("Asym. Differentiation");

    std::vector<std::string> down_names;
    down_names.push_back("Apoptosis");
    for(int j = 0; j < numChildren; j++)
    {
        down_names.push_back(children[j]->name);
    }

    // Iterate over the number of events
    for(int i = 0; i < n_events; i++)
    {
        if(verbose)
            std::cout << "Step " << (i+1) << ": ";
        int choice = choose(init_rates);
        if(verbose)
            std::cout << init_names[choice] << ' ';

        // Gained a cell, need to remove by initiating a down event
        if(choice == 0)
        {
            mitosis();
            int choice2 = choose(down_rates);
            // Cell death, or a differentiation event
            if(choice2 == 0)
                apoptosis();
            else
            {
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

// Dedifferentation event
void FixedPop::dediff()
{
    // Choose random cell to differentiate and remove
    int randomIndex = gsl_rng_uniform_int(rng, cells.total);

	Node* state = cells.getAt(randomIndex);

	// Copy down cell's information for passing on
    long int code = state->barcode;
    std::string mut = state->mutation;
    double f = state->fitness;

    cells.removeAt(randomIndex);

    // Set up differentiation pathway vectors
    std::vector<double> c_norm = normalize(zeta);
    assert(c_norm.size() == dediff_children.size());

    double r = gsl_rng_uniform(rng);
    // Choose child population to differentiate to according to random number and cumulative diff prob vector
    dediffEvents++;
    //upstreamEvents++;
    for(size_t i = 0; i < dediff_children.size(); i++)
    {
        if(r <= c_norm[i])
        {
            dediff_children[i]->add_dediffCell(code, mut, f);
            return;
        }

    }

}

// Amend adding a cell to this population since any increase in cell population must be met with a decrease
void FixedPop::add_cell(long int code, std::string mut, double f)
{
	//std::cout << "I'm calling FixedPop add_cell new" << std::endl;
    addcellEvents += 1;
    down_Event();
	upstreamEvents++;

    cells.insert(code, mut, f, 1);

    if(mut != "0")
        numMutCells++;
	//std::cout << "finished add_cell FPN" << std::endl;
}


double FixedPop::getNumEvents()
{
    return (alpha + std::accumulate(beta.begin(), beta.end(), 0.0) + std::accumulate(gamma1.begin(), gamma1.end(), 0.0) + death) * cells.total;
}
