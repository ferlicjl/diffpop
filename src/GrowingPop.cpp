/*
 * =====================================================================================
 *
 *       Filename:  GrowingPop.cpp
 *
 *    Description:  Class for an exponentially growing population
 *
 *        Version:  1.0
 *        Created:  03/01/2017 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University,
 *
 * =====================================================================================
 */

#include "GrowingPop.h"
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
#include <random>
#include <iterator>
#include <cmath>
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;


// Default Constuctor
GrowingPop::GrowingPop()
{
    name = "Default";
    parent = NULL;
	tree = NULL;
    numChildren = 0;
    numCells = 0;
    alpha = 0;
    death = 0;
    mu = 0;
    numMutCells = 0;
    maxMutation = 2;
	numMutations = maxMutation;
    upstreamEvents = 0;
    efficiency = 1.0;
    alphaEvents = 0;
    gamma1Events = 0;
	gamma2Events = 0;
    deathEvents = 0;
    betaEvents = 0;
    dediffEvents = 0;
    addcellEvents = 0;
	mutateEvents = 0;
}

// Constructor for a GrowingPop with given name
GrowingPop::GrowingPop(std::string cname) : GrowingPop()
{
    name = cname;
}

// Constructor for a GrowingPop given a name & number of cells
GrowingPop::GrowingPop(std::string cname, long int ncells) : GrowingPop(cname, ncells, 0)
{
}

// Constructor for a GrowingPop given a name, number of cells, and barcoding offset
GrowingPop::GrowingPop(std::string cname, long int ncells, long int offset) : GrowingPop(cname)
{
    numCells = ncells;
    for(long int i = 0 + offset; i < ncells + offset; i++)
    {
        //codes.push_back(i);
        //mutation.push_back(0);
        //fitness.push_back(1);

		cells.insert(i, "0", 1, 1);
    }
}

// Constructor for a GrowingPop given a name, number of cells, , barcoding offset, and a labelling efficiency
// Note: Unlabelled cells will have code = -1
GrowingPop::GrowingPop(std::string cname, long int ncells, long int offset, double eff) : GrowingPop(cname)
{
    numCells = ncells;
    efficiency = eff;
	// Generate cells if a bernoulli(eff) is a success
    for(long int i = 0 + offset; i < ncells + offset; i++)
    {
        int labelled = gsl_ran_bernoulli(rng, eff);
        if(labelled > 0)
        {
			cells.insert(i, "0", 1, 1);
        }
        else
        {
			cells.insert(-1, "0", 1, 1);
        }
    }
}

// Destructor - no dynamic memory items as of yet
GrowingPop::~GrowingPop()
{
	cells.deleteList();
}

// Add a child GrowingPop with given differentiation rate
// sets the child's parent to this GrowingPop
void GrowingPop::addChild(GrowingPop* child, double diff1_rate, double diff2_rate, double asym_rate)
{
	std::vector<double> rs = {diff1_rate, diff2_rate, asym_rate};

	if ( loadDiffRates.find(child->name) == loadDiffRates.end() ) {
		// I do NOT already have a Child population child
		loadDiffRates.emplace(child->name, rs);
		c_map.emplace(child->name, child);
	} else {
		// I already have a Child population, just add a new rate
		loadDiffRates[child->name] = rs;
}

}

// Add a dedifferentation population with given rate
// This population should be somewhere upstream of the current population
void GrowingPop::addDediff(GrowingPop* to, double rate)
{
    dediff_children.push_back(to);
    zeta.push_back(rate);
}

// Set a parent for a given GrowingPop
// Should not typically be used because the parent does not inherit the current GrowingPop as a child - no differentiation to this GrowingPop
void GrowingPop::setParent(GrowingPop* p)
{
    parent = p;
}

// Associates the GrowingPop with its DiffTree
void GrowingPop::setTree(DiffTree* t)
{
    tree = t;
}


// Set the number of events in this GrowingPop, caused by upstream populations
void GrowingPop::setUpstreamEvents(int x)
{
    upstreamEvents = x;
}

// Zero out all event counters
void GrowingPop::zeroEvents()
{
    alphaEvents = 0;
    gamma1Events = 0;
	gamma2Events = 0;
    deathEvents = 0;
    betaEvents = 0;
    dediffEvents = 0;
	mutateEvents = 0;
    addcellEvents = 0;
}

void GrowingPop::setDiff1(GrowingPop* child, double r){
	if ( loadDiffRates.find(child->name) == loadDiffRates.end() ) {
		// not found
		std::vector<double> rs = {r, 0, 0};
		loadDiffRates.emplace(child->name, rs);
		c_map.emplace(child->name, child);
	} else {
		// found in map already
		loadDiffRates[child->name][0] = r;
	}
}

void GrowingPop::setDiff2(GrowingPop* child, double r){
	if ( loadDiffRates.find(child->name) == loadDiffRates.end() ) {
		// not found
		std::vector<double> rs = {0, r, 0};
		loadDiffRates.emplace(child->name, rs);
		c_map.emplace(child->name, child);
	} else {
		// found in map already
		loadDiffRates[child->name][1] = r;
	}
}

void GrowingPop::setAsymdiff(GrowingPop* child, double r){
	if ( loadDiffRates.find(child->name) == loadDiffRates.end() ) {
		// not found
		std::vector<double> rs = {0, 0, r};
		loadDiffRates.emplace(child->name, rs);
		c_map.emplace(child->name, child);
	} else {
		// found in map already
		loadDiffRates[child->name][2] = r;
	}
}

// Calculate delta based on other parameter values in order to maintain fixed state
// NOTE: we use d as the rate in our simulation because all de-differentiation to this population will result in death
double GrowingPop::calcDelta(){
	double g1 = std::accumulate(gamma1.begin(), gamma1.end(), 0.0);
	double g2 = std::accumulate(gamma2.begin(), gamma2.end(), 0.0);
	double z = std::accumulate(zeta.begin(), zeta.end(), 0.0);
	double d = alpha*numCells + tree->gamma1To(name) + tree->gamma2To(name) + tree->betaTo(name) - numCells*(g1 + g2 + z);
	double d2 = alpha*numCells + tree->gamma1To(name) + tree->gamma2To(name) + tree->betaTo(name) - numCells*(g1 + g2 + z) + tree->dediffTo(name);

	std::cout << "\tEffective Death Rate: " << (d2 / numCells) << std::endl;

	double ret = d / numCells;

	//std::cout << "\tReturned Death Rate: " << ret << std::endl;

	return ret;
}

// NOT WORKING
// Calculate delta based on other parameter values in order to maintain fixed state
// NOTE: we use d as the rate in our simulation because all de-differentiation to this population will result in death
double GrowingPop::calcAlpha(){
	double g1 = std::accumulate(gamma1.begin(), gamma1.end(), 0.0);
	double g2 = std::accumulate(gamma2.begin(), gamma2.end(), 0.0);
	double z = std::accumulate(zeta.begin(), zeta.end(), 0.0);
	
	death = death*cells.total - tree->dediffTo(name);

	// IF death is not greater than dedifferentiation, throw error
	if(death*cells.total < tree->dediffTo(name))
		std::cout << "Error: population " << name << " delta rate not sufficient to cover dedifferntiation to " << name << std::endl;

	double a = (death*cells.total) + cells.total * (g1 + g2 + z) - tree->gamma1To(name) - tree->gamma2To(name) - tree->betaTo(name);

	//double d = alpha*numCells + tree->gamma1To(name) + tree->gamma2To(name) + tree->betaTo(name) - numCells*(g1 + g2 + z);
	//double d2 = alpha*numCells + tree->gamma1To(name) + tree->gamma2To(name) + tree->betaTo(name) - numCells*(g1 + g2 + z) + tree->dediffTo(name);

	std::cout << "\tCalculated Alpha Rate: " << (a / cells.total) << std::endl;

	double ret = a / cells.total;

	//std::cout << "\tReturned Death Rate: " << ret << std::endl;

	return ret;
}

// Override << operator - print out GrowingPop name
std::ostream &operator<<(std::ostream &os, const GrowingPop& c)
{
    os << c.name;
    return os;
}

// a nice method to print out a tree-like structure
void GrowingPop::printPretty(int indent)
{
    std::cout << std::setw(indent) << name << std::endl;
    for(int i = 0; i < numChildren; i++)
    {
        children[i]->printPretty(indent + 4);
    }
}

// Print out all barcodes for this GrowingPop
void GrowingPop::printCodes()
{
    cells.traversePrint();
}

// Used to return all event rates as a vector
std::vector<double> GrowingPop::getRates()
{
    std::vector<double> ret;
    ret.push_back(alpha);
	ret.push_back(std::accumulate(beta.begin(), beta.end(), 0.0));
    ret.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));
	ret.push_back(std::accumulate(gamma2.begin(), gamma2.end(), 0.0));
    ret.push_back(death);
    ret.push_back(std::accumulate(zeta.begin(), zeta.end(), 0.0));
    return ret;
}

// Used to return all "i-1" and "i" rates as a vector
std::vector<double> GrowingPop::getDownRates()
{
    std::vector<double> ret;
	ret.push_back(std::accumulate(beta.begin(), beta.end(), 0.0));
    ret.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));
	ret.push_back(std::accumulate(gamma2.begin(), gamma2.end(), 0.0));
    ret.push_back(death);
    ret.push_back(std::accumulate(zeta.begin(), zeta.end(), 0.0));
    return ret;
}

// Used to determine the total number of expected events in one time unit
double GrowingPop::getNumEvents()
{
	std::vector<double> rates = getRates();
	double ret = std::accumulate(rates.begin(), rates.end(), 0.0) * cells.total;

    return ret;
}

// Used to determine the total number of "i" and "i-1" events in one time unit
double GrowingPop::getNumDownEvents()
{
	std::vector<double> rates = getDownRates();
	double ret = std::accumulate(rates.begin(), rates.end(), 0.0) * cells.total;

    return ret;
}

// Provides access to CellPopulation children
std::vector<GrowingPop*> GrowingPop::getChildren()
{
    return children;
}

// Mitotic event - should occur at an 'alpha' rate per unit time
void GrowingPop::mitosis()
{
	// If we don't have any mutated cells, all are of same fitness
    if(cells.numMutated > 0)
        cells.setCurFitness(gsl_rng_uniform(rng) * cells.totalFitness);
    else
        cells.setCurIndex(gsl_rng_uniform_int(rng, cells.total));


	// Add cell to our population, we increase the cell state set in the previous step
    cells.increaseAtCur();

	// Mutation
	int mutation = gsl_ran_bernoulli(rng, mu);
	if (mutation == 1){
		mutateCur();
	}

    // Increase mitosis event counter
    alphaEvents++;
}

// Death event - should occur at the 'death' rate per unit time
void GrowingPop::apoptosis()
{
    // Choose a random cell
    int randomIndex = gsl_rng_uniform_int(rng, cells.total);

	// Remove it from population vectors
    cells.removeAt(randomIndex);

	// Increase apoptosis event counter
    deathEvents++;
}

// Differentiation event - 1-to-1
void GrowingPop::diff1()
{
    // Choose random cell to differentiate and remove
    int randomIndex = gsl_rng_uniform_int(rng, cells.total);

	Node* state = cells.getAt(randomIndex);

	// Copy down cell's information for passing on
    long int code = state->barcode;
    std::string mut = state->mutation;
    double f = state->fitness;

	// Remove cell from current population
    cells.removeAt(randomIndex);

    // Set up differentiation pathway vectors
    std::vector<double> c_norm = normalize(gamma1);
    assert(c_norm.size() == children.size());

	// Choose child population to differentiate to according to random number and cumulative diff prob vector
    double r = gsl_rng_uniform(rng);
    gamma1Events++;


	// If we haven't gotten rid of the cell yet, find which population to add to
	bool stillAdd = true;
    for(size_t i = 0; i < children.size(); i++)
    {
        if(r <= c_norm[i] && stillAdd)
        {
            children[i]->add_cell(code, mut, f);
			stillAdd = false;
            return;
        }
    }
}

// Differentiation event - 1-to-2
void GrowingPop::diff2()
{
    // Choose random cell to differentiate and remove
    int randomIndex = gsl_rng_uniform_int(rng, cells.total);

	Node* state = cells.getAt(randomIndex);

	// Copy down cell's information for passing on
    long int code = state->barcode;
    std::string mut = state->mutation;
    double f = state->fitness;

	// Remove cell from current population
    cells.removeAt(randomIndex);

    // Set up differentiation pathway vectors
    std::vector<double> c_norm = normalize(gamma2);
    assert(c_norm.size() == children.size());

	// Choose child population to differentiate to according to random number and cumulative diff prob vector
    double r = gsl_rng_uniform(rng);

	// Increase gamma2 event counter
    gamma2Events++;

	// Find which population we need to add to
    for(size_t i = 0; i < children.size(); i++)
    {
        if(r <= c_norm[i])
        {
            children[i]->add_cell(code, mut, f);
			children[i]->add_cell(code, mut, f);
            return;
        }
    }
}

// Asymmetric differentiation event
void GrowingPop::asym_diff()
{
    // Choose random cell to differentiate
    int randomIndex = gsl_rng_uniform_int(rng, cells.total);

	Node* state = cells.getAt(randomIndex);

	// Copy down cell's information for passing on
    long int code = state->barcode;
    std::string mut = state->mutation;
    double f = state->fitness;

    // Set up differentiation pathway vectors
    std::vector<double> c_norm = normalize(beta);
    assert(c_norm.size() == children.size());

	// Increase beta event counter
	betaEvents++;

	// Choose child population to differentiate to according to random number and cumulative diff prob vector
    double r = gsl_rng_uniform(rng);
    for(size_t i = 0; i < children.size(); i++)
    {
        if(r <= c_norm[i])
        {
            children[i]->add_cell(code, mut, f);
            return;
        }
    }
}

// Dedifferentation event
void GrowingPop::dediff()
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

	// Increase zeta event counter
	dediffEvents++;

	// Choose child population to differentiate to according to random number and cumulative diff prob vector
    double r = gsl_rng_uniform(rng);

    for(size_t i = 0; i < dediff_children.size(); i++)
    {
        if(r <= c_norm[i])
        {
            dediff_children[i]->add_cell(code, mut, f);
            return;
        }
    }
}

// Mutation
// Old method, slow due to looking up cell index again
void GrowingPop::mutate(long int index)
{
	// Get next mutation from DiffTree
	int mutation = tree->getNextMutant();

	// Save off state information for mutated cell
	Node* state = cells.getAt(index);
	std::string mut = state->mutation;

	// Get our fitness change from the DiffTree
	double fitIncrease;
	tree->ConstantGenerateFitness(&fitIncrease, &tree->fp, rng);

	double newFitness = state->fitness + fitIncrease;

	if(newFitness < tree->fp.lower_fitness)
		newFitness = tree->fp.lower_fitness;

	if(newFitness > tree->fp.upper_fitness)
		newFitness = tree->fp.upper_fitness ;


	// Insert a new cell with updated state (mutation and fitness change)
	cells.insert(state->barcode, state->mutation + ">" + std::to_string(mutation), newFitness, 1);

	// Remove a cell from the old state
	cells.removeAt(index);

	// Write out mutation information to file
	tree->writeMutant(name, mutation, newFitness);

	// Increase mutation event counter
	mutateEvents++;
}

// Mutation at NodeList cur
void GrowingPop::mutateCur()
{
	// Get next mutation from DiffTree
	int mutation = tree->getNextMutant();

	// Save off cell state and remove one from it
	Node* state = cells.getCur();
	std::string mut = state->mutation;
	cells.removeAtCur();


	// Get our fitness change from the DiffTree
	double fitIncrease;
	tree->ConstantGenerateFitness(&fitIncrease, &tree->fp, rng);

	double newFitness = state->fitness + fitIncrease;

	if(newFitness < tree->fp.lower_fitness){
		newFitness = tree->fp.lower_fitness;
	} else if(newFitness > tree->fp.upper_fitness){
		newFitness = tree->fp.upper_fitness ;
	}


	// Insert a new cell with updated state (mutation and fitness change)
	cells.insert(state->barcode, state->mutation + ">" + std::to_string(mutation), newFitness, 1);

	// Output mutation information
	tree->writeMutant(name, mutation, newFitness);

	// Increase mutation event counter
	mutateEvents++;
}

// Erase a cell from the population given its index
void GrowingPop::erase_cell(long int index){
	if(cells.getAt(index)->mutation != "0")
        numMutCells--;
	cells.removeAt(index);
}


// Add a code to this GrowingPop population
void GrowingPop::add_cell(long int code, std::string mut, double f)
{
    cells.insert(code, mut, f, 1);
    addcellEvents += 1;
    if(mut != "0")
        numMutCells++;
}

// Amend adding a cell to this population since any increase in cell population must be met with a decrease
// NOTE: this method is ONLY used in FixedPops, should probably move there
// NOTE: automatically results in an apoptosis event to keep population size constant
void GrowingPop::add_dediffCell(long int code, std::string mut, double f)
{
    addcellEvents += 1;
	apoptosis();

    // Add new cell to population finally
   cells.insert(code, mut, f, 1);

    if(mut != "0")
        numMutCells++;
}

// Initialize the codes before a first run
// -> ensures that unique barcodes are given to all cells
void GrowingPop::initializeCells(long int offset)
{
    cells.deleteList();

    // Generate barcode if a bernoulli(eff) is a success
    for(long int i = 0 + offset; i < numCells + offset; i++)
    {
        int labelled = gsl_ran_bernoulli(rng, efficiency);
        if(labelled > 0)
        {
            cells.insert(i, "0", 1, 1);
        }
        else
        {
            cells.insert(-1, "0", 1, 1);
        }
    }

    return;
}

// Initialize differentiation rates (beta, gamma1, gamma2)
void GrowingPop::initializeTree(){
	for (std::map<std::string,GrowingPop*>::iterator it=c_map.begin(); it!=c_map.end(); ++it){
		children.push_back(it->second);
		numChildren++;
		gamma1.push_back(loadDiffRates[it->first][0]);
		gamma2.push_back(loadDiffRates[it->first][1]);
		beta.push_back(loadDiffRates[it->first][2]);
	}
	return;
}


// Enacts the events that would occur in one-unit of time according to rates determined by the cell type
void GrowingPop::time_step(bool verbose)
{
	int n_events = getNumEvents();
	if(verbose)
		std::cout << n_events << std::endl;

    // Set up a vector to be able to choose an event type
    std::vector<double> rates;
    rates.push_back(alpha);
    rates.push_back(std::accumulate(beta.begin(), beta.end(), 0.0));
    rates.push_back(death);
    rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));
	rates.push_back(std::accumulate(gamma2.begin(), gamma2.end(), 0.0));
    rates.push_back(std::accumulate(zeta.begin(), zeta.end(), 0.0));

    // Set up simple string to print out event type name
    std::vector<std::string> event_names;
    event_names.push_back("Mitosis");
    event_names.push_back("Asym. Differentiation");
    event_names.push_back("Apoptosis");
    event_names.push_back("1-to-1 Differentiation");
	event_names.push_back("1-to-2 Differentiation");
    event_names.push_back("Dedifferentiation");

    // Iterate over the number of events
    for(int i = 0; i < n_events; i++)
    {
        if(verbose)
            std::cout << "Step " << (i+1) << ": ";

		// Choose which event occurs
        int choice = choose(rates);

        if(verbose)
            std::cout << event_names[choice] << " ";

        if(choice == 0)
            mitosis();
        else if (choice == 1)
            asym_diff();
        else if (choice == 2)
            apoptosis();
        else if (choice == 3)
            diff1();
		else if (choice == 4)
            diff2();
        else if (choice == 5)
            dediff();
    }
	if(verbose)
		std::cout << std::endl;
}

void GrowingPop::single_event(bool verbose)
{
    //verbose = true;
	int n_events = 1;
    if(verbose)
        std::cout << n_events << std::endl;

    // Set up a vector to be able to choose an event type
    std::vector<double> rates;
    rates.push_back(alpha);
    rates.push_back(std::accumulate(beta.begin(), beta.end(), 0.0));
    rates.push_back(death);
    rates.push_back(std::accumulate(gamma1.begin(), gamma1.end(), 0.0));
	rates.push_back(std::accumulate(gamma2.begin(), gamma2.end(), 0.0));
    rates.push_back(std::accumulate(zeta.begin(), zeta.end(), 0.0));

    // Set up simple string to print out event type name
    std::vector<std::string> event_names;
    event_names.push_back("Mitosis");
    event_names.push_back("ADiff");
    event_names.push_back("Apoptosis");
    event_names.push_back("Diff1-");
	event_names.push_back("Diff2-");
    event_names.push_back("Dediiff-");

    // Iterate over the number of events
    for(int i = 0; i < n_events; i++)
    {
        if(verbose)
            std::cout << "Step " << (i+1) << ": ";
		// Choose event
        int choice = choose(rates);
        if(verbose)
            std::cout << event_names[choice] << " ";
        if(choice == 0)
            mitosis();
        else if (choice == 1)
            asym_diff();
        else if (choice == 2)
            apoptosis();
        else if (choice == 3)
            diff1();
		else if (choice == 4)
            diff2();
        else if (choice == 5)
            dediff();
    }
    if(verbose)
        std::cout << std::endl;
}

// Print event rates for current GrowingPop
void GrowingPop::print_rates()
{
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "death: " << death << std::endl;
    for(int i = 0; i < numChildren; i++)
    {
		std::cout << children[i]->name << "Diff1 " << ": " << gamma1[i] << std::endl;
		std::cout << children[i]->name << "Diff2 " << ": " << gamma2[i] << std::endl;
        std::cout << children[i]->name << " Asymmetric Diff" << ": " << beta[i] << std::endl;

    }
    for(size_t i = 0; i < dediff_children.size(); i++)
    {
        std::cout << dediff_children[i]->name << " Dediff" << ": " << zeta[i] << std::endl;
    }
    std::cout << " Mutation: " << mu << std::endl;
}

// Set Alpha - mitotic event rate
void GrowingPop::setAlpha(double a)
{
    alpha = a;
}

// Set death - apoptosis event rate
void GrowingPop::setDeath(double d)
{
    death = d;
}

// Set mu - mutation event rate
void GrowingPop::setMu(double m)
{
    mu = m;
}


// Shannon Diversity Index
double GrowingPop::sdi()
{
   return cells.sdi();
}

std::string GrowingPop::diversity(){
	return cells.diversity();
}

// Returns percent of cells with barcode
double GrowingPop::coded()
{
    return cells.labelled();

}

// Sample from codes
// NOT USED!!!
std::vector<long int> GrowingPop::sample(int size)
{
    std::vector<long int> ret;
	/*
    ret = codes;

	gsl_ran_shuffle(rng, &ret[0], ret.size(), sizeof(ret[0]));

    std::vector<long int> l(ret.begin(), ret.begin() + size);
    return l;
	*/
	return ret;
}

void GrowingPop::writeToFile(std::ofstream& of, int time){
	cells.writeToFile(of, time);
}


