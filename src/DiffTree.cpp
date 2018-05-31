/*
 * =====================================================================================
 *
 *       Filename:  DiffTree.cpp
 *
 *    Description:  Class for calling functions throughout a hierarchy
 *
 *        Version:  1.0
 *        Created:  03/01/2017 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */

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
#include <deque>
#include <limits>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <ctime>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>


extern gsl_rng* rng;

// CONSTRUCTORS
// Default Constructor:
// - set things to defaults
// - generate output file codes (dates + random string)

DiffTree::DiffTree()
{
    root = NULL;
	mutationCounter = 0;
	simTime = 0;

    // Get date-time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);

    // Generate a random number in case jobs kick off in same second
    unsigned long int a = gsl_rng_uniform_int (rng, 100000);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%H%M%S",timeinfo);
    std::string str(buffer);

    // Combine date-time / random number for output name
    std::ostringstream o;
    o << "out_" << str << "_" << a;
    opath = o.str();
}

// Constructor given output directory
DiffTree::DiffTree(std::string outdir)
{
    root = NULL;
	mutationCounter = 0;
	simTime = 0;

    // Get date-time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);

    // Generate a random number in case jobs kick off in same second
    unsigned long int a = gsl_rng_uniform_int (rng, 100000);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%H%M%S",timeinfo);
    std::string str(buffer);

    std::ostringstream o;
    o << outdir << "out_" << str << "_" << a;
    opath = o.str();

	// Set default fitness distribution
	fp.fitness_distribution = "normal";
	fp.is_randfitness = true;
	fp.alpha_fitness = 0;
	fp.beta_fitness = 1;
	fp.pass_prob = 1;
	fp.upper_fitness = std::numeric_limits<double>::max();
	fp.lower_fitness = std::numeric_limits<double>::lowest();

	setFitnessDistribution();
}

// Constructor given a root CellPopulation
DiffTree::DiffTree(GrowingPop* r)
{
    // Set root
    root = r;

    // Get date-time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);

    // Generate a random number in case jobs kick off in same second
    unsigned long int a = gsl_rng_uniform_int (rng, 100000);
    strftime(buffer,sizeof(buffer),"%d-%m-%Y-%H%M%S",timeinfo);
    std::string str(buffer);

    // Combine date-time / random number for output name
    std::ostringstream o;
    o << "out_" << str << "_" << a;
    opath = o.str();
}

// Constructor given root and output path/prefix
DiffTree::DiffTree(GrowingPop* r, std::string o)
{
    root = r;
    opath = o;
}

// Constructor given root, output path/prefix, and a map of CellPopulations
DiffTree::DiffTree(GrowingPop* r, std::string o, std::map<std::string, GrowingPop*> m_)
{
    root = r;
    opath = o;
    m = m_;
}

// Constructor given an output path/prefix and map of CellPopulations
DiffTree::DiffTree(std::string o, std::map<std::string, GrowingPop*> m_)
{
    opath = o;
    m = m_;
}

// Destructors
// No dynamic memory in tree
DiffTree::~DiffTree()
{

}

// Set Fitness Distribution
void DiffTree::setFitnessDist(FitnessParameters f){
	fp.fitness_distribution = f.fitness_distribution;
	fp.is_randfitness = f.is_randfitness;
	fp.alpha_fitness = f.alpha_fitness;
	fp.beta_fitness = f.beta_fitness;
	fp.pass_prob = f.pass_prob;
	fp.upper_fitness = f.upper_fitness;
	fp.lower_fitness = f.lower_fitness;

	setFitnessDistribution();
}

// Set Root
void DiffTree::setRoot(GrowingPop* r)
{
    root = r;
}

// Zero out all events throughout the tree
// Used between observation times
void DiffTree::zeroEvents()
{
    for(auto const& x: m)
    {
        x.second->zeroEvents();
    }
}

// Set map of CellPopulations
void DiffTree::setMap(std::map<std::string, GrowingPop*> m_)
{
    m = m_;
}


// Get Methods

// Get a vector holding all of the rates for each CellPopulation in the map
std::vector<double> DiffTree::getRates()
{
    std::vector<double> rates;

    for(auto const& x: m)
    {
        rates.push_back(x.second->getNumEvents());
    }

    return rates;
}

// Get zero
bool DiffTree::checkZero()
{
    bool zero = false;

    for(auto const& x: m)
    {
        zero = zero || x.second->cells.total == 0;
    }

    return zero;
}

std::vector<GrowingPop*> DiffTree::getCellTypes()
{
    std::vector<GrowingPop*> types;
    for(auto const& x: m)
    {
        types.push_back(x.second);
    }
    return types;
}

int DiffTree::getNextMutant(){
	mutationCounter++;
	return mutationCounter;
}

double DiffTree::nextEventTime()
{
    double rand_next_time;

	// Go through tree and get a vector of event rates throughout tree
    std::vector<double> rates = getRates();

	// Sum event rates
    double tot_rate = std::accumulate(rates.begin(), rates.end(), 0.0);

	// Generate exponential based on tot_rate
    rand_next_time = gsl_ran_exponential(rng, 1 / tot_rate);

    return rand_next_time;
}


// update State
void DiffTree::updateState(double currentTime, double timeToNext)
{
	simTime = currentTime + timeToNext;
    std::vector<double> rates = getRates();
    std::vector<GrowingPop*> types = getCellTypes();
	assert(rates.size() == types.size());
    types[choose(rates)]->single_event();

    return;
}

// Run simulations for n time units
void DiffTree::simulate(int numTime, int outputFrequency)
{
	//initializeTree();
	//initializeCells();

	//std::cout << "I've initialized the tree!" << std::endl;

	bool verbose = true;

    // NEW outputs
    std::ofstream sdi_of;
	//std::ofstream types_of;
    std::ofstream pop_of;
    std::ofstream label_of;
    std::ofstream events_of;
	std::ofstream mutants_of;

    sdi_of.open(opath + "_diversity.csv");
	//types_of.open(opath + "_types.csv");
    pop_of.open(opath + "_pop.csv");
    label_of.open(opath + "_label.csv");
    events_of.open(opath + "_events.csv");
	mutants_of.open(opath + "_mut.csv");

	  //writeAll(0);
    writeDiversityHeader(sdi_of);
	//writeTypesHeader(types_of);
    writeCellHeader(pop_of);
    writeCellHeader(label_of);
    writeEventsHeader(events_of);
	writeMutantsHeader(mutants_of);
	writeAllHeaders();

    // Set up vector of observation times
    std::vector<int> obsTimes;
    for (int k = 0; k < numTime + 1; k++)
    {
        obsTimes.push_back(k);
    }
    //std::cout << "numTime: " << numTime << std::endl;

    // Set variables to keep track of our current time and which observation time comes next
    double curTime = 0;
    int curObsIndex = 0;

    // Keep track of number of events.  Why not?
    int numEvents = 0;

	int obsMod = pow(10, round(log10(numTime)-1));

    // Display some stuff - debugging
	if(verbose){
    //std::cout << "numTime: " << numTime << std::endl;
    //std::cout << "Simulation Start Time: " << curTime << std::endl;
    //std::cout << "Simulation End Time: " << obsTimes[numTime] << std::endl;
    //std::cout << "obsTimes.size(): " << obsTimes.size() << std::endl;
	}

    // Run until our currentTime is greater than our largest Observation time
    while(curTime <= obsTimes[numTime])
    {
        Rcpp::checkUserInterrupt();
        if(checkZero())
        {
            std::cout << "A population has size 0, ending simulation" << std::endl;
            std::cout << "End time: " << curTime << std::endl;
            break;
        }

        // Get the next event time
        double timeToNext = nextEventTime();
        //std::cout << "Current Time: " << curTime << std::endl;
        //std::cout << "Next: " << curTime + timeToNext << std::endl;

        // If our next event time is later than observation times,
        // Make our observations
        while((curTime + timeToNext > obsTimes[curObsIndex]))// & (curTime + timeToNext <= obsTimes[numTime]))
        {
			//writeAll(obsTimes[curObsIndex]);
            writeSDI(sdi_of, obsTimes[curObsIndex]);
			//writeTypes(types_of, obsTimes[curObsIndex]);
            writePopSize(pop_of, obsTimes[curObsIndex]);
            writeLabelled(label_of, obsTimes[curObsIndex]);
            writeEvents(events_of, obsTimes[curObsIndex]);
            //writeMutations(obsTimes[curObsIndex]);
			if((outputFrequency > 0 && obsTimes[curObsIndex] % outputFrequency == 0) || curTime == 0)
				writeAll(obsTimes[curObsIndex]);
            zeroEvents();
			if(verbose &&  obsTimes[curObsIndex] % obsMod == 0)
				std::cout << "Time " << obsTimes[curObsIndex] << " of " << numTime << std::endl;
            curObsIndex++;
            //std::cout << "Sucessfully made an observation" << std::endl;
			if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;
        }
		if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;
        //std::cout << "trying to update state!" << std::endl;
        // Update our state
        updateState(curTime, timeToNext);
       // std::cout << "Updated state: " << curTime << " next: " << timeToNext << std::endl;
        numEvents++;

        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
        //timeToNext = nextEventTime();
    }
	//std::cout << "curObsIndex: " << curObsIndex << std::endl;
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex] << std::endl;
    std::cout << "Number of total events: " << numEvents << std::endl;

	writeAll(obsTimes[curObsIndex]);
	writeSDI(sdi_of, obsTimes[curObsIndex]);
	//writeTypes(types_of, obsTimes[curObsIndex]);
    writePopSize(pop_of, obsTimes[curObsIndex]);
    writeLabelled(label_of, obsTimes[curObsIndex]);
    writeEvents(events_of, obsTimes[curObsIndex]);
    //writeMutations(obsTimes[curObsIndex]);
    sdi_of.close();
    pop_of.close();
    label_of.close();
    events_of.close();
	//types_of.close();

    std::ofstream done_file;
    done_file.open(opath+".done");
    done_file << "Done" << std::endl;
    done_file.close();
}

void DiffTree::time_steps(int n, int outputFrequency)
{
    //initializeTree();
	//std::cout << "made it here" << std::cout;
	bool verbose = true;
	int obsMod = pow(10, round(log10(n)-1));
	//initializeCells();

	//std::cout << "can't make it here though..." << std::cout;
    std::ofstream sdi_of;
    std::ofstream pop_of;
    std::ofstream label_of;
    std::ofstream events_of;
	//std::ofstream types_of;
	std::ofstream mutants_of;

    sdi_of.open(opath + "_diversity.csv");
    pop_of.open(opath + "_pop.csv");
    label_of.open(opath + "_label.csv");
    events_of.open(opath + "_events.csv");
	//types_of.open(opath + "_types.csv");
	mutants_of.open(opath + "_mut.csv");

    writeDiversityHeader(sdi_of);
    writeCellHeader(pop_of);
    writeCellHeader(label_of);
    writeEventsHeader(events_of);
	//writeTypesHeader(types_of);
	writeMutantsHeader(mutants_of);
	writeAllHeaders();
	writeAll(0);

    writeSDI(sdi_of, 0);
    writePopSize(pop_of, 0);
    writeLabelled(label_of, 0);
	//writeTypes(types_of, 0);
    //writeEvents(events_of, 0);

	//std::cout << "Traverse Frequency: " << outputFrequency << std::endl;
	std::cout << "Time 0 of " << n << std::endl;

    for(int i = 1; i <= n; i++)
    {
		simTime = i;
        Rcpp::checkUserInterrupt();
		if(verbose && i % obsMod == 0)
			std::cout << "Time " << i << " of " << n << std::endl;
        time_step();
        writeSDI(sdi_of, i);
        writePopSize(pop_of, i);
        writeLabelled(label_of, i);
        writeEvents(events_of, i);
		//writeTypes(types_of, i);
		//std::cout << "outputFrequency > 0: " << (outputFrequency > 0)  << " n % outputFrequency == 0: " << (n % outputFrequency == 0) << std::endl;
		if((outputFrequency > 0 && i % outputFrequency == 0) || i == n)
				writeAll(i);
        //writeMutations(i);
        zeroEvents();
    }

    sdi_of.close();
    pop_of.close();
    label_of.close();
    events_of.close();
	//types_of.close();
	//writeAll(n);

    std::ofstream done_file;
    done_file.open(opath+".done");
    done_file << "Done" << std::endl;
    done_file.close();
}

// OUTPUT METHODS
void DiffTree::writeCellHeader(std::ofstream& of)
{
    of << "time";

    for(auto const& x: bfs)
    {
        of << "," << m[x]->name;
    }

    of << std::endl;
}

void DiffTree::writeDiversityHeader(std::ofstream& of)
{
    of << "time";

    for(auto const& x: bfs)
    {
        of << "," << "species_" << m[x]->name << ",shannon_" << m[x]->name << ",simpson_" << m[x]->name;
    }

    of << std::endl;
}

void DiffTree::writeTypesHeader(std::ofstream& of)
{
    of << "time";

    for(auto const& x: bfs)
    {
        of << "," << m[x]->name;
    }

    of << std::endl;
}

void DiffTree::writeEventsHeader(std::ofstream& of)
{
    of << "time";
    for(auto const& x: bfs)
    {
        of << "," << m[x]->name << "_alpha";
        of << "," << m[x]->name << "_gamma1";
        of << "," << m[x]->name << "_gamma2";
        of << "," << m[x]->name << "_death";
        of << "," << m[x]->name << "_beta";
        of << "," << m[x]->name << "_addcell";
		of << "," << m[x]->name << "_dediff";
		of << "," << m[x]->name << "_mutation";
    }
    of << std::endl;
}

void DiffTree::writeMutantsHeader(std::ofstream& of)
{
    of << "mutant" << ",";
	of << "time" << ",";
	of << "population" << ",";
	of << "fitness";
    of << std::endl;
	of.close();
}

void DiffTree::writeSDI(std::ofstream& of, int time)
{
    of << time;

    for(auto const& x: bfs)
    {
        of << "," << m[x]->diversity();
        //of << 0 << "\t";
    }

    of << std::endl;
}

void DiffTree::writeTypes(std::ofstream& of, int time)
{
    of << time;

    for(auto const& x: bfs)
    {
        of << "," << m[x]->cells.numTypes();
        //of << 0 << "\t";
    }

    of << std::endl;
}

void DiffTree::writePopSize(std::ofstream& of, int time)
{
    of << time;

    for(auto const& x: bfs)
    {
        of << "," << m[x]->cells.total;
    }

    of << std::endl;
}

void DiffTree::writeLabelled(std::ofstream& of, int time)
{
    of << time;

    for(auto const& x: bfs)
    {
        of << "," << m[x]->coded();
    }


    of << std::endl;
}

void DiffTree::writeEvents(std::ofstream& of, int time)
{
    of << time;

    for(auto const& x: bfs)
    {
        of << "," << (double) m[x]->alphaEvents;
        of << "," << (double) m[x]->gamma1Events;
        of << "," << (double) m[x]->gamma2Events;
        of << "," << (double) m[x]->deathEvents;
        of << "," << (double) m[x]->betaEvents;
        of << "," << (double) m[x]->addcellEvents;
		of << "," << (double) m[x]->dediffEvents;
		of << "," << (double) m[x]->mutateEvents;
    }
    of << std::endl;
}

/*
void DiffTree::writeMutations(int time)
{

    std::ofstream all_of;


    for(auto const& node: m)
    {
        //CellPopulation* t = node->second;
        all_of.open(opath + "_" + m[x]->name + ".mut", std::fstream::in | std::fstream::out | std::fstream::app);
        all_of << time << "\t";

		std::map<int, int> mut = m[x]->cells.count_map_mutation();

        for(int i = 0; i <= m[x]->numMutations; i++)
        {
            all_of << mut[i] << "\t";
        }
        all_of << std::endl;
        all_of.close();
    }
}
*/

void DiffTree::writeAllHeaders()
{

    std::ofstream all_of;


    for(auto const& x: bfs)
    {
        //CellPopulation* t = node->second;
        all_of.open(opath + "_" + m[x]->name + "_census.csv", std::fstream::in | std::fstream::out | std::fstream::app);
		all_of << "time,barcode,mutation,fitness,count" << std::endl;
        //m[x]->writeToFile(all_of, time);
        all_of.close();
    }
}

void DiffTree::writeAll(int time)
{

    std::ofstream all_of;


    for(auto const& x: bfs)
    {
        //CellPopulation* t = node->second;
        all_of.open(opath + "_" + m[x]->name + "_census.csv", std::fstream::in | std::fstream::out | std::fstream::app);
		//all_of << "time,barcode,mutation,fitness,count" << std::endl;
        m[x]->writeToFile(all_of, time);
        all_of.close();
    }
}

void DiffTree::writeMutant(std::string pop, int mutant, double fitIncrease){
	std::ofstream all_of;
	all_of.open(opath + "_mut.csv", std::fstream::in | std::fstream::out | std::fstream::app);
	all_of << mutant << "," << simTime << "," << pop << "," << fitIncrease << std::endl;
	all_of.close();
}


// Other Methods

void DiffTree::bfsPrint()
{
    std::deque<GrowingPop*> deque;
    deque.push_back(root);

    while(deque.size() > 0)
    {

        GrowingPop* node = deque.front();
        deque.pop_front();
        std::cout << *(node) << std::endl;
        deque.insert( deque.end(), node->children.begin(), node->children.end() );
    }
}

void DiffTree::time_step()
{
	// Create a deque of CellPopulations to reset all upstream events to 0
    std::deque<GrowingPop*> de;
    de.push_back(root);
    while(de.size() > 0)
    {
        GrowingPop* node = de.front();
        de.pop_front();
        node->setUpstreamEvents(0);
        de.insert( de.end(), node->children.begin(), node->children.end() );
    }

	// Create a deque of CellPopulations to iterate through tree and push back root
    std::deque<GrowingPop*> deque;
    deque.push_back(root);
    //std::cout << "Made it to here in time_step()" << std::endl;
	// Go through tree and call time_step on each population (BFS)
    while(deque.size() > 0)
    {
        GrowingPop* node = deque.front();
        deque.pop_front();
		//if(node->name == "A1")
        node->time_step();
        deque.insert( deque.end(), node->children.begin(), node->children.end() );
    }
   // std::cout << "End of time_step()" << std::endl;
}

void DiffTree::print()
{
    root->printPretty(0);
}

void DiffTree::initializeCells()
{
    long int o = 0;
    std::deque<GrowingPop*> deque;
    deque.push_back(root);

    while(deque.size() > 0)
    {

        GrowingPop* node = deque.front();
        deque.pop_front();
        //std::cout << *(node) << " " << o << std::endl;
        node->initializeCells(o);
        o = o + node->cells.total;
        deque.insert( deque.end(), node->children.begin(), node->children.end() );
    }
}

void DiffTree::initializeTree()
{
    for(auto const& x: m)
    {
		//std::cout << "Initializing CP " << x.second->name << std::endl;
		x.second->setTree(this);
        x.second->initializeTree();
    }
}



void DiffTree::printCells()
{
    //long int o = 0;
    std::deque<GrowingPop*> deque;
    deque.push_back(root);

    while(deque.size() > 0)
    {

        GrowingPop* node = deque.front();
        deque.pop_front();
        std::cout << *node << std::endl;
        node->printCodes();
        //o = o + node->numCells;
        deque.insert( deque.end(), node->children.begin(), node->children.end() );
    }
}

// METHOD TO GO THROUGH AND CALCULATE
void DiffTree::calcDelta(){
	std::cout << "Adjusting net proliferation for fixed setting..." << std::endl;
for(auto const& x: bfs)
    {
		std::cout << m[x]->name << ":\tOld Death: " << m[x]->death <<std::endl;
		//std::cout << "HERE!" << std::endl;
		 double d = m[x]->calcDelta();
		 if(d >= 0){
			m[x]->death = d;
		 } else {
			 std::cout << "\t" << m[x]->name << " needs to increase alpha by " << -d + m[x]->death << std::endl;
			 m[x]->alpha += (-d + m[x]->death);
		 }
		//std::cout << m[x]->death << std::endl;
    }
}

// METHOD TO GO THROUGH AND CALCULATE
void DiffTree::calcAlpha(){
for(auto const& x: m)
    {
		std::cout << x.second->name << ":\tOld Alpha: " << x.second->alpha <<std::endl;
		//std::cout << "HERE!" << std::endl;
		 double a = x.second->calcAlpha();
		 if(a >= 0){
			x.second->alpha = a;
		 } else {
			 std::cout << x.second->name << " needs to increase death: " << -a + x.second->alpha << std::endl;
			 x.second->death += (-a + x.second->alpha);
		 }
		//std::cout << x.second->death << std::endl;
    }
}

// Calculate the rates to certain populations
double DiffTree::dediffTo(std::string pop){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* from = x.second;
		if(from->name != pop){
			for(size_t i = 0; i < from->dediff_children.size(); i++){
				if(from->dediff_children[i]->name == pop){
					sumRate = sumRate + from->zeta[i]*from->numCells;
				}
			}
		}
    }

	return sumRate;
}

double DiffTree::getDediff(std::string from, std::string to){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					sumRate = sumRate + frompop->zeta[i]*frompop->cells.total;
				}
			}
		}
    }

	return sumRate;
}

void DiffTree::setDediff(std::string from, std::string to, double r){
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					frompop->zeta[i] = r;
				}
			}
		}
    }
}

double DiffTree::gamma1To(std::string pop){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* from = x.second;
		if(from->name != pop){
			for(size_t i = 0; i < from->children.size(); i++){
				if(from->children[i]->name == pop){
					sumRate = sumRate + from->gamma1[i]*from->cells.total;
				}
			}
		}
    }

	return sumRate;
}

double DiffTree::getGamma1(std::string from, std::string to){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					sumRate = sumRate + frompop->gamma1[i]*frompop->cells.total;
				}
			}
		}
    }

	return sumRate;
}

void DiffTree::setGamma1(std::string from, std::string to, double r){
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					frompop->gamma1[i] = r;
				}
			}
		}
    }
}

double DiffTree::gamma2To(std::string pop){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* from = x.second;
		if(from->name != pop){
			for(size_t i = 0; i < from->children.size(); i++){
				if(from->children[i]->name == pop){
					sumRate = sumRate + 2*from->gamma2[i]*from->numCells;
				}
			}
		}
    }

	return sumRate;
}

double DiffTree::getGamma2(std::string from, std::string to){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					sumRate = sumRate + 2*frompop->gamma2[i]*frompop->cells.total;
				}
			}
		}
    }

	return sumRate;
}

void DiffTree::setGamma2(std::string from, std::string to, double r){
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					frompop->gamma2[i] = r;
				}
			}
		}
    }
}

double DiffTree::betaTo(std::string pop){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* from = x.second;
		if(from->name != pop){
			for(size_t i = 0; i < from->children.size(); i++){
				if(from->children[i]->name == pop){
					sumRate = sumRate + from->beta[i]*from->numCells;
				}
			}
		}
    }

	return sumRate;
}

double DiffTree::getBeta(std::string from, std::string to){
	double sumRate = 0;
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					sumRate = sumRate + frompop->beta[i]*frompop->cells.total;
				}
			}
		}
    }

	return sumRate;
}

void DiffTree::setBeta(std::string from, std::string to, double r){
 for(auto const& x: m)
    {
		GrowingPop* frompop = x.second;
		if(frompop->name == from){
			for(size_t i = 0; i < frompop->children.size(); i++){
				if(frompop->children[i]->name == to){
					frompop->beta[i] = r;
				}
			}
		}
    }

}

void DiffTree::setFitnessDistribution(){
	if(fp.fitness_distribution == "doubleexp"){
		this->ConstantGenerateFitness = &cdoubleexp;
	} else if(fp.fitness_distribution == "normal"){
		this->ConstantGenerateFitness = &cnormal;
	} else if(fp.fitness_distribution == "uniform"){
		this->ConstantGenerateFitness = &cuniform;
	}
}
