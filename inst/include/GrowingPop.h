/*
 * =====================================================================================
 *
 *       Filename:  GrowingPop.h
 *
 *    Description:  Class for an exponentially growing population
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <vector>
#include <map>

#include "NodeList.h"

class DiffTree;

class GrowingPop {
public:
	// Members
	std::string name;
	NodeList cells;
	int numChildren;
	int numCells;
	int upstreamEvents;
	int alphaEvents;
	int gamma1Events;
	int gamma2Events;
	int deathEvents;
	int betaEvents;
	int dediffEvents;
	int addcellEvents;
	int mutateEvents;

	// Mutation information
	int numMutations;
	int numMutCells;
	int maxMutation;

	// Initial barcode efficiency - initial fraction of cells receiving unique barcode
	double efficiency;


	// Pointers to other Cell Types
	GrowingPop* parent;
	std::vector<GrowingPop*> children;
	std::vector<GrowingPop*> dediff_children;

	// Pointer to tree
	DiffTree* tree;

	// Rates/Probabilities
	double alpha; // self-proliferation rate
	std::vector<double> beta; // Asymmetric differentiation rates
	std::vector<double> gamma2; // Symmetric differentiation rates
	std::vector<double> gamma1; // 1-to-1 differentiation rates
	std::vector<double> zeta; // De-differentiation to an upstream type
	double death; // apoptosis rate
	double mu; // mutation rate

	// Temporary structures used to load in differentiation rates
	std::map<std::string, std::vector<double>> loadDiffRates;
	std::map<std::string, GrowingPop*> c_map;

	// Constructors / Destructor
	GrowingPop();
	GrowingPop(std::string cname);
	GrowingPop(std::string cname, long int ncells);
	GrowingPop(std::string cname, long int ncells, long int offset);
	GrowingPop(std::string cname, long int ncells, long int offset, double eff);
	~GrowingPop();

	// Child methods
	void addChild(GrowingPop* child, double diff1_rate = 0.0, double diff2_rate = 0.0, double asym_rate = 0.0);

	// Dedifferentiation
	void addDediff(GrowingPop* to, double rate);
	
	// Set methods
	void setParent(GrowingPop* p);
	void setTree(DiffTree* t);
	void setAlpha(double a);
	void setDeath(double d);
	void setMu(double m);
	void setUpstreamEvents(int x);
	void zeroEvents();
	void setDiff1(GrowingPop* child, double r);
	void setDiff2(GrowingPop* child, double r);
	void setAsymdiff(GrowingPop* child, double r);

	// Calculates and adjusts net proliferation in FixedPops, coded for all populations however
	double calcDelta();
	double calcAlpha();

	// Print methods
	friend std::ostream& operator<<(std::ostream& os, const GrowingPop& c);
	void printPretty(int index);
	virtual void printCodes();
	void print_rates();

	// Get methods
	virtual std::vector<double> getRates();
	virtual double getNumEvents();
	virtual std::vector<double> getDownRates();
	virtual double getNumDownEvents();
	virtual std::vector<GrowingPop*> getChildren();

	// Cell manipulations
	void mitosis();
	void asym_diff();
	void apoptosis();
	void diff1();
	void diff2();
	virtual void dediff();
	void mutate(long int index);
	void mutateCur();
	void erase_cell(long int index);
	virtual void add_cell(long int code, std::string mut, double f);
	void add_dediffCell(long int code, std::string mut ,double f);

	// Initializations
	virtual void initializeCells(long int offset = 0);
	virtual void initializeTree();
	
	// Simulate a time-step
	virtual void time_step(bool verbose = false);
	virtual void single_event(bool verbose = false);

	// Statistics on cell population
	virtual double sdi();
	std::string diversity();
	virtual double coded();
	std::vector<long int> sample(int size);
	
	// Full output of population
	virtual void writeToFile(std::ofstream& of, int time);

};

// Helper methods
long int cellByFitness(std::vector<double> fitness);
std::vector<double> normalize(std::vector<double> input);
int choose(std::vector<double> input);
std::map<long int, int> count_map(std::vector<long int> input);
std::map<int, int> count_map(std::vector<int> input);
void printVec(std::vector<long int> v);
void printVec(std::vector<double> v);
