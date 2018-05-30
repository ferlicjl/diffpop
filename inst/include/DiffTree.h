/*
 * =====================================================================================
 *
 *       Filename:  DiffTree.h
 *
 *    Description:  Class for calling functions throughout a hierarchy
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>

#include "GrowingPop.h"
#include "constantRVFunctions.h"


class DiffTree {
public:
	// Members
	GrowingPop* root;
	std::string opath;
	int mutationCounter;

	// current simulation time
	double simTime;

	// Map from population name to object instances to all Populations
	std::map<std::string, GrowingPop*> m;
	
	// Vector to hold BFS order of populations
	std::vector<string> bfs;

	// Constructors / Destructor
	DiffTree();
	DiffTree(std::string outdir);
	DiffTree(GrowingPop* root);
	DiffTree(GrowingPop* root, std::string o);
	DiffTree(GrowingPop* root, std::string o, std::map<std::string, GrowingPop*> m_);
	DiffTree(std::string o, std::map<std::string, GrowingPop*> m_);
	~DiffTree();

	// Set
	void setRoot(GrowingPop* r);
	void zeroEvents();
	void setMap(std::map<std::string, GrowingPop*> m_);

	// Get
	std::vector<double> getRates();
	std::vector<GrowingPop*> getCellTypes();
	bool checkZero();

	//Mutation methods
	int getNextMutant();
	void writeMutant(std::string pop, int mutant, double fitIncrease);

	// Gillespie Methods
	double nextEventTime();
	void updateState(double currentTime, double timeToNext);

	// Simulate
	void simulate(int numTime, int outputFrequency = -1);
	void time_steps(int n, int outputFrequency = -1);

	// Output methods
	void writeCellHeader(std::ofstream& of);
	void writeDiversityHeader(std::ofstream& of);
	void writeTypesHeader(std::ofstream& of);
	void writeEventsHeader(std::ofstream& of);
	void writeMutantsHeader(std::ofstream& of);
	void writeSDI(std::ofstream& of, int time);
	void writePopSize(std::ofstream& of, int time);
	void writeLabelled(std::ofstream& of, int time);
	void writeEvents(std::ofstream& of, int time);
	void writeTypes(std::ofstream& of, int time);
	//void writeMutations(int time);
	void writeAll(int time);

	// Methods
	void bfsPrint();
	void time_step();
	void print();
	void initializeCells();
	void initializeTree();
	void printCells();

	// Methods to help calculating/adjusting net proliferation in FixedPops
	void calcDelta();
	void calcAlpha();
	double dediffTo(std::string pop);
	double gamma1To(std::string pop);
	double gamma2To(std::string pop);
	double betaTo(std::string pop);
	double getDediff(std::string from, std::string to);
	double getGamma1(std::string from, std::string to);
	double getGamma2(std::string from, std::string to);
	double getBeta(std::string from, std::string to);
	void setDediff(std::string from, std::string to, double r);
	void setGamma1(std::string from, std::string to, double r);
	void setGamma2(std::string from, std::string to, double r);
	void setBeta(std::string from, std::string to, double r);

	// Fitness Distribution Stuff
	FitnessParameters fp;
	void (*ConstantGenerateFitness)(double *, struct FitnessParameters*, gsl_rng*);
	void setFitnessDist(FitnessParameters f);
	void setFitnessDistribution();

};

