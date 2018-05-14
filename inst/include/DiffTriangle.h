/*
 * =====================================================================================
 *
 *       Filename:  DiffTriangle.h
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include "GrowingPop.h"
#include "Node.h"

// class derived from GrowingPop
class DiffTriangle : public GrowingPop {
public:
	// Members
	int num_levels;
	int mfactor;
	std::vector<DiffTriangle*> triChildren;
	std::vector<std::vector<Node*> > ptrs;

	bool updated;

	// Constructors / Destructor - same as base CellType
	DiffTriangle();
	DiffTriangle(std::string cname);
	DiffTriangle(std::string cname, int num_levels, int mfactor);
	//DiffTriangle(std::string cname, int num_levels, int mfactor, long int offset);
	//DiffTriangle(std::string cname, int num_levels, int mfactor, long int offset, double eff);
	~DiffTriangle();

	// Set methods
	//virtual std::vector<DiffTriangle*> getChildren();

	// Child methods
	void addTriChild(DiffTriangle* child);
	void addTriChild(DiffTriangle* child, double diff_rate);

	// Initialize a DiffTriangle
	virtual void initialize(long int offset = 0);

	// overridden method to add a cell, initiates differentiation cascade
	virtual void add_cell(long int code, std::string m, double f);


	// NodeList manipulations
	NodeList getRemoveLastLevel(int tri);
	void doubleCells(int tri);
	NodeList removeLastAndDouble(int tri);

	// Output a Triangle
	void printCodes();

	// Methods for the DiffTree to use to simulate
	double getSumRates();
	virtual std::vector<double> getRates();
	void time_step(bool verbose = false);

	// Overridden methods
	double sdi();
	double coded();
	
	virtual void writeToFile(std::ofstream& of);
};

// Old helper methods
int beginIndex(int mfac, int level);
std::vector<std::vector<long int> > splitDoubleVector(std::vector<long int> v, std::vector<double> probs);

