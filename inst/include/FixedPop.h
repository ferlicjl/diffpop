/*
 * =====================================================================================
 *
 *       Filename:  FixedPop.h
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

#include <string>
#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include "GrowingPop.h"

// Class derived from GrowingPop
class FixedPop : public GrowingPop {
public:

	// Constructors / Destructor - same as base CellType
	FixedPop();
	FixedPop(std::string cname);
	FixedPop(std::string cname, long int ncells);
	FixedPop(std::string cname, long int ncells, long int offset);
	FixedPop(std::string cname, long int ncells, long int offset, double eff);
	~FixedPop();


	// overriden methods to yield fixed population size
	void dediff();
	virtual void time_step(bool verbose = false);
	void add_cell(long int code, std::string mut, double f);
	virtual void single_event(bool verbose = false);
	
	// Up/Down events to call
	void up_Event(bool verbose = false);
	void down_Event(bool verbose = false);
	
	double getNumEvents();
};

