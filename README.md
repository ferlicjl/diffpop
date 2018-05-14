---
title: "DIFFpop"
author: "Jeremy Ferlic, Jiantao Shi, Thomas O. McDonald, and Franziska Michor"
date: "April 18, 2018"
output: pdf_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# DIFFpop

DIFFpop is an R package to simulate cell labeling experiments performed on differentiation hierarchies.  


## Dependencies

* [GNU Scientific Library](https://www.gnu.org/software/gsl/)
    + (OSX) `brew install gsl` with Homebrew or from
  [here](http://ftpmirror.gnu.org/gsl/).
    + (Windows) download and extract
  the file [local###.zip](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/)
  and create an environmental variable LIB_GSL to add the directory (see notes
  about Windows installation below for more details).
    + (Linux) install libgsl0-dev and gsl-bin.
* [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (*Windows only*) 
* [devtools](https://github.com/hadley/devtools)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
    + Issues arise with installation of recursive dependencies using
    install_github. Installing this R package first solves this issue.
    
### Important Notes about Windows installation
Rtools contains the necessary resources to compile C++ files when installing
packages in R. GSL is also required which can be downloaded from [here](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip).
After downloading, unzip to your R directory. Depending on if your computer is
32 or 64-bit, move the library files from __local###/lib/i386__ (32-bit) or
__local###/lib/x64__ (64-bit) to __local###/lib__.

To set the environmental variable LIB_GSL on a Windows 7/10 computer, go to
"Advanced system settings" in *Control Panel > System and Security > System*
and click *Environmental Variables*. Create a new system variable with

* Variable Name: __LIB_GSL__
* Variable Value: __"C:/path/to/local323"__ (include quotes)

## Recommended R packages
The following R packages are required for certain functions.

* igraph
* ggplot2
* R.utils

# Installation

To install in R, type:

```r
install.packages("devtools")
devtools::install_git("https://github.com/ferlicjl/diffpop.git")
```

Installing the library should compile all necessary functions so that DIFFpop
can be run as an R function.

# Uses

DIFFpop is an R package that uses C++ to simulate the maturation process of cells through user-defined differentiation hierarchies according to user-defined event rates.  The software is capable of simulating hierarchies in two manners: using a branching process with growing populations and using a modified Moran process with fixed population sizes.

For increasing population simulations, the system will grow or decline as a branching process according to the exact rates specified by the user.  If at any point any population becomes extinct, the simulation halts. A branching process is a stochastic process used to model the growth and composition of reproducing populations. Assumptions made in branching processes are individuals live for a random amount of time before some event. Here, we will only consider events that are analogous to cellular differentiation events, i.e. mitosis, differentiation, and apoptosis. Individuals of the same type are independent and identically distributed.  Simulation of the branching process is an application of the direct Gillespie Stochastic Simulation Algorithm.

Simulations with fixed population sizes are carried out using a modified Moran Process.  A traditional Moran Process consist of a population of cells, each belonging to a certain class, and having an individual fitness.  In each time step, a random cell is selected for division and a random cell is selected for apotosis, ensuring that the population size remains constant.  To introduce selection, cells with higher fitness values are more likely to be selected for division.  To modify the Moran Process and place it in the setting of differentiation, we have included additional cellular events.  By continuing to couple events together, we are able to maintain a constant population size.

# Software Design and Class Structures

## GrowingPop

A GrowingPop is the base class used to designate the various cell types throughout a differentiation tree.  A GrowingPop contains a list of cell states, functions to enact cellular events on those cell states, and event rates at which to perform those functions.  The hierarchical structure is maintained by pointers to upstream and downstream CellPopulations.

## FixedPop

A FixedPop is a class derived from a GrowingPop.  In order to maintain a constant population size, cellular events are coupled together; i.e. a mitosis event generating an additional cell is immediately coupled with a differentiation or death event.  Similarly, if the number of cells in the FixedPop population increases by one from upstream differentiation, a differentiation or death event of its own is immediately enacted to maintain the population level.  

## DiffTriangle

A DiffTriangle cell type is used to represent the downstream fully differentiated cells.  Cells are arranged in a triangle formation.  Cells enter the population on the highest level of the triangle, experiencing further differentiation and division to progress down the triangle.  When a new cell enters the DiffTriangle population, it causes an already existing cell on the highest level to divide and further differentiate to the next level of the triangle.  Those two cells each displace a pre-existing cell, causing them to divide and differentiate (thus generating four newly displaced cells), which in turn displace cells that further displace cells until reaching the lowest level of the triangle.  A displaced cell from the last row of the triangle can either be passed on to an offspring cell type (if there are further cell types in the hierarchy), or die out.  Importantly, DiffTriangle structures will not initiate any cellular events of their own, as differentiation waves throughout a triangle is only initiated when receiving a new cell from an upstream population.

The population size of a DiffTriangle is specified by two parameters; the first, $z$ is the number of cell divisions until full maturation or the number of levels in the triangle.  The second is called the $mfactor$, which is the number of triangles to be stacked side-by-side.  If $mfactor$ is greater than 1, then a cell entering the DiffTriangle population simply chooses at random which specific triangle to enter.

## DiffTree

A DiffTree contains pointers to all populations throughout the differentiation hierarchy.  From the DiffTree object, all cell types in the hierarchy can be accessed in either a breadth-first or depth-first manner.  Hierarchy-wide functions, such as simulating the hierarchy and recording output, are initiated by calling functions of the DiffTree.  

 |  Class Type | Population Size | Use                                  | 
 | ----------+---------------+------------------------------------ | 
 |  GrowingPop | Dynamic | dynamically sized population with exponential event waiting times | 
 |  FixedPop | Constant | constant size homogenous population  | 
 |  DiffTriangle | Constant | constant size population with $z$ levels of maturity | 


## NodeList Data Structure

The cells of each population are maintained by a NodeList.  A NodeList is a doubly-linked list of Nodes.  Each Node keeps track of a particular cell state, defined by the combination of a barcode, mutation status, and fitness value, as well as a count of how many cells belong to that particular state.  In addition, DiffTriangle Nodes also contain data to record which triangle and at which level a cell resides.  Mutation status is a string listing which mutations have occurred in that state.  Information about particular mutations can be found in the mutation information file.  The NodeList keeps track of the total number of cells in the list and overall fitness for the population.  Methods exist to insert or remove cells from the list, maintaining a left-balanced list by cell count -- that is, cell states with higher cell counts are found to the left of the list.  This approach allows for more efficient indexing by cell count, particularly as diversity in the compartment decreases or dominant clones arise. 

\begin{center}
\includegraphics{nodelist1}\

NodeList automatically balances by count:\newline
    \begin{verbatim}insert(barcode 2, mutation 0, fitness 1, count 15)\end{verbatim}

\includegraphics{nodelist2}\
\end{center}

#Event Types and Parameters

Cellular events in DIFFpop are enacted according to their accompanying parameter rates, in units of number of events per cell per time unit.  

 | Parameter	 | Variable Type	 | Description                                                    | 
 | -----------+---------------+-------------------------------------------------------------- | 
$\alpha$ (alpha) | 	double | 	mitotic self-renewal rate | 
$\beta$ (beta) |  double | 	asymmetric differentiation rate to downstream cell type | 
$\gamma_1$ (gamma1) | 	double | 	mitosis-independent (one-to-one) differentiation rate to downstream cell type | 
$\gamma_2$ (gamma2) | 	double | 	mitosis-dependent (one-to-two) differentiation rate to downstream cell type | 
$\zeta$ (zeta) | 	double | 	de-differentiation rate to upstream cell type | 
$\delta$ (delta) | 	double | 	apoptosis rate | 
$\mu$ (mu) | 	double | 	probability of mutation per mitotic event | 
*Note Population specific parameters are specified using subscripts, i.e. $\alpha_{(LT-HSC)}$ is the mitotic self renewal rate for the LT-HSC population and $\gamma_{1(LT-HSC,ST-HSC)}$ is the one-to-one differentiation rate from the LT-HSC population to the ST-HSC population.


Events can be split into three categories based on how they affect the population size of the compartment: 

* “i-1” events result in a one-cell deficit
    + differentiation ($\gamma_1$/$\gamma_2$)
    + de-differentiation ($\zeta$) 
    + apoptosis ($\delta$)
* “i” events maintain the population size
    + asymmetric differentiation ($\beta$)
* “i+1” events result in a one-cell surplus
    + mitosis ($\alpha$)

Mutations in DIFFpop occur only during mitosis events.  Each mitosis event results in a new mutation with probability $\mu$. Thus, in the MPP population, the rate of mitosis events resulting in a new mutation is $\alpha_{(MPP)}\mu_{(MPP)}$ and the rate of mitosis events resulting in no mutation is $\alpha_{(MPP)}(1-\mu_{(MPP)})$.  Mutations accumulate according to the infinite allele assumption, that is, a new mutation leads to a new allele that has yet to be seen in the population.

In the FixedPop setting, all de-differentiation events will result in an apoptosis event in the receiving population.  This is necessary to avoid circular equations when calculating adjustments to net proliferation in order to maintain a constant population size.  


# Simulation Overview

After the tree structure has been specified in R, the task of simulating is handed off to a C++ backend.  Before events can be enacted, the tree structure and other user-specified parmaeters must be read from the input files. 

## Initializing a Simulation

1.	Read tree structure from input files
2.	Initialize cell list for each population
3.	Write information file
4.	Output tree structure

## Simulation via Gillespie Algorithm

1.	Initiate simulation at time 0
2.	Generate the time to the next event
    a.	Iterate through population hierarchy, collecting total event rate from each population
    b.	Generate a time from exponential distribution, parameter equal to sum of event rates from all populations
3.	Make observations
4.	Update state
    a.	Choose which population enacts event
    b.	Enact single event on that population
5.	If simulation time remains, repeat from step 2

## Simulation via modified Moran Process

1.	For each time unit:
    a.	Iterate through tree in breadth-first manner
    b.	For each population:
        i.	Get total number of “i-1” and “i” events [$\gamma_1$, $\gamma_2$, $\delta$, $\zeta$, and $\beta$] = $j$
        ii.	Enact $j$ events choosing which event according to the event rates. If event results in a one-cell deficit, enact a mitosis [$\alpha$] event

## Output Files

Each run of the simulation will be given a unique file prefix, consisting of the date and time when the simulation initiated, followed by a random integer. 

1.	Population size (out_31-05-2018_123456.pop)
    +	Size of each population at each observation time
2.	Shannon diversity index (out_31-05-2018_123456.sdi)
    +	Shannon diversity index is calculated at each observation time for each population
3.	Fraction of labelled cells (out_31-05-2018_123456.label)
    +	Fraction of cells that contain unique barcode/label at each observation time for each population
4.	Event rates (out_31-05-2018_123456.events)
    +	Number of events that occurred between each observation time
5.	Mutation summary (out_31-05-2018_123456.mut)
    +	Time, compartment location, and additional fitness at which a new mutation arises


# Maintaining a constant population size

In order to maintain a constant population size, a relationship must exist between the event rates of the compartment.  Specifically, those event rates that result in an excess of cells in the compartment [“i+1” rates] must be balanced with those event rates that result in a deficit of cells [“i-1” rates].  

Let $\alpha_{(x)}$ denote the mitotic self-renewal ($\alpha$) rate of population $x$. \newline
Let $\gamma_{1(x,y)}$ denote the one-to-one differentiation ($\gamma_1$) rate from population $x$ to population $y$.\newline
Let $\gamma_{2(x,y)}$ denote the one-to-two differentiation ($\gamma_2$) rate from population $x$ to population $y$.\newline
Let $\beta_{(x,y)}$ denote the asymmetric differentiation ($\beta$) rate from population $x$ to population $y$.\newline
Let $\zeta_{(x,y)}$ denote the one-to-one de-differentiation ($\zeta$) rate from population $x$ to population $y$.\newline
Let $\delta_{(x)}$ denote the cell death ($\delta$) rate of population $x$.\newline
Let $n_{(x)}$ be the size of population $x$.\newline

Then, for any compartment $A$,

\begin{align*}
n_{(A)}\alpha_{(A)} &+ \sum_{\text{pop }i \neq A}n_{(i)}\left( \gamma_{1(i,A)} + 2\gamma_{2(i,A)} + \beta_{(i,A)} + \zeta_{(i,A)}\right) = n_{(A)}\left[ \delta_{(A)} + \sum_{\text{pop }i\neq A}\left( \gamma_{1(A,i)} +\gamma_{2(A,i)} + \zeta_{(A,i)}\right) \right] \\
\implies \\
\delta_{(A)} &= \frac{n_{(A)}\left[\alpha_{(A)}- \sum_{\text{pop }i\neq A}\left( \gamma_{1(A,i)} + \gamma_{2(A,i)} + \zeta_{(A,i)} \right)\right] + \sum_{\text{pop }i \neq A}n_{(i)}\left( \gamma_{1(i,A)} + 2\gamma_{2(i,A)} + \beta_{(i,A)} + \zeta_{(i,A)}\right)}{n_{(A)}}
\end{align*}

For the population size to remain constant, this first line of the equation must hold.  That is, events that increase the population size (self-renewal and influx from other populations) must be balanced by events that decrease the population (cell death and differentiation).  In the modified Moran Process, we force this to hold by automatically calculating delta for each population.  If this calculated delta value is positive, we simply set the effective death rate equal to this value.  If this calculated delta value is negative, we increase the alpha rate of the population by this value.
 
# Fitness Distribution
	
Throughout the differentiation hierarchy, whenever a new clone arises due to mutation, a change in fitness can be drawn from a random distribution.  The parameters of that distribution can be specified by the user in R.  

 | Parameter Name | Type |  Description                                 | 
 | --------------+----+-------------------------------------------- | 
 | fitness_distribution |  string  |  "doubleexp", "normal", "uniform" | 
 | alpha_fitness |  double |  alpha parameter for fitness distribution | 
 | beta_fitness |  double  |  beta parameter for fitness distribution | 
 | pass_prob |  boolean  |  probability that mutation does not result in a fitness change | 
 | upper_fitness |  double  |  upper bound on fitness | 
 | lower_fitness |  double  |  lower bound on fitness | 

If the distribution function selected is normal, fitness additions are drawn from a $N(alpha\_fitness, beta\_fitness)$ distribution.  If the distribution function selected is uniform, fitness additions are drawn from a $U(alpha\_fitness, beta\_fitness)$ distribution.  If the distribution function selected is double exponential, alpha_fitness refers to the rate parameter of an exponetial distribution for the positive range and beta_fitness refers to the rate parameter of an exponential distribution for the negative range.


# Using DIFFpop in R

DIFFpop consists of three steps to run simulations of differentiation hierarchies.  The first step is to build and describe the tree structure using the appropriate R functions. Populations of cells are created using their specific type.  Users must give each population a unique name, as well as an initial population size.  Optionally, users may specify an initial cell barcoding frequency, this is the proportion of initial cells that receive a unique barcode.  If this parameter is not set, no unique barcodes will be created for the population.  To specify the transitions between populations, the addEdge function is used, along with the correct parameters: the initiating (from) population, the receiving (to) population, event type as a string ("alpha", "beta", "gamma1", "gamma2", "delta", "zeta", "mu"), and event rate.  For events involving only one population ("alpha", "delta", or "mu"), set that population as both the initiating and receiving population.  The last step is to specify which population is the root of the differentiation hiearchy, that is, which population is the furthest upstream, using the setRoot function in R.

After the tree has been completely specified, the writeTree function is used to create the necessary input files for the C++ backend.  

The last step is to use either the simulateTree or simulateFixedTree function to kick off the C++ backend simulation.    

## Simulation Parameters

The following parameters are available to the user in the simulateTree/simulateFixedTree function.  

 | Parameter	 | Variable Type	 | Description                                                    | 
 | -----------+---------------+-------------------------------------------------------------- | 
 | nObs | 	integer | 	number of time units to simulate | 
 | traverseFrequency |  integer |  how often to output full census of populations | 
 | indir | 	string | 	directory location of input files | 
 | outdir | 	string | 	directory location for output files | 
 | seed | 	numeric |  optional seed for the random number generator | 

Observations are made and output files updated at every integer time unit through $nObs$.  In addition, full printouts of each population are made every $traverseFrequency$ time unit.  The $indir$ directory informs the C++ backend where the input files are located for the differentiation hierarchy and $outdir$ specifies a particular directory to place all output files.  Optionally, the user can specify a numeric $seed$ for the GSL random number generator used throughout the simulation.
