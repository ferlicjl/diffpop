% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{GrowingPop}
\alias{GrowingPop}
\title{GrowingPop}
\usage{
GrowingPop(tree, name, size, label = 0)
}
\arguments{
\item{tree}{DiffTree object to add GrowingPop to}

\item{name}{population name}

\item{size}{initial population size}

\item{label}{initial probabilty that a cell receives a unique barcode}
}
\description{
Designates a GrowingPop, a structure used to simulate differentiation via a branching process
}
\details{
A GrowingPop is the base class used to designate the various cell types throughout a differentiation tree.  A GrowingPop contains a list of cell states, functions to enact cellular events on those cell states, and event rates at which to perform those functions.  The hierarchical structure is maintained by pointers to upstream and downstream CellPopulations.
}
\examples{
\dontrun{
GrowingPop(myTree, "Population A", 1000, 0.5)
}
}
