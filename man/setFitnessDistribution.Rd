% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{setFitnessDistribution}
\alias{setFitnessDistribution}
\title{setFitnessDistribution}
\usage{
setFitnessDistribution(tree, distribution = "normal", alpha_fitness = 0,
  beta_fitness = 1, pass_prob = 1, upper_fitness = NA,
  lower_fitness = 0)
}
\arguments{
\item{tree}{DiffTree object for which to set fitness distribution}

\item{distribution}{random distribution from which to draw ("normal", "doubleexp", "uniform")}

\item{alpha_fitness}{alpha parameter for fitness distribution}

\item{beta_fitness}{beta parameter for fitness distribution}

\item{pass_prob}{probability that new clone developed a passenger mutation that does not affect its fitness}

\item{upper_fitness}{upper limit on clone fitness}

\item{lower_fitness}{lower limit on clone fitness}

\item{is_random}{boolean variable set to TRUE if random fitness distribution is specified}
}
\description{
Used to designate distribution from which changes in fitness for a new clone are drawn
}
\examples{
\dontrun{
setFitnessDistribution(tree = myTree, distribution = "uniform", alpha_fitness = 0, beta_fitness = 1, pass_prob = 0, upper_fitness = 5, lower_fitness = 0)
}
}
