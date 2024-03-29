% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{est_rbsvar}
\alias{est_rbsvar}
\title{Estimates rbsvar model by sampling from the posterior using DE-MC algorithm (Documentation incomplete)}
\usage{
est_rbsvar(
  model,
  N,
  n = 2,
  K = 2,
  m0 = NA,
  rel_gamma = 1,
  output = NULL,
  new_chain = FALSE,
  parallel_chains = FALSE,
  parallel_likelihood = FALSE,
  max_cores = NA,
  progress_bar = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{model}{A list outputted by \code{init_rbsvar()}.}

\item{N}{...}

\item{n}{...}

\item{K}{...}

\item{m0}{...}

\item{rel_gamma}{...}

\item{output}{...}

\item{new_chain}{...}

\item{parallel_chains}{...}

\item{parallel_likelihood}{...}

\item{max_cores}{...}

\item{progress_bar}{...}

\item{verbose}{...}
}
\value{
A list that containing (among other things) the chain(s) sampled from the posterior. May be passed for example to \code{irf()}.
}
\description{
Estimates rbsvar model by sampling from the posterior using DE-MC algorithm (Documentation incomplete)
}
