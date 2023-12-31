% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/community_detection.R
\name{mcode.post.process}
\alias{mcode.post.process}
\title{MCODE Post-Processing Procedure}
\usage{
mcode.post.process(
  graph,
  vertex.weight,
  haircut,
  fluff,
  fdt = 0.8,
  set.complex.g,
  seen
)
}
\arguments{
\item{graph}{An `igraph` object representing the input network.}

\item{vertex.weight}{A numeric vector containing weights assigned to each vertex in the graph.}

\item{haircut}{Logical; determines whether to trim singly-connected nodes from the clusters.}

\item{fluff}{Logical; determines whether to expand clusters by including bordering nodes.}

\item{fdt}{A numeric value setting the threshold for the fluff procedure; default is 0.8.}

\item{set.complex.g}{A list where each element is a numeric vector representing the indices of
vertices in a detected cluster.}

\item{seen}{A numeric vector indicating which vertices have been processed or visited.}
}
\value{
A list where each element is a numeric vector of vertex indices, representing refined
        clusters after post-processing steps have been applied.
}
\description{
This function refines the detected clusters from the MCODE algorithm. It filters out smaller
clusters (size <= 2) and, based on given parameters, may expand ('fluff') clusters or trim
singly-connected nodes ('haircut') from them. This post-processing aims to improve the
relevance and robustness of the detected clusters.
}
\examples{
\dontrun{
# Assuming the existence of a graph, vertex weights, and other necessary inputs:
# graph, vertex.weight, haircut, fluff, set.complex.g, and seen.
refined_clusters <- mcode.post.process(graph, vertex.weight, haircut, fluff, fdt=0.8, set.complex.g, seen)
print(length(refined_clusters)) # Output the number of refined clusters.
}
}
