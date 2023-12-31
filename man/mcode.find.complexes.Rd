% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/community_detection.R
\name{mcode.find.complexes}
\alias{mcode.find.complexes}
\title{MCODE Find Complexes Procedure}
\usage{
mcode.find.complexes(graph, neigh, vertex_weight, D)
}
\arguments{
\item{graph}{An `igraph` object representing the input network.}

\item{neigh}{A list where each element contains the neighbors of each vertex in the graph.}

\item{vertex_weight}{A numeric vector containing the weights for each vertex in the graph.}

\item{D}{A numeric value representing the density parameter used in the identification of complexes.}
}
\value{
A list with two components:
        - `COMPLEX`: A list where each element is a numeric vector of vertex indices forming an identified complex.
        - `seen`: A numeric vector marking the vertices that have been processed.
}
\description{
This function discovers multiple densely connected clusters (complexes) within
a given graph. It iterates over all vertices, utilizing the `mcode.find.complex`
helper function, to efficiently uncover all potential complexes based on vertex
weights and a specified density parameter.
}
\examples{
\dontrun{
# Assuming that you have a graph and its related information like neigh, vertex_weight, and D.

result <- mcode.find.complexes(graph, neigh, vertex_weight, D)
print(length(result$COMPLEX)) # Output the number of identified complexes
}
}
