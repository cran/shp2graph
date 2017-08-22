\name{Degree.list}
\alias{Degree.list}
\title{Degree (In-degree and Out-degree) of nodes}
\description{
A function that returns degrees of nodes from provided \dQuote{\link{nodelist}} and \dQuote{\link{edgelist}}; 
while in-degrees and out-degrees are returned if edges are directed).
}
\usage{
Degree.list(nodelist, edgelist, Directed=F)
}
\arguments{
 \item{nodelist}{A \dQuote{\link{nodelist}} object}
 \item{edgelist}{An \dQuote{\link{edgelist}} object}
 \item{Directed}{TRUE if edges are directed; FALSE, otherwise} 
}
\value{
Lists of different contents are returned for undirected and directed edges 
respectively:
For undirected graph:
  \item{DegreeL}{An integer vector of degrees for each node in the given
                \dQuote{\link{nodelist}}}
For directed type:
  \item{InDegreeL}{An integer vector of In-degrees for each node 
                  in the given \dQuote{\link{nodelist}}}
  \item{OutDegreeL}{An integer vector of Out-degrees for each node 
                    in the given \dQuote{\link{nodelist}}}  
}
\note{
This function returns differently for undirected and directed 
networks, where both In-degrees and Out-degrees are calculated for a directed network,
and only degrees are returned for an undirected networks.
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{degree, node}
