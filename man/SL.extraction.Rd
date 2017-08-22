\name{SL.extraction}
\alias{SL.extraction}
\title{Extract self-loops form given network}
\description{
A function to extract self-loops away in given network.
}
\usage{
SL.extraction(nodelist, edgelist, eadf=NULL,Directed=F, DegreeL=NULL, 
              InDegreeL=NULL, OutDegreeL=NULL,Nexception=NULL,Eexception=NULL)
}
\arguments{
  \item{nodelist}{A \dQuote{\link{nodelist}} object;}
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{eadf}{ A data frame of attributes corresponding to all the edges;}
  \item{Directed}{TRUE if edges are directed, FALSE otherwise;}  
  \item{DegreeL}{DegreeL}{An integer vector of degrees for each node in the given
                \dQuote{\link{nodelist}}, and it could be ignored if edges are directed}
  \item{InDegreeL}{An integer vector of In-degrees for each node 
                  in the given \dQuote{\link{nodelist}}, and it could be ignored 
                  if edges are undirected}
  \item{OutDegreeL}{An integer vector of Out-degrees for each node 
                  in the given \dQuote{\link{nodelist}}, and it could be ignored 
                  if edges are undirected}
  \item{Nexception}{A vector of node IDs concerned as exceptions, and all the 
                   edges with these nodes included won't be processed;}
  \item{Eexception}{A vector of edge IDs concerned as exceptions, and all these 
                    edges won't be processed;}   
}
\value{
Two types of list returned for undirected and directed edges, respectively:

For \dQuote{undirected} edges:
  \item{newNodelist}{New node list with self-loops removed;}
  \item{newEdgelist}{New edge list with self-loops removed;}
  \item{newEadf}{New attribute data frame for the returned edgelist;}
  \item{DegreeL}{New degree vector cooresponding to the newly returned node list;}
  
For \dQuote{directed} network:
  \item{newNodelist}{New node list with self-loops removed;}
  \item{newEdgelist}{New edgelist with self-loops removed;}
  \item{newEadf}{Data frame of attributes for the newly returned edge list;}
  \item{InDegreeL}{New In-degree vector cooresponding to the newly returned node list;}
  \item{OutDegreeL}{New Out-degree vector cooresponding to the newly returned node list;}  
}
\note{
An edge is recognised as a structure of self-loop when it starts and ends at the 
same point (location).
}
\seealso{\code{\link{PN.amalgamation}}, \code{\link{ME.simplification}}}
\author{Binbin Lu \email{lubinbin220@gmail.com}}

\keyword{self-loop}
