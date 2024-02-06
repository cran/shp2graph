\name{ME.simplification}
\alias{ME.simplification}
\title{Simplify multiple edges in a network}
\description{
A function to simplfy multiple-edge into one representative edge, where a 
multiple-edge structure refers to a set of edges sharing with the same pair of nodes.
}
\usage{
ME.simplification(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, 
                  DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL,Nexception=NULL,
                  Eexception=NULL)
}
\arguments{
  \item{nodelist}{A \dQuote{\link{nodelist}} object;}
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{eadf}{ A data frame of attributes corresponding to all the edges;}
  \item{ea.prop}{A vector (of the length equalling to the number of edge attributes in
  \dQuote{eadf}) with values of 1, 2, 3 or 4, and each specific
  value represent a function to assign the attributes of the edited edges, see also \link{Redef.functions}:
  
  1->sum(v)
  
  2->min(v)
  
  3->max(v)
  
  4->mean(v)
  }
                 
  \item{Directed}{TRUE if edges are directed, FALSE otherwise;}  
  \item{DegreeL}{An integer vector of degrees for each node in the given
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
  \item{newNodelist}{New node list with multiple-edges simplified;}
  \item{newEdgelist}{New edge list with multiple-edges simplified;}
  \item{newEadf}{New attribute data frame for the returned edgelist;}
  \item{DegreeL}{New degree vector cooresponding to the newly returned node list;}
  
For \dQuote{directed} network:
  \item{newNodelist}{New node list with multiple-edges simplified;}
  \item{newEdgelist}{New edgelist with multiple-edges simplified;}
  \item{newEadf}{Data frame of attributes for the newly returned edge list;}
  \item{InDegreeL}{New In-degree vector cooresponding to the newly returned node list;}
  \item{OutDegreeL}{New Out-degree vector cooresponding to the newly returned node list;}   
}
\note{
Edges are recognised as a structure of multiple-edge when they share the same pair
of nodes. Note that the order of the nodes should be also the same when edges are
directed.
}
\seealso{\code{\link{SL.extraction}}, \code{\link{PN.amalgamation}}}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}

\keyword{simplified graph}
