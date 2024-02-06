\name{Directed}
\alias{Directed}
\title{Orientate all the edges in a given edgelist}
\description{
A function to orientate each edge according to the given vector.
}
\usage{
Directed(edgelist, direction.v=rep(0,length(edgelist[,1])), eadf=NULL)
}
\arguments{
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{direction.v}{A vector (of the length equalling to the number of edges in
  the given \dQuote{\link{edgelist}}) with values of 1 (TRUE) or 0 (FALSE), 1 (TRUE)
  indicates a directed (one-way) edge, while 0 (FALSE) means an undirected (double-way) edge}
  \item{eadf}{ A data frame of attributes corresponding to all the edges;}                             
}
\value{
A list consisted of:
  \item{newEdgelist}{An \dQuote{\link{edgelist}} object with directed edges}
  \item{newEadf}{A data frame of attributes for the new \dQuote{\link{edgelist}}}
}
\details{
Within a road network, some road segments might be one-way while the rest 
are double-way. This suitation makes it complex to define directed or undirected
 edges in a graph. This function is to orientate each edge according to the given 
 vector, \dQuote{direction.v}: 1 (TRUE) indicates one-way, while 0 (FALSE) represents 
 double-way. All the double-way (undirected) edges in the given \dQuote{\link{edgelist}} are 
 redefined as two directed edges (e.g. a double-way (undirected) edge(nid1, nid2) is redefined 
 as two one-way (directed) edges, (nid1, nid2) and (nid2, nid1)). The one-way (directed) edges 
 are left in their directed forms. In other words, all the edges returned by this function will 
 be directed, and used to be construct a directed graph.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{directed graph}
