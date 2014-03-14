\name{Directed}
\alias{Directed}
\title{Orientate each edge in a given edgelist}
\description{
This function is to orientate each edge according to the given vector.
}
\usage{
Directed(edgelist, direction.v=rep(0,length(edgelist[,1])), eadf=NULL)
}
\arguments{
  \item{edgelist}{An edgelist object, see \dQuote{edgelist}}
  \item{direction.v}{A vector (of length equal to number of edges in edgelist)
                    with values of 1 or 0, see details below;}
  \item{eadf}{Attribute data frame for all the edges;}                             
}
\value{
A list consisted of:
  \item{newEdgelist}{A new \dQuote{edgelist} with directed edges}
  \item{newEadf}{A new attribute data frame for the new \dQuote{edgelist}}
}
\details{
In a road network data set, some road segments might be one-way while the rest 
are double-way. This suitation makes it complex for defining edges in the graph.
This property could be specified in \dQuote{direction.v}: 1 means one-way and 0 
represents double-way. All the double-way edges in edgelist are clarified as two 
directed edges (i.e. edge(nid1, nid2) and edge(nid1, nid2) and edge(nid2, nid1)), 
and duplicate the corresponding attributes. The one-way edges could be regarded 
as directed.
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{directed, graph}
