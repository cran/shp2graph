\name{ME.simplification}
\alias{ME.simplification}
\title{Simplify multiple edges in a network}
\description{
This function simplfies multiple-edge into one representative edge, where
multiple-edge refers to a set of edges sharing with the same pair of nodes.
}
\usage{
ME.simplification(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, 
                  DegreeL=NULL, InDegreeL=NULL, OutDegreeL=NULL,Nexception=NULL,
                  Eexception=NULL)
}
\arguments{
  \item{nodelist}{A \dQuote{nodelist} object}
  \item{edgelist}{An \dQuote{edgelist} object}
  \item{eadf}{Attribute data frame for all the edges, of which the number of 
              rows equals to the number of edges;}
  \item{ea.prop}{a vector (of lenght equal to the number of columns in 
                 \dQuote{eadf}) consisted of 1, 2, 3 or 4, and one value 
                 cooresponds to a kind of redefination
                 function, see \link{Redef.functions};}
  \item{Directed}{A logical parameter to specify whether edges are directed or not;}  
  \item{DegreeL}{A vector (of length equal to the number of nodes) of degrees for 
                 the nodes, and it should be NULL if \dQuote{Directed} is TRUE;}
  \item{InDegreeL}{A vector (of length equal to the number of nodes) of 
                   in-degrees for the nodes;}
  \item{OutDegreeL}{A vector of out-degrees for the nodes;}
  \item{Nexception}{A vector of node IDs considered as exceptions, and all the 
                   nodes included won't be processed;}
  \item{Eexception}{A vector of edge IDs considered as exceptions, and all the 
                    edges included won't be processed;}  
}
\value{
Two types of list for \dQuote{undirected} and \dQuote{directed} edges respectively:

For \dQuote{undirected} network:
  \item{newNodelist}{The new nodelist with multiple-edges simplified;}
  \item{newEdgelist}{The new edgelist with multiple-edges simplified;}
  \item{newEadf}{The new attribute data frame for the returned edgelist;}
  \item{DegreeL}{The new degree vector cooresponding to the returned nodelist;}
  
For \dQuote{directed} network:
  \item{newNodelist}{The new nodelist with multiple-edges simplified;}
  \item{newEdgelist}{The new edgelist with multiple-edges simplified;}
  \item{newEadf}{The new attribute data frame for the returned edgelist;}
  \item{InDegreeL}{The new in-degree vector cooresponding to the returned 
                   \dQuote{nodelist};}
  \item{OutDegreeL}{The new out-degree vector cooresponding to the returned 
                   \dQuote{nodelist};}   
}
\note{
If input edges are directed, the judgement of multiple-edge will be based on its 
in-degree and out-degree equalling to 1, then InDegreeL and OutDegreeL should be 
specified; Otherwise, the criterion will be the degree equalling to 2 and then 
DegreeL is the significant.
}
\seealso{\code{\link{SL.extraction}}, \code{\link{PN.amalgamation}}}
\author{Binbin Lu \email{lubinbin220@gmail.com}}

\keyword{spatial, graph}
