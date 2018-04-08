\name{PN.amalgamation}
\alias{PN.amalgamation}
\alias{Minus.DegreeL}
\alias{Update.nodelist}
\alias{Update.edgelist}
\alias{extend.eadf}
\alias{is.exception} 
\title{Amalgamate edges connected by a pseudo-node}
\description{
A function to amalgamate edges connected by a pseudo-node.
}
\usage{
PN.amalgamation(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, 
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
  \item{newNodelist}{New node list with pseudo-nodes removed;}
  \item{newEdgelist}{New edge list with pseudo-nodes removed;}
  \item{newEadf}{New attribute data frame for the returned edgelist;}
  \item{DegreeL}{New degree vector cooresponding to the newly returned node list;}
  
For \dQuote{directed} network:
  \item{newNodelist}{New node list with pseudo-nodes removed;}
  \item{newEdgelist}{New edgelist with pseudo-nodes removed;}
  \item{newEadf}{Data frame of attributes for the newly returned edge list;}
  \item{InDegreeL}{New In-degree vector cooresponding to the newly returned node list;}
  \item{OutDegreeL}{New Out-degree vector cooresponding to the newly returned node list;}  
}
\note{
Node in a directed network is recognised a a pseudo-node with in-degree 
and out-degree equalling to 1; while it could be a pseudo-node with degree 
equalling to 2 in a undirected network.
}
\examples{
data(ORN)
rtNEL<-readshpnw(ORN.nt, ELComputed=TRUE)
res.sl<-SL.extraction(rtNEL[[2]],rtNEL[[3]])
res.me<-ME.simplification(res.sl[[1]],res.sl[[2]],DegreeL=res.sl[[4]]) 
res.pn<-PN.amalgamation(res.me[[1]],res.me[[2]],DegreeL=res.me[[4]])
ptcoords<-Nodes.coordinates(res.pn[[1]])
#plot(ORN.nt)
#points(ptcoords, col="green")
#plot(ORN.nt)
#points(Nodes.coordinates(rtNEL[[2]]), col="red")
}
\seealso{\code{\link{SL.extraction}}, \code{\link{ME.simplification}}}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{pseudo, node}
