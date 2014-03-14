\name{PN.amalgamation}
\alias{PN.amalgamation}
\alias{Minus.DegreeL}
\alias{Update.nodelist}
\alias{Update.edgelist}
\alias{extend.eadf}
\alias{is.exception} 
\title{Amalgamate edges connected by a pseudo-node}
\description{
This function is to amalgamate edges connected by a pseudo-node, i.e. to get rid 
of pseudo-nodes in a network.
}
\usage{
PN.amalgamation(nodelist, edgelist, eadf=NULL, ea.prop=NULL, Directed=F, 
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
For undirected type:
  \item{newNodelist}{The new nodelist with pseudo-nodes removed;}
  \item{newEdgelist}{The new edgelist with edges connected by a pseudo-node 
                     amalgamated;}
  \item{newEadf}{The new attribute data frame for the returned 
                    \dQuote{edgelist};}
  \item{DegreeL}{The new degree vector cooresponding to the returned 
                   \dQuote{nodelist};}
For directed type:
  \item{newNodelist}{The new \dQuote{nodelist} with pseudo-node removed;}
  \item{newEdgelist}{The new \dQuote{edgelist} with edges connected by a pseudo-node 
                     amalgamated;}
  \item{newEadf}{The new attribute data frame for the returned \dQuote{edgelist};}
  \item{InDegreeL}{The new in-degree vector cooresponding to the returned 
                  \dQuote{nodelist};}
  \item{OutDegreeL}{The new out-degree vector cooresponding to the returned 
                   \dQuote{nodelist};}  
}
\note{
If input edges are directed, a node is recognised as a pseudo-node with in-degree 
and out-degree equalling to 1, then InDegreeL and OutDegreeL are required; 
Otherwise,  a node is recognised as a pseudo-node with degree equalling to 2.
}
\examples{
data(ORN)
rtNEL<-readshpnw(rn, ELComputed=TRUE)
res.sl<-SL.extraction(rtNEL[[2]],rtNEL[[3]])
res.me<-ME.simplification(res.sl[[1]],res.sl[[2]],DegreeL=res.sl[[4]]) 
res.pn<-PN.amalgamation(res.me[[1]],res.me[[2]],DegreeL=res.me[[4]])
ptcoords<-Nodes.coordinates(res.pn[[1]])
plot(rn)
points(ptcoords, col="green")
plot(rn)
points(Nodes.coordinates(rtNEL[[2]]), col="red")
}
\seealso{\code{\link{SL.extraction}}, \code{\link{ME.simplification}}}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{pseudo, node}
