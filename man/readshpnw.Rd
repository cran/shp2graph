\name{readshpnw}
\alias{readshpnw}
\title{Read a network from a \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object}
\description{
A function to read a spatial network from a \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object,
and resolve it into a  \dQuote{\link{nodelist}} and \dQuote{\link{edgelist}} for the following conversion.
}
\usage{
readshpnw(ntdata, ELComputed=FALSE, longlat=FALSE, Detailed=FALSE, ea.prop=NULL)
}
\arguments{
  \item{ntdata}{A \dQuote{SpatialLinesDataFrame} or \dQuote{SpatialLines} object}
  \item{ELComputed}{If TRUE, the length of each edge will be calculated and returned}
  \item{longlat}{If TRUE, distances on an ellipse with WGS84 parameters will be returned}
  \item{Detailed}{If TRUE, all the vertices within the polylines of a spatial network
                  will be recognised as nodes; if FALSE, only two endpoints of 
                  each polyline are treated as nodes}
  \item{ea.prop}{A vector (of the length equalling to the number of edge attributes in
                 in the network \dQuote{ntdata}) with values of 1 or 0:
                 
                 1: the corresponding column in the data frame will be kept for 
                 the new network;
                 
                 0: the corresponding column in the data frame will not be kept.} 
}
\value{
A list consisted of:
  \item{Detailed}{TRUE if the output is under a \dQuote{Detailed} mode;}
  \item{nodelist}{A \dQuote{\link{nodelist}} object}
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{Edgelength}{If \dQuote{ELComputed} is TRUE, lengths of all the edges will be 
                    returned as a numeric vector; otherwise it will be NULL.}
  \item{Eadf}{A data frame of attributes of the returned edges, of which the structure is
   a data frame with [EdgeID,(attributes inherited from the input network data)]}
  \item{nodexlist}{A vector contains X-coordinates of all the nodes}
  \item{nodeylist}{A vector contains Y-coordinates of all the nodes}  
}
\details{
This function plays the first step to convert a spatial network (in a \dQuote{SpatialLines} 
or \dQuote{SpatialLinesDataFrame} object) into an \dQuote{igraph} object by returning
a \dQuote{\link{nodelist}} and \dQuote{\link{edgelist}}. 
}
\examples{
data(ORN)
rtNEL<-readshpnw(ORN.nt)
nl<-rtNEL[[2]]
el<-rtNEL[[3]]
#Compute edge length
rtNEL<-readshpnw(ORN.nt, ELComputed=TRUE)
edgelength<-rtNEL[[4]]
eadf<-rtNEL[[5]]
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{spatial, graph}
