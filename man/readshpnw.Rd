\name{readshpnw}
\alias{readshpnw}
\title{Read a network from a \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object}
\description{
This function is to split nodes and edges from a \dQuote{SpatialLines} or 
\dQuote{SpatialLinesDataFrame} object, and return \dQuote{nodelist} and \dQuote{edgelist}.
}
\usage{
readshpnw(data=list(), ELComputed=FALSE, longlat=FALSE, Detailed=FALSE, ea.prop=NULL)
}
\arguments{
  \item{data}{A \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object;}
  \item{ELComputed}{if TRUE, the edge length will be computed}
  \item{longlat}{If TRUE, distance will be calculated on an ellipse with WGS84 
                 parameters;}
  \item{Detailed}{Default is FALSE, only two endpoints (starting and ending) of 
                  each polyline are recognised as nodes and collected in the 
                  \dQuote{nodelist}; If TRUE, all the endpoints will be recognised
                  as nodes and collected.}
  \item{ea.prop}{If data is a \dQuote{SpatialLinesDataFrame} object and \dQuote{Detailed} 
  is TRUE, the ea.prop shoul be a given as a vector (of lenght equal to the number 
  of columns in data.frame(data)) with values 0 or 1, for defining the rules of 
  re-attributing the new edges: 0 means that equalization is used for the arrtibute
  inherition from the original data;  1 implies that weighted-mean based on the 
  edge length is adopted.}  
}
\value{
A list consisted of:
  \item{Detailed}{TRUE if the output is under a \dQuote{Detailed} mode, and 
                  \dQuote{edgelist} will have a different structure;}
  \item{nodelist}{A \dQuote{nodelist} object}
  \item{edgelist}{An \dQuote{edgelist} object}
  \item{Edgelength}{A vector (of length equal to the number of edges) of edge 
                    lengths if \dQuote{ELComputed} is TRUE;}
  \item{Eadf}{A data frame of edge attributes, [EdgeID,...(items extracted from 
             the \dQuote{SpatialLinesDataFrame} object )...]}
  \item{nodexlist}{A vector contains X-coordinates of all the nodes}
  \item{nodeylist}{A vector contains Y-coordinates of all the nodes}  
}
\details{
This function is the first step to convert a newtwork data into an object of 
\dQuote{igraph}. With a given \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame},
it produces a {nodelist}, \dQuote{edgelist}, and data frame for edge attributes. 
}
\examples{
data(ORN)
rtNEL<-readshpnw(rn)
nl<-rtNEL[[2]]
el<-rtNEL[[3]]
#Compute edge length
rtNEL<-readshpnw(rn, ELComputed=TRUE)
edgelength<-rtNEL[[4]]
eadf<-rtNEL[[5]]
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{spatial, graph}
