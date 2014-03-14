\name{nodelist}
\alias{nodelist}
\title{An object defined for extracted nodes from a network}
\description{
This is a structure for nodes extracted from the spatial network, and each node 
cooresponds to one row in the object; for each row, the structure is defined 
as [NodeID,coordinate(X,Y)].
}
\details{
This list is the immediate result of extracting nodes from a \dQuote{SpatialLines} 
or \dQuote{SpatialLinesDataFrame} object, see \link{readshpnw}. As a transition 
object, it is a necessary input for many functions in this package, 
like \link{Degree.list}, \link{ME.simplification}, \link{Nodes.coordinates}. 
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{node}
