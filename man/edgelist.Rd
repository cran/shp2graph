\name{edgelist}
\alias{edgelist}
\title{A structure for edge information of a network}
\description{
This is an intergradation of edges from spatial data to graph data, and each 
edge cooresponds to one row in the object, of which the row structure is
designed as [EdgeID,NodeID(from),NodeID(to)].
}
\details{
This object is the immediate result of establishing edges between nodes from a 
\dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object, see 
\link{readshpnw}. As a transition object, it is a necessary input for many 
functions in this package. 
}
\note{
If the parameter \dQuote{Detailed} specified in \link{readshpnw} is TRUE, all 
the endpoints of polyines will be extracted as nodes, then the converted 
graph will have the same spatial details with the transformed \dQuote{SpatialLines} 
or \dQuote{SpatialLinesDataFrame} object. To retrieve the original attributes in 
the \dQuote{SpatialLinesDataFrame} object, the original edge ID is also kept and 
the row structure will be [EdgeID,eid,NodeID(from),NodeID(to)], in which EdgeID 
refers to the new edge id while eid represents the original edge ID.  
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{edge}
