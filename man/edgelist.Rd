\name{edgelist}
\alias{edgelist}
\title{A designed structure to denote edges from a spatial network}
\description{
\dQuote{edgelist} is an interchange structure of edges from a spatial network into 
an \dQuote{igraph} object. It is a three-column matrix, of which each row is 
designed as [EdgeID,NodeID(from),NodeID(to)].
}
\details{
Both \dQuote{\link{nodelist}} and \dQuote{\link{edgelist}} are interchange structures
of nodes and edges extracted from a spatial network in a \dQuote{SpatialLines} or 
\dQuote{SpatialLinesDataFrame} object. They are always concerned together, and returned
by the function \link{readshpnw} as initial results for the following steps.
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
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{edge}
