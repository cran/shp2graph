\name{nodelist}
\alias{nodelist}
\title{A designed structure to denote nodes from a spatial network}
\description{
\dQuote{nodelist} is an interchange structure of nodes from a spatial network into 
an \dQuote{igraph} object. It is a two-column data frame, of which each row is 
designed as [NodeID,coordinate(X,Y)].
}
\details{
Both \dQuote{\link{nodelist}} and \dQuote{\link{edgelist}} are interchange structures
of nodes and edges extracted from a spatial network in a \dQuote{SpatialLines} or 
\dQuote{SpatialLinesDataFrame} object. They are always concerned together, and returned
by the function \link{readshpnw} as initial results for the following steps.
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{node list}
