\name{nel2igraph}
\alias{nel2igraph}
\title{Produce an \dQuote{igraph} object}
\description{
A function to produce an \dQuote{igraph} object with the  \dQuote{\link{nodelist}} and 
\dQuote{\link{edgelist}}, which could be returned by the function \link{readshpnw}.
}
\usage{
nel2igraph(nodelist, edgelist, weight = NULL, eadf = NULL, Directed = FALSE)
}
\arguments{
  \item{nodelist}{A \dQuote{\link{nodelist}} object}
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{weight}{A numberic vector to weight all the edges in the \dQuote{\link{edgelist}}, 
                of which the length equals to the number of edges;}
  \item{eadf}{ A data frame of attributes corresponding to all the edges;}
  \item{Directed}{TRUE if edges are directed, FALSE otherwise;}   
}
\value{
  \item{gr}{An \dQuote{igraph} object} 
}
\details{
1. The weighting vector, \dQuote{weight}, will be used as default for any weigted
calculations with edges in the \dQuote{igraph} object.

2. The coordinate of each node is attached as attributes \dQuote{X} and \dQuote{Y}, 
which could be retrived via the function \dQuote{\code{get.vertex.attribute}} from 
the package \pkg{igraph}.
}
\examples{
data(ORN)
rtNEL<-readshpnw(ORN.nt, ELComputed=TRUE)
#Add the edge length as the weight for graph
igr<-nel2igraph(rtNEL[[2]],rtNEL[[3]],weight=rtNEL[[4]])
plot(igr, vertex.label=NA, vertex.size=2,vertex.size2=2)
#plot(ORN.nt)
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{igraph}