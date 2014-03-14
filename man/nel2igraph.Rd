\name{nel2igraph}
\alias{nel2igraph}
\title{Produce a \dQuote{igraph} object}
\description{
This function is to produce an \dQuote{igraph} object using the extracted 
nodelist and edgelist from a spatial network.
}
\usage{
nel2igraph(nodelist, edgelist, weight = NULL, eadf = NULL, mode = "all")
}
\arguments{
  \item{nodelist}{A \dQuote{nodelist} object}
  \item{edgelist}{An \dQuote{edgelist} object}
  \item{weight}{A numberic vector for weighting all the edges in the edgelist, 
                of which the length equals to the number of edges;}
  \item{eadf}{Attribute data frame for all the edges, of which the number of 
              rows equals to the number of edges;}
  \item{mode}{Character scalar, it specifies whether the graph to create is 
  undirected (\dQuote{all} or \dQuote{total}) or directed; and in the latter case, whether it 
  contains the outgoing (\dQuote{out}) or the incoming (\dQuote{in}) neighbors of the vertices.}  
}
\value{
  \item{gr}{An object of \dQuote{igraph}} 
}
\details{
1. The vector \dQuote{weight} will be the default weights of the edges for any 
relative computation carried out on the produced \dQuote{igraph} object.

2. The coordinates of vertices are attached as attributes \dQuote{X} and \dQuote{Y}, 
which could be retrived with function \dQuote{get.vertex.attribute} from 
package \pkg{igraph}.
}
\examples{
data(ORN)
rtNEL<-readshpnw(rn, ELComputed=TRUE)
#Add the edge length as the weight for graph
igr<-nel2igraph(rtNEL[[2]],rtNEL[[3]],weight=rtNEL[[4]])
coord<-cbind(get.vertex.attribute(igr, "X"), get.vertex.attribute(igr, "Y"))
plot.igraph(igr)
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{igraph}