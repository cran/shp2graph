\name{Degree.list}
\alias{Degree.list}
\title{Degree (In-degree and Out-degree)of nodes}
\description{
The function returns degrees of nodes from provided nodelist and edgelist; 
in-degrees and out-degrees are returned if edges are directed).
}
\usage{
Degree.list(nodelist, edgelist, Directed=F)
}
\arguments{
 \item{nodelist}{A \dQuote{nodelist} object}
 \item{edgelist}{An \dQuote{edgelist} object}
 \item{Directed}{TRUE value if \dQuote{edgelist} are directed; otherwise, FALSE} 
}
\value{
Lists of different contents are returned for undirected and directed edges 
respectively:
For undirected type:
  \item{DegreeL}{An integer vector of degrees cooresponding to  each node in 
                \dQuote{nodelist}}
For directed type:
  \item{InDegreeL}{An integer vector of In-degrees cooresponding to  each node 
                  in \dQuote{nodelist}}
  \item{OutDegreeL}{An integer vector of Out-degrees cooresponding to  each node 
                    in \dQuote{nodelist}}  
}
\note{
The outputs of this function are different between undirected and directed 
networks, actually DegreeL can also be computed by DegreeL=InDegreeL+OutDegreeL.
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{degree, node}
