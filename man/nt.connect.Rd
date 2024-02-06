\name{nt.connect}
\alias{nt.connect}
\title{Check the connectivity of a given network}
\description{
A function to check the connectivity of a given network.
}
\arguments{
  \item{nt}{A \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object.}  
}
\usage{
nt.connect(nt)
}
\details{
In this function, all the nodes are traversed, and classified dynamically via a 
rule, that any pair of different nodes fall into a same category if they could 
be reached from each other. A map will be plotted spontaneously with differet categories 
(i.e. self-connected parts) in distinct colors. Finally, the self-connected part with
the most nodes will be returned as a new \dQuote{SpatialLinesDataFrame} object. 
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{network topology}
