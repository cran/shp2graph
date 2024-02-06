\name{Redef.functions}
\alias{Redef.functions}
\title{A collection of functions for redefining attributes of edited edges}
\description{
A function to provide a collection of functions to redefine attributes of edited
 edges in the functions \link{PN.amalgamation} and \link{ME.simplification}.
}
\usage{
Redef.functions(v, typ=1)
}
\arguments{
 \item{v}{An input vector for the specified function;}
 \item{typ}{ A value equalling to 1, 2, 3 or 4 to specify a redefinition function:
  
  1->sum(v)
  
  2->min(v)
  
  3->max(v)
  
  4->mean(v)}
}
\note{
It is easy to include more specific methods for attribute redefinition via 
incorporating the corresponding functions. In other words, you can define more 
methods by numbering them with 5, 6, etc.  
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{attribute redefinition}
