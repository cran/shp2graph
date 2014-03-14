\name{Redef.functions}
\alias{Redef.functions}
\title{A collection of functions for redefining attributes of new edges}
\description{
The function is used to redifine attributes of new edges when the optimization 
functions \link{PN.amalgamation} and \link{ME.simplification} are conducted.
}
\usage{
Redef.functions(v, typ=1)
}
\arguments{
 \item{v}{An input vector for the specified function;}
 \item{typ}{ a value from \dQuote{1, 2, 3, 4} to specify the redefinition method:
 1->sum(v), 2->min(v),3->max(v),4->mean(v)}
}
\note{
With specific cases, it is easy to extend this collection with some other functions.
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
