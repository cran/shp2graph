\name{LNNT}
\alias{LNNT}
\alias{LN.nt}
\docType{data}
\title{London Road network data (SpatialLinesDataFrame)}
\description{
  The London Road network data set is produced by the UK Ordnance Survey (OS) in 2001. 
}
\usage{data(LNNT)}
\format{
 A \dQuote{SpatialLinesDataFrame} object named \dQuote{LN.nt}.
  }
 \details{
  In this data set, roads are categorized into four types (nt_RoadTyp):
  
  Mo: motorway
  
  a: A-road
  
  b: B-road
  
  Mi: minor road
  
  Speed limit for each road segment is also incoporated in the attribute \dQuote{Speed}
  (Transport for London, 2011)
 }
\references{
  Lu, B., M. Charlton, P. Harris & A. S. Fotheringham (2014) Geographically weighted 
  regression with a non-Euclidean distance metric: a case study using hedonic house 
  price data. International Journal of Geographical Information Science, 28, 660-681.
  
  Transport for London, (2011) Intelligent Speed Adaptation [online]. 
  http://www.tfl.gov.uk/corporate/projectsandschemes/7893.aspx#For_developers. 
  Available from: http://www.tfl.gov.uk/corporate/projectsandschemes/7893.aspx#For_developers.
}
\keyword{Road network data,Speed,travel time, example}
