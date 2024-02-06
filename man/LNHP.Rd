\name{LNHP}
\alias{LNHP}
\alias{LN.prop}
\docType{data}
\title{London house price and hedonic data (SpatialPointsDataFrame)}
\description{
  This data set is sampled from a house price data set provided by the Nationwide 
  Building Society of the UK, and combined with various hedonic variables.
}
\usage{data(LNHP)}
\format{
 A \dQuote{SpatialPointsDataFrame} object named \dQuote{LN.prop}.
  }
\details{
  In this data set, attributes are interpreted as follows:
  
  PURCHASE: sale price in GBP
  
  FLOORSZ: floor size of the property in square metres;
  
  BATH2: 1 if the property has 2 or more bathrooms, 0 otherwise; 
  
  BEDS2: 1 if the property has 2 or more bedrooms, 0 otherwise;
	
  CENTHEAT: 1 if the property has central heating, 0 otherwise;
  
  GARAGE1: 1 if the property has one or more garages, 0 otherwise;
  
  BLDPWW1: 1 if the property was built prior to 1914, 0 otherwise;
  
  BLDINTW: is 1 if the property was built between 1914 and 1939, 0 otherwise;
  
  BLD60S: 1 if the property was built between 1960 and 1969, 0 otherwise;
  
  BLD70S: 1 if the property was built between 1970 and 1979, 0 otherwise;
  
  BLD80S: 1 if the property was built between 1980 and 1989, 0 otherwise;
  
  BLD90S: 1 if the property was built between 1990 and 2000, 0 otherwise;
  
  TYPEDETCH: 1 if the property is detached (i.e. it is stand-alone), 0 otherwise; 
  
  TYPETRRD: 1 if the property is in a terrace of similar houses, 0 otherwise;
  
  TYPEFLAT: 1 if the property is a flat or apartment, 0 otherwise;
  
  PROF: percentage of the workforce in professional or managerial occupations in 
        the census enumeration district in which the house is located.
  
  UNEMPLOY: percentage of unemployment in  the census enumeration district
            in which the house is located. 
  
  X: X coordinate.
  
  Y: Y coordinat.
 }
\references{
  Lu, B., M. Charlton, P. Harris & A. S. Fotheringham (2014) Geographically weighted 
  regression with a non-Euclidean distance metric: a case study using hedonic house 
  price data. International Journal of Geographical Information Science, 28, 660-681.
}
