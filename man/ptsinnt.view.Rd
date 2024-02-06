\name{ptsinnt.view}
\alias{ptsinnt.view}
\title{Visualize the integration between a point data set and spatila network}
\description{
A specific function to visualizethe integration between a point data set and spatila 
network, i.e. results returned by the function \code{\link{points2network}}
}
\usage{
ptsinnt.view(ntdata, nodelist, pointsxy, CoorespondIDs, VElist=NULL)
}
\arguments{
  \item{ntdata}{A \dQuote{SpatialLinesDataFrame} or \dQuote{SpatialLines} object}
  \item{nodelist}{A \dQuote{\link{nodelist}} object}
  \item{pointsxy}{A two-column matrix of point coordinates (X, Y)}
  \item{CoorespondIDs}{A vector of the cooresponding node ID for each point in \dQuote{pointsxy}}
  \item{VElist}{A list of virtual edges if added, i.e. approaches 3 and 4 are used}
}
\examples{
data(ORN)
#pts<-spsample(ORN.nt, 100, type="random")
#ptsxy<-coordinates(pts)[,1:2]
#ptsxy<-cbind(ptsxy[,1]+0.008,ptsxy[,2]+0.008)
#Mapping each point to the nearest node in the network/graph
#res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, mapping.method=1)
#Visualize the results without virtual edges 
#ptsinnt.view(ntdata=ORN.nt, nodelist=res[[1]], pointsxy=ptsxy, 
#             CoorespondIDs=res[[3]])
#Visualize the results with virtual edges
#res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, mapping.method=3,
#                    ea.prop=rep(0,37))
#VElist<-res[[7]]
#ptsinnt.view(ntdata=ORN.nt, nodelist=res[[1]], pointsxy=ptsxy, 
#             CoorespondIDs=res[[3]], VElist=VElist)
}
\seealso{
\code{\link{points2network}}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{visulization}
