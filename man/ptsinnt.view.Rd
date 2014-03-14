\name{ptsinnt.view}
\alias{ptsinnt.view}
\title{Visualize the result of integrating a set of data points with given network}
\description{
The function is to visualize the result of integrating a set of data points with 
a given network. It is an immediate visualization tool of 
the result from function \link{points2network} to give user an impression 
how the points connected with the given network.
}
\usage{
ptsinnt.view(ntdata, nodelist, pointsxy, CoorespondIDs, VElist=NULL)
}
\arguments{
  \item{ntdata}{A \dQuote{SpatialLines} or \dQuote{SpatialLinesDataFrame} object;}
  \item{nodelist}{An nodelist object, see \link{nodelist};}
  \item{pointsxy}{A two columns vector of X-Y coordinates of the given points 
                  set;}
  \item{CoorespondIDs}{A vector(of the length equal to the number of points) of 
                       cooresponding node IDs of each point;}
  \item{VElist}{A list of virtual edges if existed}
}
\examples{
data(ORN)
pts<-spsample(rn, 100, type="random")
ptsxy<-coordinates(pts)[,1:2]
ptsxy<-cbind(ptsxy[,1]+0.008,ptsxy[,2]+0.008)
#Mapping each point to the nearest node in the network/graph
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=1)
#Visualize the results without virtual edges 
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy, 
             CoorespondIDs=res[[3]])
#Visualize the results with virtual edges
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=3,
                    ea.prop=rep(0,37))
VElist<-res[[7]]
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy, 
             CoorespondIDs=res[[3]], VElist=VElist)
}
\seealso{
\code{\link{points2network}}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
