\name{points2network}
\alias{points2network}
\alias{Nearest.nodes}
\alias{footpoint.nodes}
\alias{Add.nodes}
\alias{virtualedge.nn}
\alias{point.in.bbox}
\title{Integrate a point data set into a network}
\description{
A function to integrate an individual point data set into a network under different
rules (see details below).
}
\usage{
points2network(ntdata,pointsxy,approach=1,ELComputed=FALSE,longlat=F,
               Detailed=F, ea.prop=NULL)
}
\arguments{                                                                   
  \item{ntdata}{A \dQuote{SpatialLinesDataFrame} or \dQuote{SpatialLines} object}
  \item{pointsxy}{A two-column matrix of point coordinates (X, Y)}
  \item{approach}{specified by an integer ranging from 1 to 4 to define the approach
  for integration:
  
1: Mapping each point to the nearest node in the network/graph

2: Mapping each point to the nearest point (add them as nodes if they are not)
   on the network
   
3: Add a new edge(virtual edge) between each point and the nearest node

4: Add a new edge(virtual edge) between each point and the nearest point 
}
  \item{ELComputed}{If TRUE, the length of each edge will be calculated and returned}
  \item{longlat}{If TRUE, distances on an ellipse with WGS84 parameters will be returned}
  \item{Detailed}{If TRUE, all the vertices within the polylines of a spatial network
   will be recognised as nodes; if FALSE, only two endpoints of each polyline are treated 
                 as nodes}
  \item{ea.prop}{A vector (of the length equalling to the number of edge attributes in
  in the network \dQuote{ntdata}) with values of 1 or 0:
                 
                 1: the corresponding column in the data frame will be kept for 
                 the new network;
                 
                 0: the corresponding column in the data frame will not be kept.}  
}
\value{
A list consisted of:
  \item{nodelist}{A \dQuote{\link{nodelist}} object}
  \item{edgelist}{An \dQuote{\link{edgelist}} object}
  \item{CoorespondIDs}{A vector of the cooresponding node ID for each point in \dQuote{pointsxy}}
  \item{nodexlist}{A vector contains X-coordinates of all the nodes}
  \item{nodeylist}{A vector contains Y-coordinates of all the nodes} 
  \item{Eadf}{A data frame of attributes of the returned edges, of which the structure is
   a data frame with [EdgeID,(attributes inherited from the input network data)]}
  \item{VElist}{A list of virtual edges if added, i.e. approaches 3 and 4 are used}
  \item{Edgelength}{If \dQuote{ELComputed} is TRUE, lengths of all the edges will be 
  returned as a numeric vector; otherwise it will be NULL.}
}
\seealso{\code{\link{ptsinnt.view}}}
\examples{
\dontrun{
data(ORN)
pts<-spsample(ORN.nt, 100, type="random")
ptsxy<-coordinates(pts)[,1:2]
ptsxy<-cbind(ptsxy[,1]+0.008,ptsxy[,2]+0.008)
#Mapping each point to the nearest node in the network/graph
res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, approach=1)
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy,
             CoorespondIDs=res[[3]])
#Mapping each point to the nearest point (add them as nodes if they are not) on 
#the network
res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, approach=2,ea.prop=rep(0,37))
ptsinnt.view(ntdata=ORN.nt, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]]) 
#Add a new edge(Virtual edge) between each point and the nearest node
res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, approach=3,ea.prop=rep(0,37))
VElist<-res[[7]]
ptsinnt.view(ntdata=ORN.nt, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]], 
             VElist=VElist)
#Add a new edge(Virtual edge) between each point and the nearest point
res<-points2network(ntdata=ORN.nt,pointsxy=ptsxy, approach=4,ea.prop=rep(0,37))
VElist<-res[[7]]
ptsinnt.view(ntdata=ORN.nt, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]],
             VElist=VElist) 
}
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}}
\keyword{point-network projection}
