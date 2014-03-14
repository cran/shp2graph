\name{points2network}
\alias{points2network}
\alias{Nearest.nodes}
\alias{footpoint.nodes}
\alias{Add.nodes}
\alias{virtualedge.nn}
\alias{point.in.bbox}
\title{Integrate a point data set with a network}
\description{
In order to associate a  point data set with a given network, this function finds a 
corresponding node for each point in the data set under a certain rule(see details).
}
\usage{
points2network(ntdata,pointsxy,mapping.method=1,ELComputed=FALSE,longlat=F,
               Detailed=F, ea.prop=NULL)
}
\arguments{                                                                   
  \item{ntdata}{a \dQuote{SpatialLinesDataFrame} object}
  \item{pointsxy}{A two column matrix contains X-Y coordinates of the points.}
  \item{mapping.method}{Mapping methods between the given points 
                        and network will be used, see the details below}
  \item{ELComputed}{if TRUE, the edge length will be computed}
  \item{longlat}{if TRUE, use distances on an ellipse with WGS84 parameters}
  \item{Detailed}{if TRUE, all the points of polylines will be regarded as nodes
                 in the network; if FALSE, only endpoints of polylines are treated 
                 as nodes}
  \item{ea.prop}{A vector contains 1 or 0 values (the length equals to the number 
                 of arributes of the \dQuote{ntdata}):
                 
                 1: the corresponding attribute will be kept and re-orgnised for 
                 the new network;
                 
                 0: the corresponding attribute will not be kept.}  
}
\value{
A list consisted of:
  \item{nodelist}{A \dQuote{nodelist} object}
  \item{edgelist}{An \dQuote{edgelist} object}
  \item{CoorespondIDs}{A vector contains cooresponding node IDs for each point}
  \item{nodexlist}{A vector contains X-coordinates of all the nodes}
  \item{nodeylist}{A vector contains Y-coordinates of all the nodes} 
  \item{Eadf}{a data frame of edge attributes, dataframe(EdgeID,...
             (attributes extracted from original file)...)}
  \item{VElist}{A list of virtual edges if existed}
  \item{Edgelength}{If ELComputed is TRUE, Edgelength will be a vector consisted
   of edge lengths cooresponding to each edge, else it will be NULL}
   
}
\details{
The value of mapping.method can be 1,2,3 or 4, which respectively means:

1: Mapping each point to the nearest node in the network/graph

2: Mapping each point to the nearest point (add them as nodes if they are not)
   on the network
   
3: Add a new edge(virtual edge) between each point and the nearest node

4: Add a new edge(virtual edge) between each point and the nearest point
}
\seealso{\code{\link{ptsinnt.view}}}
\examples{
\dontrun{
data(ORN)
pts<-spsample(rn, 100, type="random")
ptsxy<-coordinates(pts)[,1:2]
ptsxy<-cbind(ptsxy[,1]+0.008,ptsxy[,2]+0.008)
#Mapping each point to the nearest node in the network/graph
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=1)
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy,
             CoorespondIDs=res[[3]])
#Mapping each point to the nearest point (add them as nodes if they are not) on 
#the network
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=2,ea.prop=rep(0,37))
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]]) 
#Add a new edge(Virtual edge) between each point and the nearest node
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=3,ea.prop=rep(0,37))
VElist<-res[[7]]
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]], 
             VElist=VElist)
#Add a new edge(Virtual edge) between each point and the nearest point
res<-points2network(ntdata=rn,pointsxy=ptsxy, mapping.method=4,ea.prop=rep(0,37))
VElist<-res[[7]]
ptsinnt.view(ntdata=rn, nodelist=res[[1]], pointsxy=ptsxy, CoorespondIDs=res[[3]],
             VElist=VElist) 
}
}
\author{Binbin Lu \email{lubinbin220@gmail.com}}
\keyword{spatial, graph}
