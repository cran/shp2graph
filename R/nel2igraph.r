nel2igraph<-function(nodelist, edgelist, weight = NULL, eadf = NULL, mode = "all")
{
   nodes <- nodelist[, 1]
   Ne <- length(edgelist[, 1])
   Nn <- length(nodes)
   edL = vector("list", length = Nn)
   for (i in 1:Nn) {
        edL[[i]] <- numeric(0)
    }
    names(edL) <-as.character(nodes)
    if (!is.null(weight)) {
        if (length(weight) != Ne && is.numeric(weight))
            stop("Please give right edge weight, which must be numeric and the same length as edges elment")
    }
    if (!is.null(eadf)) {
        if (length(eadf[, 1]) != Ne)
            stop("The eadf must be numeric and the same length as edges elment")
    }
  # total=Ne
   #pb <- txtProgressBar(min = 0, max = total, style = 3)
   if (mode=="all")
       Directed=F
   else
       Directed=T
   if (Directed)
   {
       for (i in 1:Ne) {
            oidx <- which(nodes == edgelist[i, 2])
            eidx <- which(nodes == edgelist[i, 3])
            edL[[oidx]]<- c(edL[[oidx]], eidx-1)
    #    setTxtProgressBar(pb, i)
        }
   }
   else {
       for (i in 1:Ne) {
            oidx <- which(nodes == edgelist[i, 2])
            eidx <- which(nodes == edgelist[i, 3])
            edL[[oidx]]<- c(edL[[oidx]], eidx-1)
            edL[[eidx]]<- c(edL[[eidx]], oidx-1)
    #    setTxtProgressBar(pb, i)
        }
   }
   gr<-graph.adjlist(edL,mode=mode)
   gr <- set.vertex.attribute(gr,"X", index=(1:Nn)-1, Nodes.coordinates(nodelist)[,1])
   gr <- set.vertex.attribute(gr,"Y", index=(1:Nn)-1, Nodes.coordinates(nodelist)[,2])
   #total=Ne
   #pb <- txtProgressBar(min = 0, max = total, style = 3)
   #gr<-set.edge.attribute(gr, "weight", value=0)
   if (!is.null(eadf)){
     eanms<-colnames[eadf]
     for (i in 1:length(eanms)) gr<-set.edge.attribute(gr, name=eanms[i], value=0)
   }
   gr.es<-E(gr)
   for (i in 1:Ne)
   {
     oidx <- which(nodes == edgelist[i, 2])-1
     eidx <- which(nodes == edgelist[i, 3])-1
     if (Directed)
     {
        gr<-set.edge.attribute(gr, "weight", index= gr.es[oidx%->%eidx], value=weight[i])
       #if (!is.null(eadf)) {
           #eanms<-colnames[eadf]
     }
     else
     {
       gr<-set.edge.attribute(gr, "weight", index= gr.es[oidx%--%eidx], value=weight[i])
     }
   }
    gr
}