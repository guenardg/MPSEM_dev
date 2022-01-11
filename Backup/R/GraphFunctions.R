#
### MPSEM package R functions and C wrappers
#

#

#
add.edge <- function(x,from,to,edge=list(),label=NULL) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(length(from) != length(to))
    stop("Number of origins(",length(from),") mismatch that of destinations (",length(to),").")
  if(!is.list(edge))
    stop("Values for edges must be provided as a list.")
  if(length(edge))
    for(i in 1L:length(edge))
      if(length(edge[[i]]) != length(from))
        stop("Edge property '",names(edge)[i],"' has length ",length(edge[[i]])," but the ",length(from)," edges are to be created.")
  if(!is.null(label)) {
    if(is.character(label)) {
      if(length(label) != length(from))
        stop(length(label)," labels are provided, but ",length(from)," are required.")
    } else {
      stop("Labels should be of type character.")
    }
  } else {
    label = as.character(attr(x,"ev")[1L]+(1L:length(from)))
  }
  if(is.character(from)) {
    safe <- from
    from <- match(from,attr(x,"vlabel"))
    if(any(is.na(from)))
      stop("Unknown origin vertices (",paste(safe[which(is.na(from))],collapse=","),").")
  } else {
    if(any(from > attr(x,"ev")[2L]))
      stop("Unknown origin vertices (",paste(from[from > attr(x,"ev")[2L]],collapse=","),").")
  }
  if(is.character(to)) {
    safe <- to
    to <- match(to,attr(x,"vlabel"))
    if(any(is.na(to)))
      stop("Unknown destination vertices (",paste(safe[which(is.na(to))],collapse=","),").")
  } else {
    if(any(to > attr(x,"ev")[2L]))
      stop("Unknown destination vertices (",paste(to[to > attr(x,"ev")[2L]],collapse=","),").")
  }
  x$edge[[1L]] <- c(x$edge[[1L]],from)
  x$edge[[2L]] <- c(x$edge[[2L]],to)
  for (i in names(x$edge)) {
    if(i != "") {
      x$edge[[i]] <- if(is.null(edge[[i]])) c(x$edge[[i]],rep(NA,length(from))) else c(x$edge[[i]],edge[[i]])
    }
  }
  for (i in names(edge)) {
    x$edge[[i]] <- if(is.null(x$edge[[i]])) c(rep(NA,attr(x,"ev")[1]),edge[[i]]) else x$edge[[i]]
  }
  attr(x,"ev")[1L] <- attr(x,"ev")[1L]+length(from)
  attr(x,"elabel") <- c(attr(x,"elabel"),label)
  return(x)
}
#
rm.edge <- function(x,id) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"elabel"))
    if(any(is.na(id)))
      stop("Unknown edge(s) (",paste(safe[which(is.na(id))],collapse=","),").")
  } else {
    if(any(id > attr(x,"ev")[1L]))
      stop("Unknown edge(s) (",paste(id[id > attr(x,"ev")[1L]],collapse=","),").")
  }
  for (i in 1L:length(x$edge))
    x$edge[[i]] <- x$edge[[i]][-id]
  attr(x,"ev")[1L] <- attr(x,"ev")[1L]-length(id)
  attr(x,"elabel") <- attr(x,"elabel")[-id]
  return(x)
}
#
rm.vertex <- function(x,id) {
  if(class(x) != "graph")
    stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"vlabel"))
    if(any(is.na(id)))
      stop("Unknown vertex(es) (",paste(safe[which(is.na(id))],collapse=","),").")
  } else {
    if(any(id > attr(x,"ev")[2L]))
      stop("Unknown vertex(es) (",paste(id[id > attr(x,"ev")[2L]],collapse=","),").")
  }
  x <- rm.edge(x,id=which(!is.na(match(x$edge[[1L]],id)) | !is.na(match(x$edge[[2L]],id))))
  mask <- rep(NA,attr(x,"ev")[2L])
  mask[-id] <- 1L:(attr(x,"ev")[2L]-length(id))
  x$edge[[1L]] <- mask[x$edge[[1L]]] ; x$edge[[2L]] <- mask[x$edge[[2L]]]
  for (i in names(x$vertex))
    x$vertex[[i]] <- x$vertex[[i]][-id]
  attr(x,"ev")[2L] <- attr(x,"ev")[2L]-length(id)
  attr(x,"vlabel") <- attr(x,"vlabel")[-id]
  return(x)
}
#
collapse.vertex <- function(x,id) {
  if(class(x) != "graph")
  stop("Parameter x must be of class graph.")
  if(is.character(id)) {
    safe <- id
    id <- match(id,attr(x,"vlabel"))
    if(any(is.na(id)))
      stop("Unknown vertex(es) (",paste(safe[which(is.na(id))],collapse=","),").")
  } else {
    if(any(id > attr(x,"ev")[2L]))
      stop("Unknown vertex(es) (",paste(id[id > attr(x,"ev")[2L]],collapse=","),").")
  }
  for(i in id) {
    up <- which(!is.na(match(x$edge[[2L]],i))) ; lup <- length(up)
    down <- which(!is.na(match(x$edge[[1L]],i))) ; ldown <- length(down)
    if(!(lup&ldown))    # If the vertex is not an intermediary between other vertex, simply remove it with its edges.
      x <- rm.vertex(x,i)
    else {
      from <- x$edge[[1L]][up] ; to <- x$edge[[2L]][down]
      from <- rep(from,each=ldown) ; to <- rep(to,lup)
      # Prevents the edge that already exist to be recreated.
      lstrip <- lup*ldown ; strip <- rep(FALSE,lstrip)
      for (j in 1L:lstrip)
        strip[j] <- any((from[j] == x$edge[[1L]]) & (to[j] == x$edge[[2L]]))
      if(all(strip))    # If all the intermediary connections already exist, simply remove the vertex with its edges.
        x <- rm.vertex(x,i)
      else {
        if(!is.null(attr(x,"elabel"))) {
          if(!is.null(attr(x,"vlabel")))
            newlab <- paste(attr(x,"vlabel")[from[!strip]],attr(x,"vlabel")[to[!strip]],sep="->")
          else
            newlab <- paste("V#",from[!strip],"->V#",to[!strip],sep="")
        } else
          newlab <- NULL
        if(!is.null(x$edge$length)) {
          ll <- list(length=rep(x$edge$length[up],each=ldown)[!strip] + rep(x$edge$length[down],each=lup)[!strip])
        } else
          ll <- list()
        x <- add.edge(x,from[!strip],to[!strip],ll,newlab)
        x <- rm.vertex(x,i)
      }
    }
  }
  return(x)
}
#
Phylo2DirectedGraph <- function(tp) {
  if(!is.rooted(tp))
    warning("The tree is not rooted. Direction taken from the first edge in the list.")
  if(is.null(tp$node.label))
    tp$node.label <- paste("n",1:tp$Nnode,sep="")
  x <- pop.graph(n=tp$Nnode+length(tp$tip.label),label=c(tp$tip.label,tp$node.label),
                 vertex=list(species=c(rep(TRUE,length(tp$tip.label)),rep(FALSE,tp$Nnode))))
  x <- add.edge(x,from=tp$edge[,1L],to=tp$edge[,2L],
                label=c(paste("E",1L:nrow(tp$edge),sep="")),edge=list(distance=tp$edge.length))
  if(!is.null(tp$root.edge))
    warning("The root edge of the tree has been omitted from the phylogenetic graph.")
  return(x)
}
#
