## **************************************************************************
##
##    (c) 2010-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Directed Graph - Utility Functions **
##
##    This file is part of MPSEM
##
##    MPSEM is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    MPSEM is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with MPSEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Graph Utility Functions
#' 
#' @description A suite of graph utility functions.
#' 
#' @name graph-utils
#' 
#' @aliases getOrigin getConnected getTerminal reorderGraph graphDist isTree
#' isDivergent isLinear
#' 
#' @param x A \code{\link{graph-class}} object.
#' @param order An integer vector of the vertex indices.
#' @param shuffleEdge A Boolean. Whether to randomly shuffle the order that the
#' edges are stored in the \code{\link{graph-class}} object (\code{FALSE}).
#' 
#' @details A origin vertex is one having only outgoing edge(s) and no incoming
#' edge, whereas a terminal vertex is one having only incoming edge(s) and no
#' outgoing edge. A non-connected vertex has no edge, whereas a connected vertex
#' may have incoming edge(s), outgoing edge(s), or both.
#' 
#' Reordering a graph with a \code{processOrder} attribute will come with a
#' recalculation of the process order, whereas doing so on a graph with a
#' \code{dist} attribute cause the pairwise distance matrix to also be
#' reordered.
#' 
#' @return
#' \describe{
#'   \item{getOrigin}{A vector of integer.}
#'   \item{getConnected}{A vector of integer.}
#'   \item{getNonConnected}{A vector of integer.}
#'   \item{getTerminal}{A vector of integer.}
#'   \item{reorderGraph}{A \code{\link{graph-class}} object.}
#'   \item{graphDist}{A pairwise distance matrix such as the one obtained from
#'   function \code{\link[stats]{dist}}.}
#'   \item{isTree}{A \code{logical} stipulating whether the graph is a tree.}
#'   \item{isDivergent}{A \code{logical} stipulating whether the graph has
#'   divergence.}
#'   \item{isLinear}{A \code{logical} stipulating whether the graph is a linear
#'   sequence of vertices.}
#' }
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @seealso \code{\link{graph-class}}.
#' 
#' @examples ## Create and example graph with 10 vertices and 16 edges:
#' pop.graph(
#'   n = 10,
#'   vertex = list(
#'     species = rep(TRUE,10),
#'     x = c(2,3,2,4,3,4,2,1,1,0),
#'     y = c(-2,1,2,0,-0.5,-2,0,-1,1,0)
#'   ),
#'   label = sprintf("V%d",1:10)
#' ) %>%
#'   add.edge(
#'     from = c(10,10,9,9,8,8,3,7,7,10,2,2,5,1,4,5),
#'     to = c(9,8,3,7,7,1,2,2,5,2,1,4,4,4,6,6),
#'     edge = list(distance=c(1,1,1,1,1,1,1,1,1,4,2,1,1,3,1,1)),
#'     label = sprintf("E%d",1:16)
#'   ) -> x
#' x
#' 
#' getOrigin(x)         ## The graph has a single origin vertex.
#' 
#' getConnected(x)      ## All the vertices
#' getNonConnected(x)   ## are connected.
#' 
#' getTerminal(x)       ## The graph has a single terminal vertex.
#' 
#' isTree(x)            ## The graph is not a tree.
#' isDivergent(x)       ## The graoh has divergences.
#' isLinear(x)          ## The graph is not a linear vertex sequence.
#' 
#' ## The average pairwise distances between the vertices:
#' graphDist(x)
#' 
#' ## Reordering of the vertices:
#' xr <- reorderGraph(x,c(5:1,8,6,7,10,9))
#' xr
#' 
#' getOrigin(xr)     ## Same origin vertex, but at a different index.
#' getTerminal(xr)   ## Same terminal vertex, but at a different index.
#' graphDist(xr)     ## Same distances, but in a different order.
#' 
NULL
#' 
#' @describeIn graph-utils
#' 
#' Get Origin Vertex
#' 
#' Obtain the origin vert(ex/ices) of a directed graph; an origin vertex is one
#' with no incoming edge.
#' 
#' @export
getOrigin <- function(x) {
  
  nv <- attr(x, "ev")[2L]
  
  tmp <- rep(FALSE, nv)
  tmp[x$edge[[1L]]] <- TRUE
  tmp[x$edge[[2L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- attr(x, "vlabel")[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Connected Vertex
#' 
#' Obtain the connected vert(ex/ices) of a graph.
#' 
#' @export
getConnected <- function(x) {
  
  nv <- attr(x, "ev")[2L]
  
  tmp <- rep(FALSE, nv)
  tmp[x$edge[[1L]]] <- TRUE
  tmp[x$edge[[2L]]] <- TRUE
  
  tmp <- which(tmp)
  names(tmp) <- attr(x, "vlabel")[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Non-connected Vertex
#' 
#' Obtain the non-connected connected vert(ex/ices) of a graph.
#' 
#' @export
getNonConnected <- function(x) {
  
  nv <- attr(x, "ev")[2L]
  
  tmp <- rep(TRUE, nv)
  tmp[x$edge[[1L]]] <- FALSE
  tmp[x$edge[[2L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- attr(x, "vlabel")[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Get Terminal Vertex
#' 
#' Obtain the terminal vert(ex/ices) of a directed graph; a terminal vertex is
#' one with no outgoing edge.
#' 
#' @export
getTerminal <- function(x) {
  
  tmp <- rep(TRUE, attr(x,"ev")[2L])
  tmp[x$edge[[1L]]] <- FALSE
  
  tmp <- which(tmp)
  names(tmp) <- attr(x,"vlabel")[tmp]
  
  tmp
}
#' 
#' @describeIn graph-utils
#' 
#' Reorder Vertices
#' 
#' Reorder the vertices of a directed graph.
#' 
#' @export
reorderGraph <- function(x, order, shuffleEdge = FALSE) {
  
  ev <- attr(x, "ev")
  
  if(length(order) != ev[2L])
    stop("Argument 'order' is not of the correct size.")
  if(is.character(order))
    order <- match(order,attr(x,"vlabel"))
  if(!all(order %in% 1L:ev[2L]))
    stop("Incorrect vertex reference(s) in 'order'")
  if(!all(1L:ev[2L] %in% order))
    stop("Missing / duplicated vertex reference(s) in 'order'")
  
  pop.graph(
    n = ev[2L],
    vertex = lapply(x$vertex, function(x, o) x[o], o=order),
    label = attr(x, "vlabel")[order]
  ) -> out
  
  if(shuffleEdge) {
    
    so <- sample(ev[1L], replace=FALSE)
    
    add.edge(
      x = out,
      from = match(x$edge[[1L]], order)[so],
      to = match(x$edge[[2L]], order)[so],
      edge = lapply(x$edge[-(1L:2L)], function(x, o) x[o], o=so),
      label = attr(x, "elabel")[so]
    ) -> out
  } else
    add.edge(
      x = out,
      from = match(x$edge[[1L]], order),
      to = match(x$edge[[2L]], order),
      edge = x$edge[-(1L:2L)],
      label = attr(x, "elabel")
    ) -> out
  
  ## Is a process order was defined:
  if(!is.null(attr(x,"processOrder"))) {
    attr(out,"processOrder") <- integer(ev[2L])
    attr(out,"processOrder")[order] <- attr(x,"processOrder")
  }
  
  ## This reordering of the distances has to be verified:
  if(!is.null(attr(x,"dist"))) {
    
    diss <- attr(x,"dist")
    
    for(i in 1L:(ev[2L] - 1L))
      attr(x,"dist")[dst_idx(ev[2L],order[i],order[(i+1L):ev[2L]])] ->
      diss[dst_idx(ev[2L],i,(i+1L):ev[2L])]
    
    attr(diss,"Labels") <- attr(diss,"Labels")[order]
    attr(out,"dist") <- diss
  }
  
  out
}
#' 
#' @describeIn graph-utils
#' 
#' Graph Distance Matrix
#' 
#' Obtain a matrix of the (average) graph distance among the vertices.
#' 
#' @export
graphDist <- function(x) {
  
  if(is.null(x$edge$distance))
    stop("Edges must have a 'distance' property.")
  
  if(length(getOrigin(x)) > 1L)
    stop("This procedure is only suitable for single-origin graphs.")
  
  if(length(getNonConnected(x)))
    stop("All vertices must be connected.")
  
  ev <- attr(x,"ev")
  
  ord <- attr(x,"processOrder")
  
  if(is.null(ord))
    ord <- getProcessOrder(x)
  
  diss <- numeric(ev[2L]*(ev[2L] - 1L)/2L)
  
  for(k in 2L:ev[2L]) {
    
    ## This is the ascendant vertices (as indices in 'ord'):
    asc <- match(x$edge[[1L]][x$edge[[2L]] == ord[k]], ord)
    
    ## transfer the distances associated with the ascendant vertices into the
    ## dissimilarity matrix:
    tstp <- x$edge$distance[x$edge[[2L]] == ord[k]]
    diss[dst_idx(ev[2L], ord[asc], ord[k])] <- tstp
    
    ## This is the remaining vertices:
    others <- (1L:k)[-c(asc,k)]
    
    if(length(others)) {
      
      dd <- diss[dst_idx(ev[2L], ord[asc[1L]], ord[others])] + tstp[1L]
      if(length(asc) > 1L) {
        for(i in 2L:length(asc))
          dd <- dd + (diss[dst_idx(ev[2L], ord[asc[i]], ord[others])] + tstp[i])
        dd <- dd/length(asc)
      }
      diss[dst_idx(ev[2L], ord[k], ord[others])] <- dd
    }
  }
  
  structure(
    diss,
    Size = ev[2L],
    Labels = attr(x,"vlabel"),
    Diag = FALSE,
    Upper = FALSE,
    method = "patristic",
    class = "dist",
    call = match.call()
  )
}
#' 
#' @describeIn graph-utils
#' 
#' Tree Test
#' 
#' Testing whether the graph is a tree.
#' 
#' @export
isTree <- function(x) {
  
  ev <- attr(x,"ev")
  
  for(i in 1L:ev[2L])
    if(sum(x$edge[[2L]] == i) > 1L)
      return(FALSE)
  
  TRUE
}
#' 
#' @describeIn graph-utils
#' 
#' Divergence Test
#' 
#' Testing whether the graph has divergence.
#' 
#' @export
isDivergent <- function(x) {
  
  ev <- attr(x,"ev")
  
  for(i in 1L:ev[2L])
    if(sum(x$edge[[1L]] == i) > 1L)
      return(TRUE)
  
  FALSE
}
#' 
#' @describeIn graph-utils
#' 
#' Linearity Test
#' 
#' Testing whether the graph is a linear sequence.
#' 
#' @export
isLinear <- function(x) isTree(x) && !isDivergent(x)
#'
