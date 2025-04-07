## **************************************************************************
##
##    (c) 2010-2024 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Class and method for directed graphs **
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
#' Class and Method for Directed Graphs
#' 
#' @description Class and methods to handle MPSEM graphs.
#' 
#' @docType class
#' 
#' @name graph-class
#' 
#' @aliases graph
#' 
#' @param x,object A \code{graph-class} object.
#' @param col Color of the arrows representing the edges. Default is
#' \code{"grey80"}.
#' @param bg Background colors for the origin, intermediate, and terminal
#' vertices, respectively. Default is \code{c("red","black","blue")}.
#' @param pch Plotting 'character' (see \code{\link[graphics]{points}}) for the
#' details. Default is \code{21} (a solid dot with a colored background).
#' @param length Length of the edge of the arrow head (see
#' \code{\link[graphics]{arrows}}). The default value is \code{0.05}.
#' @param pt.cex The relative point size used for the vertex markers. The
#' default value is \code{0.75}.
#' @param value A vector or \code{\link{data.frame}} containing the values to be
#' given to the \code{graph-class} object.
#' @param ... Additional parameters to be passed to other functions of methods.
#' 
#' @details Prints user-relevant information about the graph: number of edges
#' and vertices, edge and vertex labels, addition edge properties and vertex
#' properties.
#' 
#' The plot method verifies whether there are vertex properties 'type', 'x', and
#' 'y'. If any of them is missing, an automatic procedure tries to make a
#' worthwhile plot of the vertices. The result may or may not by suitable to the
#' user. In the latter case, it is possible to edit the graph manually using
#' function \code{\link{graphModplot}}.
#' 
#' @format A \code{graph-class} object contains:
#' \describe{
#'   \item{ edge }{ A list whose first two unnamed members are the indices of
#'   the origin and destination vertices. Additional members must be named and
#'   are additional edge properties (e.g. length). }
#'   \item{ vertex }{ A list that optionally contains vertex properties, if any
#'   (or an empty \code{\link{list}} if none). }
#' }
#' 
#' @return
#' \describe{
#'   \item{print}{\code{NULL}.}
#'   \item{plot}{A \code{\link{graph-class}} object with the vertex properties
#'   'type', 'x', and 'y' added.}
#' }
#' Both methods return invisibly.
#' 
#' @author \packageAuthor{MPSEM}
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution 4: 1120-1131
#' 
#' @seealso \code{\link{graph-functions}}.
#' 
#' @examples ## Create an exemplary graph:
#' pop.graph(
#'   n = 13,
#'   vertex = list(
#'     species = rep(TRUE,13),
#'     type = c(2,2,3,1,2,2,2,2,2,2,3,3,3),
#'     x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5),
#'     y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5)
#'   ),
#'   label = sprintf("V%d",1:13)
#' ) %>%
#'   add.edge(
#'     from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,6,6,9,10,10),
#'     to = c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,3,13,10,12,11),
#'     edge = list(
#'       distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,
#'                    4.8,3.2,3.5,4.4,2.5,3.4,4.3,3.1,2.2)
#'     ),
#'     label = sprintf("E%d",1:19)
#'   ) -> x
#' 
#' ## The print method:
#' x
#' 
#' ## Obtaining some information about the graph:
#' getOrigin(x)          ## It has a single origin called V4.
#' getNonConnected(x)    ## It has no unconnected vertex and thus.
#' getConnected(x)       ## 13 connected vertices: V1-V13.
#' graphDist(x)          ## Its pairwise mean path distance matrix.
#' 
#' ## The plot method uses the vertex coordinates and vertex types, when available
#' plot(x, length=0.1)
#' 
#' ## A copy of 'x':
#' x2 <- x
#' 
#' ## Strip the plotting information from the copy:
#' x2$vertex$type <- NULL
#' x2$vertex$x <- NULL
#' x2$vertex$y <- NULL
#' 
#' ## The plot method now attempts to guess a suitable point disposition:
#' x3 <- plot(x2, length=0.1)
#' 
#' ## the plotting information can be retrieved as follows:
#' x3$vertex$type
#' x3$vertex$x
#' x3$vertex$y
#' 
NULL
#' 
#' @describeIn graph-class
#' 
#' Print Graph
#' 
#' A print method for graph-class objects.
#' 
#' @method print graph
#' 
#' @importFrom utils head tail
#' 
#' @export
print.graph <- function (x, ...) {
  
  ev <- attr(x, "ev")
  
  cat("\nA graph with",ev[2L], "vertices and", ev[1L], "edges\n")
  cat("----------------------------------\n")
  
  vl <- attr(x, "vlabel")
  
  if(!is.null(vl)) {
    
    cat("Vertex labels: ")
    
    if(length(vl) > 10L) {
      cat(paste(c(head(vl,5L), paste("... +", length(vl) - 10L,"more ..."),
                  tail(vl,3L)), collapse=", "))
    } else cat(paste(vl, collapse=", "))
    
    cat("\n")
  }
  
  el <- attr(x, "elabel")
  
  if(!is.null(el)) {
    
    cat("Edge labels: ")
    
    if(length(el) > 10L) {
      cat(paste(c(head(el, 5L), paste("... +", length(el) - 10L,"more ..."),
                  tail(el, 3L)), collapse=", "))
    } else cat(paste(el, collapse=", "))
    
    cat("\n")
  }
  
  if(length(attr(x$vertex, "names") > 0L)) {
    cat("Vertex information: ", paste(attr(x$vertex, "names"), collapse = ", "),
        "\n")
  } else cat("No available vertex information\n")
  
  if(length(attr(x$edge, "names") > 2L)) {
    cat("Available edge information: ",
        paste(attr(x$edge, "names")[-(1L:2L)], collapse = ", "), "\n")
  } else
    cat("No edge information available\n")
  
  cat("\n")
  
  invisible(NULL)
}
#' 
#' @describeIn graph-class
#' 
#' Graph Length
#' 
#' Get the number of vertices in a graph.
#' 
#' @method length graph
#' 
#' @export
length.graph <- function(x) attr(x,"ev")[2L]
#' 
#' @describeIn graph-class
#' 
#' Extract Coordinates
#' 
#' Extracts the display coordinates of a \code{graph-class} object.
#' 
#' @method coordinates graph
#' 
#' @export
coordinates.graph <- function(x) {
  
  if(!is.null(x$vertex$x)) {
    
    out <- data.frame(x=x$vertex$x, row.names = attr(x,"vlabel"))
    
    if(!is.null(x$vertex$y)) {
      
      out$y <- x$vertex$y
      
      if(!is.null(x$vertex$z))
        out$z <- x$vertex$z
    }
    
    out
  } else
    NULL
}
#' 
#' @describeIn graph-class
#' 
#' Set Coordinates
#' 
#' Set the display coordinates of a \code{graph-class} object.
#' 
#' @method coordinates<- graph
#' 
#' @export
`coordinates<-.graph` <- function(x, value) {
  
  if(is.null(value)) {
    
    x$vertex$x <- NULL
    x$vertex$y <- NULL
    x$vertex$z <- NULL
    
  } else {
    
    ev <- attr(x,"ev")
    
    if(NROW(value) != ev[2L]) {
      
      if(NCOL(value) > 1L) {
        
        if(NROW(value) > 1L)
          warning("The number of rows of 'value' is ", NROW(value),
                  ", but the number of vertices is ", ev[2L])
        
        idx <- rep(1L:NROW(value), length.out=ev[2L])
        x$vertex$x <- value[idx,1L]
        x$vertex$y <- value[idx,2L]
        if(NCOL(value) > 2L)
          x$vertex$z <- value[idx,3L]
        
      } else {
        
        value <- unlist(value)
        
        if(length(value) > 1L)
          warning("The length of 'value' is ", length(value),
                  ", but the number of vertices is ", ev[2L])
        
        x$vertex$x <- rep(value, length.out=ev[2L])
      }
    } else {
      if(NCOL(value) > 1L) {
        
        x$vertex$x <- value[,1L]
        x$vertex$y <- value[,2L]
        
        if(NCOL(value) > 2L)
          x$vertex$z <- value[,3L]
        
      } else
        x$vertex$x <- unlist(value, use.names=FALSE)
    }
  }
  
  x
}
#' 
#' @describeIn graph-class
#' 
#' Extract Labels
#' 
#' Extracts the (vertex) labels of a graph.
#' 
#' @method labels graph
#' 
#' @export
labels.graph <- function(object, ...) attr(object,"vlabel")
#' 
#' @describeIn graph-class
#' 
#' Set Labels
#' 
#' Sets the (vertex) labels of a graph.
#' 
#' @method labels<- graph
#' 
#' @export
`labels<-.graph` <- function(object, value) {
  
  if(length(value) != attr(object,"ev")[2L])
    stop("The length of the labels (", length(value), ") does not match the ",
         "number of vertices (", attr(object,"ev")[2L], ")")
  
  attr(object,"vlabel") <- value
  
  object
}
#' 
#' @describeIn graph-class
#' 
#' Extract Names
#' 
#' Extracts the names of the nodal values associated with the vertex of a graph.
#' 
#' @method names graph
#' 
#' @export
names.graph <- function(x) names(x$vertex)
#' 
#' @describeIn graph-class
#' 
#' Set Names
#' 
#' Sets the names of the nodal values associated with the vertex of a graph.
#' 
#' @method names<- graph
#' 
#' @export
`names<-.graph` <- function(x, value) {
  
  if(!is.character(value))
    stop("The 'value' must be of character")
  
  if(length(value) != length(x$vertex))
    stop("The length of the labels (", length(x$vertex),
         ") does not match the number of vertices (", attr(x,"ev")[2L], ")")
  
  names(x$vertex) <- value
  
  x
}
#' 
#' 
#' 
#' @describeIn graph-class
#' 
#' Plot Graph
#' 
#' A plotting method for graph-class objects.
#' 
#' @method plot graph
#' 
#' @importFrom graphics arrows points
#' 
#' @export
plot.graph <- function(x, col = "grey80", bg = c("red","black","blue"),
                       pch = 21L, length = 0.05, pt.cex = 0.75, ...) {
  
  if(is.null(x$vertex$x) || is.null(x$vertex$y) || is.null(x$vertex$type))
    x <- getVertexCoordinate(x)
  
  par <- par(no.readonly=TRUE)
  on.exit(par(par))
  
  par(mar=c(1,1,1,1))
  
  plot(NA, xlim=range(x$vertex$x), ylim=range(x$vertex$y), type="n",
       axes=FALSE, ...)
  
  arrows(
    x0 = x$vertex$x[x$edge[[1L]]],
    x1 = x$vertex$x[x$edge[[2L]]],
    y0 = x$vertex$y[x$edge[[1L]]],
    y1 = x$vertex$y[x$edge[[2L]]],
    col = col,
    length = length,
    ...
  )
  
  points(
    x = x$vertex$x,
    y = x$vertex$y,
    pch = pch,
    bg = bg[x$vertex$type],
    cex = pt.cex,
    ...
  )
  
  invisible(x)
}
#' 
#' @describeIn graph-class
#' 
#' Transformation to a Tree
#' 
#' An \code{\link[ape]{as.phylo}} method to transforms a graph-class object into
#' a "phylo" class object.
#' 
#' @method as.phylo graph
#' 
#' @importFrom ape as.phylo
#' 
#' @export
as.phylo.graph <- function(x, ...) {
  
  ev <- attr(x,"ev")
  
  po <- attr(x,"processOrder")
  if(is.null(po))
    po <- getProcessOrder(x)
  x <- reorderGraph(x, po)
  ## attr(x,"processOrder")
  
  o <- getOrigin(x)
  
  if(length(o) > 1L)
    stop("The graph has to have a single origin to be transformable ",
         "into a tree.")
  
  ## The graph is coerced into as tree, if necessary.
  if(!isTree(x)) {
    
    d <- attr(x,"dist")
    if(is.null(d))
      d <- graphDist(x)
    
    do <- numeric(ev[2L])
    do[-o] <- d[dst_idx(ev[2L],o)]
    
    ## i=1L
    for(i in 1L:ev[2L]) {
      wh <- which(x$edge[[2L]] == i)
      if(length(wh) > 1L) {
        wh <- wh[-which.min(do[x$edge[[1L]][wh]])]
        
        as.data.frame(
          lapply(x$edge, function(x) x[wh]),
          row.names = attr(x,"elabel")[wh],
        ) -> df
        colnames(df)[1L:2L] <- c("from","to")
        
        
        attr(x,"discarded") <- rbind(attr(x,"discarded"), df)
        
        x <- rm.edge(x, wh)
      }
    }
    
    warning("The graph had to be transformed into a tree.")
  }
  
  sw <- integer(ev[2L])
  tip <- getTerminal(x)
  sw[tip] <- 1L:length(tip)
  sw[-tip] <- (length(tip) + 1L):ev[2L]
  
  structure(
    list(
      edge = cbind(sw[x$edge[[1L]]], sw[x$edge[[2L]]]),
      tip.label = attr(x,"vlabel")[tip],
      node.label = attr(x,"vlabel")[-tip],
      Nnode = ev[2L] - length(tip),
      edge.length = x$edge$distance
    ),
    class = "phylo"
  )
}
#' 
