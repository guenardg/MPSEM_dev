
library(MPSEM)

## rm(list=ls())

load(file="gr_all.rds")

isLinear(gr_linear)
isLinear(gr_dich)
isLinear(gr_dst)

isDivergent(gr_linear)
isDivergent(gr_dich)
isDivergent(gr_dst)

isTree(gr_linear)
isTree(gr_dich)
isTree(gr_dst)

as.phylo(gr_linear)
as.phylo(gr_dich)
as.phylo(gr_dst)

coordinates(gr_linear)
coordinates(gr_dich)
coordinates(gr_dst)

coordinates(gr_linear) <- data.frame(1L:attr(gr_linear,"ev")[2L])
coordinates(gr_dich) <- data.frame(1L:attr(gr_dich,"ev")[2L],
                                   1L:attr(gr_dich,"ev")[2L])
coordinates(gr_dst) <- data.frame(1L:attr(gr_dst,"ev")[2L],
                                  1L:attr(gr_dst,"ev")[2L],
                                  1L:attr(gr_dst,"ev")[2L])

coordinates(gr_linear)
coordinates(gr_dich)
coordinates(gr_dst)

x <- gr_linear
x <- gr_dich
x <- gr_dst
coordinates(x)
## coordinates(x) <- NULL
## coordinates(x)

length(x)

labels(x)

labels(x) <- sprintf("N%d",1:length(x))
x

names(x)

names(x) <- c("species","chld","xx")
x

i = c("species","xxxx")
i = c(TRUE,FALSE,TRUE,TRUE)

`[[.graph` <- function(x, i) {
  
  if(is.character(i)) {
    i <- match(i, names(x))
    i <- i[!is.na(i)]
  } else if(is.logical(i)) {
    i <- which(i)
    i <- i[!is.na(i)]
  } else
    i <- i[!is.na(i) && !(i > )]
  
  
  
  
  out <- x
  
  
  
  out$vertex <- out$vertex[i]
  
  out
}

x[[i=c("species","xxxx")]]$vertex


  if(is.character(i)) {
    
    ii <- i
    i <- match(i, names(x))
    
    if(any(is.na(i)))
      return(NULL)
  } else if(is.logical(i)) {
    
    i <- which(i)
  }
  
  if(any(i > length(x$vertex)))
    stop("Subscript(s) out of bound: ",
         paste(i[i > length(x$vertex)], collapse = ", "))
  
  x$vertex <- x$vertex[i]
  
  x
}

x[["species"]]
x[["chld"]]
x[["xx"]]
x[["yy"]]

`[[<-.graph` <- function(x, i, value) {
  
  if(is.null(value)) {
    x$vertex[i] <- NULL
    return(x)
  }
  
  if(inherits(value, "graph")) {
    if(length(value) != length(x))
      stop("The graph objects have different lengths ('value': ", length(value)
           , ", whereas 'x': ", length(x), ")")
    if(!all(labels(value) == labels(x)))
      warning("Label mismatch between graph objects")
    
      
    value <- as.data.frame(value$vertex, row.names=labels(value))
  }
    
  else if(is.numeric(value)) 
    value <- as.data.frame(value, row.names = names(value))
  
  
  if(NROW(value) != length(x))
    stop("The number of values (", NROW(value),
         ") does not correspond to the number of vertices (", length(x) ,")")
  
  if(NCOL(value) != length(i))
    stop("The number of variables (", NCOL(value),
         ") does not correspond to the number of indices (", length(i) ,")")
  
  if(NCOL(value) > 1L && !is.list(value)) {
    stop("Multiple variables have to be provided as a list or data frame")
  } else
    value <- is.data.frame(value, row)
  
  
  
  if(length(i) > 1L) {
    for(j in 1L:length(i))
      x$vertex[[i[j]]] <- value[[j]]
  } else
    x$vertex[[i]] <- value
  
  x
}

x[["xxx"]] <- x[["xx"]]
x[["xxx"]] <- NULL

x[[2L]] <- NULL
x[[10L]] <- NULL
x[["zut"]] <- NULL
x[[c("zut","flute")]] <- data.frame(1L:100L,100L:1L)

x[[3L:9L]] <- NULL





x$vertex

`[.graph` <- function(x, i, j, drop=TRUE) {
  
  value <- as.data.frame(x$vertex, row.names = labels(x))
  
  if(!missing(j)) {
    
    if(is.character(j)) {
      
      jj <- j
      j <- match(j, names(x))
      
      if(any(is.na(j)))
        stop("Unknown vertex descriptor(s): ",
             paste(jj[is.na(j)], collapse=", "))
      
    }
    
    if(any(j > ncol(value)))
      stop("Subscript(s) out of bound: ", paste(j, collapse = ", "))
  }
  
  if(!missing(i)) {
    
    if(is.character(i)) {
      
      ii <- i
      i <- match(i, labels(x))
      
      if(any(is.na(i)))
        stop("Unknown vertex or vertices: ",
             paste(ii[is.na(i)], collapse=", "))
      
    }
    
    if(any(i > nrow(value)))
      stop("Subscript(s) out of bound: ", paste(i, collapse = ", "))
  }
  
  value[if(!missing(i)) i, if(!missing(j)) j, drop]
}

x[,"xx"]




y <- rnorm(100,0,1)

plot.graph <- function(x, y = NULL, bg = c("red","black","blue"), pch = 21L,
                       length = 0.05, pt.cex = 0.75, phylo = FALSE, ...) {
  
  ev <- attr(x,"ev")
  
  ## If vertex coordinates exist, use them to plot the graph.
  if(!(is.null(x$vertex$x) || is.null(x$vertex$y))) {
    
    plot(NA, xlim=range(x$vertex$x), ylim=range(x$vertex$y), type="n", ...)
    
    arrows(
      x0 = x$vertex$x[x$edge[[1L]]],
      x1 = x$vertex$x[x$edge[[2L]]],
      y0 = x$vertex$y[x$edge[[1L]]],
      y1 = x$vertex$y[x$edge[[2L]]],
      length = length,
      ...
    )
    
    points(
      x = x$vertex$x,
      y = x$vertex$y,
      pch = pch,
      bg = if(!is.null(x$vertex$type)) bg[x$vertex$type] else bg[2L],
      cex = pt.cex,
      ...
    )
    
    return(invisible(x))
  }
  
  po <- attr(x,"processOrder")
  if(is.null(po))
    po <- MPSEM:::getProcessOrder(x)
  
  if(isLinear(x)) {
    
    x$vertex$x <- numeric(ev[2L])
    x$vertex$x[po] <- cumsum(c(0,x$edge$distance[match(x$edge[[1L]],po)]))
    
    if(!is.null(y)) {
      x$vertex$y <- y
      if(length(y) != ev[2L]) {
        x$vertex$y <- rep(y, length.out=ev[2L])
        warning("The number of values in 'y' (",length(y),
                ") was not equal to the number of vertices (",ev[2L],")")
      }
    } else
      x$vertex$y <- rep(0,ev[2L])
    
    x$vertex$type <- rep(2L,ev[2L])
    x$vertex$type[getOrigin(x)] <- 1L
    x$vertex$type[getTerminal(x)] <- 3L
    
    plot(x=x$vertex$x[po], y=x$vertex$y[po], type="l", ...)
    
    points(
      x = x$vertex$x[po],
      y = x$vertex$y[po],
      pch = pch,
      bg = bg[x$vertex$type[po]],
      cex = pt.cex,
      ...
    )
    
    return(invisible(x))
  }
  
  if(phylo && isTree(x)) {
    
    phy <- as.phylo(x)
    
    plot(phy, type="cladogram", show.node.label = TRUE, ...)
    
    if(!is.null(y))
      warning("Argument 'y' is not handeled when 'phylo = TRUE'")
    
    return(invisible(phy))
  }
  
  if(isDivergent(x)) {
    
    imat <- InflMat(x)
    if(!is.null(x$edge$distance)) 
      imat <- t(t(imat) * sqrt(x$edge$distance))
    rawc <- svd(scale(imat, scale = FALSE), nu = 2L, nv = 0L)$u
    a <- atan2(rawc[,2L],rawc[,1L])
    a[order(a)] <- seq(-pi, pi, length.out=ev[2L] + 1L)[1L:ev[2L]]
    
    x$vertex$x <- cos(a)
    x$vertex$y <- sin(a)
    
    x$vertex$type <- rep(2L,ev[2L])
    x$vertex$type[getOrigin(x)] <- 1L
    x$vertex$type[getTerminal(x)] <- 3L
    
    plot(NA, xlim=range(x$vertex$x), ylim=range(x$vertex$y), type="n", ...)
    
    arrows(
      x0 = x$vertex$x[x$edge[[1L]]],
      x1 = x$vertex$x[x$edge[[2L]]],
      y0 = x$vertex$y[x$edge[[1L]]],
      y1 = x$vertex$y[x$edge[[2L]]],
      length = length,
      ...
    )
    
    points(
      x = x$vertex$x,
      y = x$vertex$y,
      pch = pch,
      bg = if(!is.null(x$vertex$type)) bg[x$vertex$type] else bg[2L],
      cex = pt.cex,
      ...
    )
    
    return(invisible(x))
  }
  
  warning("Nothing could be plotted.")
  
  return(invisible(x))
}


  
  
  if(is.null(x$vertex$x) || is.null(x$vertex$y))
    x <- MPSEM:::getVertexCoordinate(x)
  
  
  
  
  
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




plot(
  as.phylo(gr_dich_r),
  show.node.label = TRUE,
  tip.color = "red"
)

plot(
  as.phylo(gr_dst),
  type = "cladogram",    ## phylogram, cladogram, fan, unrooted, radial, tidy
  show.node.label = TRUE,
  tip.color = "red"
)







x <- gr_dich_r








cbind(sw[x$edge[[1L]]], sw[x$edge[[2L]]])


## str(rtree(10L))



set.seed(2182955)

## A linear evolutionary sequence with random edge lengths between 2 and 5:
randomGraph(
  NV = 100,
  NC = function(...) 1,
  NP = function(...) 1,
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(...) NULL,
  ts_min = 2,
  ts_max = 5
) -> gr_lin

## attributes(gr_lin)

## Simulate a trait:
simulateTrait(
  x = gr_lin,
  tem = tem[[1]],
  state = 2,
  value = 50,
  a = 1
) -> simTrait

## Display the trait values with the optima
x <- cumsum(c(0,gr_lin$edge$distance))
plot(x=x, y=simTrait$value, type="l", las=1)
lines(x=x, y=c(30,50,80)[simTrait$state], lty=3)

## Simulate an other trait (Brownian motion):
simulateTrait(
  x = gr_lin,
  tem = tem[[2]],
  value = 10,
  a = 1
) -> simTrait

## Show the trait values:
plot(x=x, y=simTrait$value, type="l", las=1)

## A distance-based network:
N <- 100
coords <- cbind(x=runif(N,-1,1), y=runif(N,-1,1))
rownames(coords) <- sprintf("N%d",1:N)
dst <- dist(coords)
gr_dst <- dstGraph(d=dst, th=0.35, origin=15)

## attributes(gr_dst)
graphDist(gr_dst)


## Simulate a trait (Brownian motion) on the distance-based network:
simulateTrait(
  x = gr_dst,
  tem = tem[[2]],
  value = 0,
  a = 1
) -> simTrait

## Display the results of the simulation:
gp <- par(no.readonly = TRUE)

par(mar=c(5,5,1,1))
plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1, las=1, xlab="X", ylab="Y")
points(x=coords[,"x"], y=coords[,"y"], pch=21, cex=abs(simTrait$value)/3,
       bg=gray(0.5 - 0.5*sign(simTrait$value)))

par(gp)





## Function prototype:

simulateTrait(
  x = gr_ex,
  tem = traitEvolSim(
    name = "Trait 1",
    sigma = 1.5,
    alpha = 0.15,
    optima = c(30,50,80),
    transition = matrix(c(NA,0.1,0.0,0.1,NA,0.1,0.0,0.1,NA), 3, 3)
  ),
  state = 2,
  value = 50,
  a = 1
)

simulateTrait(
  x = gr_ex,
  tem = traitEvolSim(
    name = "Trait 2",
    sigma = 2.5
  ),
  value = 0,
  a = 1
)

simulateTrait(
  x = gr_ex,
  tem = traitEvolSim(
    name = "Trait 3",
    alpha = 0.05,
    optima = 15
  ),
  value = 15,
  a = 1
)

simulateTrait(
  x = gr_ex,
  tem = traitEvolSim(
    name = "Trait 4",
    sigma = 2.0,
    alpha = 0.25,
    optima = -25
  ),
  value = -25,
  a = 1
)








simulateSequence(
  x = gr_big,
  Q = DNArate(
    model = "HKY85",
    piGap = 0.25,
    deletionRate = 0.1,
    insertionRate = 0.1,
    pi = c(0.4, 0.4, 0.1, 0.1),
    par = 0.25
  ),
  sqn = drawDNASequence(100, piGap = 0.25, pi = c(0.4,0.4,0.1,0.1)),
  rate = drawEvolRate(100, gamma.shape = 3, gamma.scale = 2e-02),
  a=1
) -> seq

concatenate(seq) %>%
  show.sequence

concatenate(seq, discard="-") %>%
  show.sequence








## This is the matrix holding the sequences:
seq <- matrix(raw(), Ngeneration, Nnucleotide)

## Drawing the initial sequence:
seq[1,] <- drawDNASequence(Nnucleotide, piGap = 0.25, pi = c(0.4,0.4,0.1,0.1))

## Each site has its own mean evolution rate, which are drawn as follows:
erate <- drawEvolRate(Nnucleotide, gamma.shape = 5, gamma.scale = 5e-03)

## Using the Hasegawa, Kishino, and Yano (1985) model:
DNArate(
  model = "HKY85",
  piGap = 0.25,
  deletionRate = 0.1,
  insertionRate = 0.1,
  pi = c(0.4, 0.4, 0.1, 0.1),
  par = 0.25
) -> Q

## Instantiating a molecular evolution models for each site using the single
## change rate matrix, a constant time step of 1 (Ma), and individual mean
## evolution rates drawn previously:
em <- list()
for(j in 1:Nnucleotide)
  em[[j]] <- molEvolSim(Q = Q, step = 1, rho = erate[j])

## These for loops call the $evolve() function to evolve each site for each
## generation as follows:
for(i in 2:Ngeneration)
  for(j in 1:Nnucleotide)
    seq[i, j] <- em[[j]]$evolve(seq[i - 1, j])

## Sequences with the gaps (perfect alignment):
concatenate(seq) %>%
  show.sequence

## Sequences with the gaps removed (prior to multiple sequence alignment):
concatenate(seq, discard="-") %>%
  show.sequence

## Clean up:
rm(Q, em1, tr, Ngeneration, Nnucleotide, seq, erate, em, i, j)

