## **************************************************************************
##
##    (c) 2010-2025 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    ** Phylogenetic Eigenvector Maps (PEM) Main Functions **
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
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with MPSEM. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Phylogenetic Eigenvector Maps
#' 
#' @description Functions to calculate and manipulate Phylogenetic Eigenvector
#' Maps (PEM), which are sets of eigenfunctions describing the structure of a 
#' phylogenetic graph. Each computation function is briefly described in 
#' section \code{Functions} below.
#' 
#' @name PEM-functions
#' 
#' @param x A \code{\link{graph-class}} object or a model matrix of auxiliary
#' trait values to be used alongside the PEM eigenfunctions for modelling the
#' response trait(s) (see argument \code{y} below).
#' @param w A \code{\link{graph-class}} object.
#' @param d A numeric vector of the evolutionary distances (\code{PEMweights})
#' or a character string specifying the row name of the edge table (located in
#' the \code{\link{graph-class}} object's attribute \code{edge}) where the
#' evolutionary distances (edge lengths) can be found.
#' @param a The steepness parameter describing whether changes occur, on
#' average: progressively long edges (a close to 0) or abruptly at vertices (a
#' close to 1; default: \code{0}).
#' @param psi Relative evolution rate along the edges. This parameter only
#' becomes relevant when multiple values are assigned to different portions of
#' the phylogeny (default: \code{1}).
#' @param sp A character string giving the name of a Boolean (i.e., type
#' \code{\link{logical}}) vertex property specifying which of the vertices are
#' species with known traits (see \code{\link{graph-class}} for the details).
#' @param tol A numeric singular value threshold above which a singular vector
#' is retained.
#' @param object A \code{\link{PEM-class}} object.
#' @param y A numeric vector (single trait) or a numeric matrix of response
#' traits.
#' @param lower Lower limit for the \sQuote{L-BFGS-B} optimization algorithm
#' implemented in \code{\link{optim}}.
#' @param upper Upper limit for the \sQuote{L-BFGS-B} optimization algorithm
#' implemented in \code{\link{optim}}.
#' @param tree First parameter of function \code{getGraphLocations}:
#' Phylogenetic tree object with class \sQuote{phylo} (package \link[ape]{ape})
#' containing all species (model and target) used in the study.
#' @param target Name of the target species to extract using the tree
#' \code{tree}.
#' @param gsc The output of \code{getGraphLocations}.
#' 
#' @return The returned value depends on the function:
#' \describe{
#' \item{InflMat}{A binary influence matrix of the graph with as many rows as
#' its number of vertices and as many columns as its number of edges.}
#' \item{PEMweights}{A set of numeric value to be used as weights during PEM
#' calculation.}
#' \item{PEM.build}{A \code{\link{PEM-class}} object.}
#' \item{PEM.updater}{A \code{\link{PEM-class}} object.}
#' \item{PEM.fitSimple}{A \code{\link{PEM-class}} object with embedded fitting
#' parameters.}
#' \item{PEM.forcedSimple}{...}
#' \item{getGraphLocations}{...}
#' \item{getAncGraphLocations}{...}
#' }
#' 
#' @author \packageAuthor{MPSEM} --
#' Maintainer: \packageMaintainer{MPSEM}
#' 
#' @references
#' Guénard, G., Legendre, P., and Peres-Neto, P. 2013. Phylogenetic eigenvector
#' maps: a framework to model and predict species traits. Methods in Ecology 
#' and Evolution. 4: 1120--1131
#' 
#' Makarenkov, V., Legendre, L. & Desdevise, Y. 2004. Modelling phylogenetic
#' relationships using reticulated networks. Zoologica Scripta 33: 89--96
#' 
#' Blanchet, F. G., Legendre, P. & Borcard, D. 2008. Modelling directional
#' spatial processes in ecological data. Ecological Modelling 215: 325--336
#' 
#' @seealso \code{\link{PEM-class}}
#' 
#' @importFrom ape is.rooted drop.tip
#' @importFrom stats optim na.omit
#' @importFrom MASS ginv
#' 
#' @examples
#' ## Synthetic example
#' 
#' ## This example describes the phyogeny of 7 species (A to G) in a tree with 6
#' ## nodes, presented in Newick format, read by function
#' ## read.tree of package ape.
#' 
#' t1 <- read.tree(text=paste(
#'   "(((A:0.15,B:0.2)N4:0.15,C:0.35)N2:0.25,((D:0.25,E:0.1)N5:0.3,",
#'   "(F:0.15,G:0.2)N6:0.3)N3:0.1)N1:0.1;",sep=""))
#' t1
#' summary(t1)
#' 
#' ## Turn the tree into a graph
#' x <- as.graph(t1)
#' 
#' ## Calculate the (binary) influence matrix; E1 to E12 are the tree edges
#' ## Edge E12 comes from the tree origin
#' InflMat(x)
#' InflMat(x)[x$species,]
#' 
#' ## A suite of weighting function profiles:
#' seq(0,1.5,0.01) %>%
#'   plot(y=PEMweights(., a=0), ylim=c(0,1.7), type="l", xlab="distance",
#'        ylab="weight")
#' 
#' seq(0,1.5,0.01) %>%
#'   lines(y=PEMweights(., a=0.5), col="red")
#' 
#' seq(0,1.5,0.01) %>%
#'   lines(y=PEMweights(., a=0.5, psi=1.5), col="green")
#' 
#' seq(0,1.5,0.01) %>%
#'   lines(y=PEMweights(., a=0.9), col="blue")
#' 
#' ## Building phylogenetic eigenvector maps
#' PEM1 <- PEM.build(x)
#' PEM2 <- PEM.build(x, a=0.2)
#' PEM3 <- PEM.build(x, a=1)
#' PEM4 <- PEM.updater(PEM3, a=0.5)
#' 
#' ## Print summary statistics about PEM1
#' print(PEM1)
#' 
## ## Extract the eigenvectors (species A--G, 6 eigenvectors)
## as.data.frame(PEM4)
#' 
#' ## Example of a made up set of trait values for the 7 species
#' y <- c(A=-1.1436265,B=-0.3186166,C=1.9364105,D=1.7164079,E=1.0013993,
#'        F=-1.8586351,G=-2.0236371)
#' 
#' ## Estimate a single steepness parameter for the whole tree
#' PEM.fitSimple(
#'   y = y,
#'   w = x,
#'   d = "distance",
#'   sp = "species",
#'   lower = 0,
#'   upper = 1
#' ) -> PEMfs1
#' 
#' ## Optimisation results:
#' PEMfs1$optim
#' 
#' ## Force neutral evolution over the whole tree:
#' PEMfrc1 <- PEM.forcedSimple(y=y, w=x, d="distance", sp="species", a=0)
#' 
#' ## Steepness parameter forced on each individual edge:
#' edge(PEMfrc1$x)$a
#' 
#' 
#' ## Graph locations for target species X, Y, and Z not found in the original
#' ## data set
#' read.tree(
#'   text = paste(
#'     "((X:0.45,((A:0.15,B:0.2)N4:0.15,(C:0.25,Z:0.2)NZ:0.1)N2:0.05)NX:0.2,",
#'     "(((D:0.25,E:0.1)N5:0.05,Y:0.25)NY:0.25,(F:0.15,G:0.2)N6:0.3)N3:0.1)N1;",
#'     sep=""
#'   )
#' ) -> tree
#' 
#' tree
#' 
#' ## Summary of the structure of the tree:
#' summary(tree)
#' 
#' grloc <- getGraphLocations(tree, target=c("X","Y","Z"))
#' 
#' grloc
#' 
#' PEM.fitSimple(
#'   y = y,
#'   w = grloc$x,
#'   d = "distance",
#'   sp = "species",
#'   lower = 0,
#'   upper = 1
#' ) -> PEMfs2
#' 
#' PEMfs2
#' 
#' ## Same as for PEMfs1$optim
#' PEMfs2$optim
#' 
#' ## Get the PEM scores from the species graph locations:
#' PEMsc1 <- Locations2PEMscores(PEMfs2, grloc)
#' lm1 <- lm(y ~ V_2 + V_3 + V_5, data=PEMfs2)
#' 
#' ## Making prdictions for the species in locations `grloc`
#' ## using linear model `lm1`:
#' ypred <- predict(object=PEMfs2, targets=grloc, lmobject=lm1, interval="none")
#' 
#' ## Removing species X, Y, and Z from the tree in `tpAll`:
#' treeModel <- drop.tip(tree, c("X","Y","Z"))
#' 
#' ## Plot the results
#' layout(t(c(1,1,2)))
#' par(mar=c(6,2,2,0.5)+0.1)
#' plot(treeModel, show.tip.label=TRUE, show.node.label=TRUE, root.edge = TRUE,
#'      srt = 0, adj=0.5, label.offset=0.08, font=1, cex=1.5, xpd=TRUE)
#' edgelabels(paste("E", 1:nrow(treeModel$edge), sep=""),
#'            edge=1:nrow(treeModel$edge), bg="white", font=1, cex=1)
#' points(x=0.20,y=2.25,pch=21,bg="black")
#' lines(x=c(0.20,0.20,0.65), y=c(2.25,0.55,0.55), xpd=TRUE, lty=2)
#' text("X",x=0.69, y=0.55, xpd=TRUE, font=1, cex=1.5)
#' points(x=0.35, y=4.5,pch=21,bg="black")
#' lines(x=c(0.35,0.35,0.6), y=c(4.5,5.47,5.47), xpd=TRUE, lty=2)
#' text("Y", x=0.64, y=5.47, xpd=TRUE, font=1, cex=1.5)
#' points(x=0.35, y=3, pch=21, bg="black")
#' lines(x=c(0.35,0.35,0.55), y=c(3,3.5,3.5), xpd=TRUE, lty=2)
#' text("Z", x=0.59, y=3.5, xpd=TRUE, font=1, cex=1.5)
#' text(c("NX","NY","NZ"), x=c(0.20,0.35,0.35), y=c(2.25,4.5,3)+0.3*c(1,-1,-1),
#'      font=1, cex=1)
#' add.scale.bar(length=0.1, cex=1.25)
#' par(mar=c(3.75,0,2,2)+0.1)
#' plot(x=y, y=1:7, ylim=c(0.45,7), xlim=c(-4,4), axes=FALSE, type="n", xlab="")
#' axis(1, label=c("-4","-2","0","2","4"), at=c(-4,-2,0,2,4))
#' abline(v=0)
#' 
#' ## Plot the observed values
#' points(x=y, y=1:7, xlim=c(-2,2), pch=21, bg="black")
#' text("B)", x=-3.5, y=7, cex=1.5, xpd=TRUE)
#' text("Trait value", x=0, y=-0.5, cex=1.25, xpd=TRUE)
#' 
#' ## Plot the predicted values
#' points(x=ypred, y=c(0.5,5.5,3.5), pch=23, bg="white", cex=1.25)
#' 
#' ## Estimate the ancestral trait values
#' ANCloc <- getAncGraphLocations(x)
#' PEMfsAnc <- PEM.fitSimple(y=y, x=NULL, w=ANCloc$x, d="distance",
#'                           sp="species", lower=0, upper=1)
#' PEMfsAnc$optim
#' 
#' ## Get the PEM scores from the species graph locations:
#' PEManc1 <- Locations2PEMscores(PEMfsAnc, ANCloc)
#' 
#' ## Making predictions for the ancestral species whose locations are found in
#' ## `ANCloc` using the linear model `lm1`:
#' y_anc <- predict(object=PEMfsAnc, targets=ANCloc, lmobject=lm1,
#'                  interval="confidence")
#' 
NULL
#' 
#' @describeIn PEM-functions
#' 
#' PEM Weighting
#' 
#' A power function to obtain the edge weights used during PEM calculation.
#' 
#' @export
PEMweights <- function (d, a = 0, psi = 1) {
  
  nd <- length(d)
  a <- rep(a, length.out = nd)
  psi <- rep(psi, length.out = nd)
  
  .C(
    "PEMweightC",
    as.double(d),
    as.integer(nd),
    as.double(a),
    as.double(psi),
    res = double(nd),
    PACKAGE = "MPSEM"
  )$res
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Building
#' 
#' Calculates a PEM with parameters given by arguments a and psi.
#' 
#' @export
PEM.build <- function(x, d = "distance", sp = "species", a = 0, psi = 1,
                      tol = .Machine$double.eps^0.5) {
  
  if(!inherits(x, "graph"))
    stop("Parameter 'x' must be of class 'graph'")
  
  nv <- nrow(x)
  ne <- nedge(x)
  edge <- edge(x)
  
  if(is.null(edge[[d]]))
    stop("There is no property '",d,"' to be used as edges lengths.")
  
  if(is.null(x[[sp]]))
    stop("There is no property '",sp,"' to indicate species vertices.")
  
  nsp <- sum(x[[sp]])
  
  ### All graphs should be single-rooted
  a <- rep(a, length.out=ne)
  psi <- rep(psi, length.out=ne)
  out <- list(x=x, sp=x[[sp]])
  
  matrix(
    .C("InflMatC",
       ne,
       nv,
       as.integer(edge[[1L]]),
       as.integer(edge[[2L]]),
       B = integer(nv*ne),
       PACKAGE = "MPSEM"
    )$B,
    nrow = nv,
    ncol = ne
  ) -> out[["B"]]
  
  c(
    out,
    .C(
      "PEMbuildC",
      ne = ne,
      nsp = nsp,
      Bc = as.double(out$B[out$sp,]),
      means = double(ne),
      dist = as.double(edge[[d]]),
      a = as.double(a),
      psi = as.double(psi),
      w = double(ne),
      BcW = double(nsp*ne),
      PACKAGE = "MPSEM"
    )
  ) -> out
  
  attr(out$Bc, "dim") <- c(nsp, ne)
  attr(out$BcW, "dim") <- c(nsp, ne)
  
  dimnames(out$Bc) <- dimnames(out$BcW) <-
    list(rownames(x)[x[[sp]]], rownames(edge))
  
  out <- c(out, La.svd(out$BcW, nsp, nsp))
  
  sel <- out$d >= tol
  out$d <- out$d[sel]
  out$u <- out$u[,sel,drop=FALSE]
  out$vt <- out$vt[sel,,drop=FALSE]
  rownames(out$vt) <- colnames(out$u) <- paste("V", 1L:sum(sel), sep="_")
  rownames(out$u) <- rownames(x)[x[[sp]]]
  colnames(out$vt) <- rownames(edge)
  
  structure(
    out,
    class = "PEM"
  )
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Update
#' 
#' Update a PEM with new parameters given by arguments a and psi.
#' 
#' @export
PEM.updater <- function(object, a, psi = 1, tol = .Machine$double.eps^0.5) {
  
  if(!inherits(object, "PEM"))
    stop("Parameter 'x' must be of class 'PEM'")
  
  nsp <- object$nsp
  ne <- object$ne
  object$a <- rep(a, length.out=ne)
  object$psi <- rep(psi, length.out=ne)
  
  .C(
    "PEMupdateC",
    as.integer(ne),
    as.integer(nsp),
    as.double(object$Bc),
    as.double(object$dist),
    as.double(object$a),
    as.double(object$psi),
    w=as.double(object$w),
    BcW = as.double(object$BcW),
    PACKAGE = "MPSEM"
  ) -> out
  
  object$w <- out$w
  object$BcW[] <- out$BcW
  
  out <- La.svd(object$BcW, nsp, nsp)
  sel <- out$d > tol
  object$d <- out$d[sel]
  object$u[] <- out$u[, sel, drop = FALSE]
  object$vt[] <- out$vt[sel, , drop = FALSE]
  
  object
}
#' 
#' @describeIn PEM-functions
#' 
#' Fitting a PEM to Data while Estimating Global Steepness
#' 
#' Fits a PEM to a data set estimating the selection (steepness) parameter using
#' gradient descent. The selection and evolution rate (psi = 1) are assumed to
#' be homogeneous for the whole phylogenetic network.
#' 
#' @export
PEM.fitSimple <- function(y, w, x = NULL, d = "distance", sp = "species",
                          lower = 0, upper = 1, tol = .Machine$double.eps^0.5) {
  
  if(!inherits(w, "graph"))
    stop("Parameter 'w' must be of class 'graph'")
  
  object <- PEM.build(w, d=d, sp=sp, a=lower)
  
  nsp <- sum(w[[sp]])
  
  if(NROW(y) != nsp)
    stop("The number of rows in 'y' must match the number of species in ",
         "phylogenetic graph 'w'")
  
  p <- NCOL(y)
  
  if(!is.matrix(y))
    y <- matrix(y, ncol=1L, dimnames=list(names(y),"y"))
  
  f <- function(par, y, x, object, nsp, p, tol) {
    
    object <- PEM.updater(object, a=par[1L], psi=1, tol=tol)
    Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
    logdetS <- sum(log(object$d^2))
    
    if(is.null(x)) {
      BX <- mean(y)
      res <- y - BX
    } else {
      BX <- ginv(t(x) %*% Sinv %*% x) %*% t(x) %*% Sinv %*% y
      res <- y - x %*% BX
      BX <- c(mean(res), BX)
      res <- res - BX[1L]
    }
    
    dev <- 0
    for(i in 1L:p) {
      S2 <- as.numeric(t(res[,i,drop=FALSE]) %*% Sinv %*% res[,i,drop=FALSE]/nsp)
      dev <- dev + nsp + (nsp - 1L)*log(2*pi) + nsp*log(S2) + logdetS
    }
    
    dev
  }
  
  optim(
    par = 0,
    fn = f,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    y = y,
    x = x,
    object = object,
    nsp = nsp,
    p = p,
    tol = tol
  ) -> opt
  
  if(opt$convergence)
    warning("No optimum found... Message from optim() - ",opt$message,
            ". Status = ",opt$convergence)
  
  object <- PEM.updater(object, a=opt$par[1L], psi=1, tol=tol)
  
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  
  if(is.null(x)) {
    BX <- mean(y)
    res <- y - BX
  } else {
    BX <- MASS::ginv(t(x) %*% Sinv %*% x) %*% t(x) %*% Sinv %*% y
    res <- y - x %*% BX
    BX <- c(mean(res), BX)
    res <- res - BX[1L]
  }
  
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE]) %*% Sinv %*% res[,i,drop=FALSE]/nsp)
  
  object$S2 <- S2
  object$y <- y
  object$optim <- opt
  attr(object$x,"edge")[["a"]] <- rep(opt$par[1L], nrow(attr(object$x,"edge")))
  
  object
}
#' 
#' @describeIn PEM-functions
#' 
#' Fitting a PEM to Data while Forcing Global Steepness
#' 
#' Fits a PEM to a data set forcing a user-provided selection (steepness)
#' parameter. The selection and evolution rate (psi = 1) are assumed to be
#' homogeneous for the whole phylogenetic network.
#' 
#' @export
PEM.forcedSimple <- function(y, w, x = NULL, d = "distance", sp = "species",
                             a = 0, psi = 1, tol = .Machine$double.eps^0.5) {
  
  if(!inherits(w, "graph"))
    stop("Parameter 'w' must be of class 'graph'")
  
  object <- PEM.build(w, d=d, sp=sp, a=a, psi=psi, tol=tol)
  
  nsp <- sum(w[[sp]])
  
  if(NROW(y) != nsp)
    stop("The number of rows in 'y' must match the number of species in ",
         "phylogenetic graph 'w'")
  
  p <- NCOL(y)
  
  if(!is.matrix(y))
    y <- matrix(y, ncol=1L, dimnames=list(names(y), "y"))
  
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  
  if(is.null(x)) {
    BX <- mean(y)
    res <- y - BX
  } else {
    BX <- MASS::ginv(t(x) %*% Sinv %*% x) %*% t(x) %*% Sinv %*% y
    res <- y - x %*% BX
    BX <- c(mean(res), BX)
    res <- res - BX[1L]
  }
  
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE]) %*% Sinv %*% res[,i,drop=FALSE]/nsp)
  
  object$S2 <- S2
  object$y <- y
  object$optim <- list()
  object$optim$par <- c(a, psi)
  
  edge(object$x)$a <- rep(a, length.out=nedge(object$x))
  edge(object$x)$psi <- rep(psi, length.out=nedge(object$x))
  
  object
}
#' 
#' @describeIn PEM-functions
#' 
#' Get Phylogenetic Graph Locations
#' 
#' Takes a phylogenetic tree and a list of species to be removed, and produce a
#' phylogenic graph without these species together with the locations of the
#' removed species on that graph (i.e., the location where the removed species
#' would be found should they be inserted again in the phylogenetic graph).
#' 
#' @export
getGraphLocations <- function(tree, target) {
  
  if(!inherits(tree, "phylo"))
    stop("Argument 'tree' must be of class 'phylo'")
  
  oldnlab <- tree$node.label
  tree$node.label <- newnlab <- paste("N", 1L:tree$Nnode, sep="")
  tpmodel <- drop.tip(tree, target)
  tpmodel$root.edge <- tree$root.edge
  xmodel <- as.graph(tpmodel)
  Bmodel <- InflMat(xmodel)
  
  matrix(
    data = NA,
    nrow = length(target),
    ncol = ncol(Bmodel),
    dimnames = list(target, colnames(Bmodel))
  ) -> loc
  
  dtt <- rep(NA, length(target))
  
  ## Target species need to be striped one-by-one.
  if(length(target))
    for(i in 1L:length(target)) {
      
      tptargetonly <-
        if(length(target) == 1L) tree else drop.tip(tree, target[-i])
      tptargetonly$root.edge <- tree$root.edge
      
      ## The vertices names.
      Vnames <- c(tptargetonly$tip.label, tptargetonly$node.label)
      
      ## The index of the vertex immediately below the target:
      VBidx <- which(tptargetonly$edge[,2L] == which(Vnames == target[i]))
      
      ## VB: the vertex immediately below the target:
      VB <- tptargetonly$edge[VBidx,1L]
      
      ## The name of VB:
      VBName <- Vnames[VB]
      
      ## Distance between VB and VX:
      dBX <- tptargetonly$edge.length[VBidx]
      
      ## VA: the vertex below the one immediately below and...
      VA <- tptargetonly$edge[tptargetonly$edge[,2L] == VB,1L]
      
      ## its name. That vertex may not exist:
      VAName <- Vnames[VA]
      
      ## Distance between VB and VA (if exists):
      dAB <- tptargetonly$edge.length[tptargetonly$edge[,2L] == VB]
      
      ## The vertices above the vertex immediately below and...
      VC <- tptargetonly$edge[tptargetonly$edge[,1L] == VB,2L]
      
      ## their names. 2 or more (in case of tri+ chotomy):
      VCName <- Vnames[VC]
      
      ## Distances above VB:
      dBC <- tptargetonly$edge.length[tptargetonly$edge[,1L] == VB]
      
      VC <- VC[VCName!=target[i]]
      dBC <- dBC[VCName!=target[i]]
      VCName <- Vnames[VC]
      
      ## Strip the target itself from the VC list.
      if(length(VA) == 0L) {
        ## If the target species is beyond the root:
        loc[i,] <- 0
        ## In all case: location = the root.
        dtt[i] <- dBX + min(dBC)
        ## Distance between the off-root species and the remaining ones.
      } else {
        ## If the target species is not beyond the root:
        dtt[i] <- dBX
        ## In all cases: dtt == the distance between the target and the vertex
        ## immediatly below.
        if(length(VC) > 1L) {
          ## When VB is more than dichotomic (several Vertices C):
          loc[i,] <- Bmodel[VBName,]*edge(xmodel)$distance
          ## Coordinates are those of the still-existing vertex B.
        } else {
          ## When VB collapses (it was dichotomic):
          loc[i,] <- Bmodel[VAName,]*edge(xmodel)$distance
          ## Copy the coordinates of vertex VA.
          loc[i,Bmodel[VAName,] != Bmodel[VCName,]] <- dAB
          ## Change the length of that edge connecting VA to the single VC
          ## (0 for VA) to dAB.
        }
      }
    }
  
  if(!is.null(oldnlab)) {
    
    oldnlab <- oldnlab[match(rownames(xmodel),newnlab)]
    rownames(xmodel)[!is.na(oldnlab)] <- na.omit(oldnlab)
  }
  
  list(
    x = xmodel,
    locations = loc,
    LCA2target = dtt
  )
}
#' 
#' @describeIn PEM-functions
#' 
#' Get Ancestral Species Location
#' 
#' Get the location on the phylogenetic graph of the immediate ancestors for a
#' list of species. The species of the list remain in the resulting phylogenetic
#' graph. This function is useful for estimating the ancestral state of a trait.
#' 
#' @export
getAncGraphLocations <- function(x, tree) {
  if(missing(x) && !missing(tree))
    x <- as.graph(tree)
  loc <- InflMat(x)[!x$species,]*edge(x)$distance
  return(list(x = x, locations = loc, LCA2target = numeric(nrow(loc))))
}
#' 
#' @describeIn PEM-functions
#' 
#' PEM Score Calculation
#' 
#' Calculates the scores of an extant or ancestral species on a phylogenetic
#' eigenvector map (i.e., its value on the eigenvectors of the map) from its
#' location on the phylogenetic graph used to build that map.
#' 
#' @export
Locations2PEMscores <- function(object, gsc) {
  
  if(object$ne != ncol(gsc$locations))
    stop("Numbers of edge coordinates mismatch: gsc$locations has ",
         ncol(gsc$locations)," columns while the phylogenetic graph has ",
         object$ne," edges.")
  
  if(is.null(object$y) || is.null(object$optim))
    stop("PEM has no attached data/optimization parameters.")
  
  ntgt <- nrow(gsc$locations)
  nsv <- length(object$d)
  
  out <- list()
  matrix(
    .C(
      "PEMLoc2Scores",
      as.integer(object$ne),
      as.double(object$means*object$w),
      as.integer(ntgt),
      as.double(gsc$locations),
      as.double(object$a),
      as.double(object$psi),
      as.integer(nsv),
      object$d,
      object$vt,
      scores=double(nsv*ntgt),
      PACKAGE = "MPSEM"
    )$scores,
    nrow = ntgt,
    ncol = nsv,
    dimnames = list(rownames(gsc$locations), colnames(object$u))
  ) -> out[["scores"]]
  
  VarianceFactor <- PEMvar(gsc$LCA2target, a=object$optim$par[1L], psi=1)
  matrix(object$S2, nrow=length(object$S2), ncol=1L) %*%
    matrix(VarianceFactor, nrow=1L, ncol=length(VarianceFactor)) -> VarianceFactor
  dimnames(VarianceFactor) <- list(colnames(object$y), rownames(gsc$locations))
  out[["VarianceFactor"]] <- VarianceFactor
  
  out
}
#' 
