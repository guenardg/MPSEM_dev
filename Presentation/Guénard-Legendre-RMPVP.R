
## Examples

## rm(list=ls())

## Examples from the MPSEM package

library(MPSEM)
library(magrittr)
library(MASS)
library(RSQLite)
source("RMPVP-aux.R")

if(FALSE) {
  SQLite() %>%
    dbConnect(
      "OUsim.sqlite"
    ) -> OUsim
  ## dbListTables(OUsim)
  ## dbDisconnect(OUsim)
  
  dbWriteTable(
    conn = OUsim,
    name = "time",
    value = data.frame(time = seq(0,5.0,0.01))
  )
  ## dbReadTable(OUsim,"time")
  
  data.frame(
    alpha = c(0,0.1,0.2,0.5,1,2,5,10),
    sigma = 1,
    from = 0,
    opt = 0
  ) %>%
    dbWriteTable(
      conn = OUsim,
      name = "cnd",
      value = .
    )
  ## dbReadTable(OUsim,"cnd")
  
  OUsim %>%
    dbExecute(
      "CREATE TABLE res (
         cnd               INTEGER NOT NULL,
         sim               NITEGER NOT NULL,
         time              INTEGER NOT NULL,
         value             REAL,
         FOREIGN KEY(cnd)  REFERENCES cnd(rowid),
         FOREIGN KEY(time) REFERENCES time(rowid)
       )"
    )
  
  OUsim %>%
    dbExecute(
      "CREATE TABLE kdmap (
         cnd              INTEGER NOT NULL,
         map              BLOB,
         FOREIGN KEY(cnd) REFERENCES cnd(rowid)
       )"
    )
  
  OUsim %>%
    dbGetQuery(
      "SELECT ROWID, * FROM cnd"
    ) -> cnd
  
  OUsim %>%
    dbGetQuery(
      "SELECT ROWID, * FROM time"
    ) -> x
  
  ## i=1L
  for(i in cnd$rowid) {
    ## j=1L
    for(j in 1L:1000L)
      OUEvolve(
        x = x$time,
        from = cnd$from[i],
        opt = cnd$opt[i],
        alpha = cnd$alpha[i],
        sigma = cnd$sigma[i]
      ) %>%
      data.frame(
        cnd = i,
        sim = j,
        time = x$rowid,
        value = .
      ) %>%
      {sprintf("(%d,%d,%d,%f)",.$cnd,.$sim,.$time,.$value)} %>%
      paste(collapse=",\n") %>%
      sprintf(
        "INSERT INTO res (cnd,sim,time,value) VALUES %s",
        .
      ) %>%
      dbExecute(OUsim,.)
  }
  rm(i,j)
  ##
  dbExecute(OUsim, "CREATE INDEX res_cnd ON res(cnd)")
  dbExecute(OUsim, "CREATE INDEX res_sim ON res(sim)")
  dbExecute(OUsim, "CREATE INDEX res_time ON res(time)")
  
  OUsim %>%
    dbGetQuery(
      sprintf(
        "SELECT value
         FROM res
         WHERE cnd = %d AND sim = %d
         ORDER BY time",
        1L,
        1L
      )
    )
    
    
    
    
    rng <- max(abs(OUsims$res[[i]]))*c(-1,1)
    
    kde2d(
      x = rep(OUsims$x,OUsims$nsim),
      y = as.numeric(OUsims$res[[i]]),
      h = c(1,2),
      n = c(200,200),
      lims = c(range(OUsims$x), rng)
    ) -> OUsims$kdmap[[i]]
  }
}








  ## j=1L
  for(j in 1L:nsim) {
    image(tmp, las=1L, col=grey(seq(1,0,-0.001)))
    lines(x=x, y=res[[i]][,j])
    if(is.null(locator(1L))) break
  }
    
  
   
  
}
##




par(mar=c(4.5,4.5,2,2))
plot(NA, xlim=range(y[,1L]), ylim=max(abs(y[,-1L]))*c(-1,1), xlab="Time",
     ylab="Trait value", las=1L)
## i=2L
for(i in 2L:ncol(y)) {
  lines(x=y[,1L], y=y[,i])
}




 %>%
  cbind(
    x = .,
    y0 = OUEvolve(., from=0, opt=0, alpha=0, sigma=1),
    y01 = OUEvolve(., from=0, opt=0, alpha=0.1, sigma=1),
    y02 = OUEvolve(., from=0, opt=0, alpha=0.2, sigma=1),
    y05 = OUEvolve(., from=0, opt=0, alpha=0.5, sigma=1),
    y1 = OUEvolve(., from=0, opt=0, alpha=1, sigma=1),
    y2 = OUEvolve(., from=0, opt=0, alpha=2, sigma=1),
    y5 = OUEvolve(., from=0, opt=0, alpha=5, sigma=1),
    y10 = OUEvolve(., from=0, opt=0, alpha=10, sigma=1)
  ) -> y
##
par(mar=c(4.5,4.5,2,2))
plot(NA, xlim=range(y[,1L]), ylim=max(abs(y[,-1L]))*c(-1,1), xlab="Time",
     ylab="Trait value", las=1L)
## i=2L
for(i in 2L:ncol(y)) {
  lines(x=y[,1L], y=y[,i])
}



## x = seq(0,2.1,0.1)
## from = 0.5
## opt = 0
## alpha = 1
## sigma = 1



edges <- c(B=3.1, C=2.7, D=3.4)
opt <- c(B=0, C=-1, D=1)
alpha <- c(B=0, C=1, D=1)
sigma <- c(B=3, C=2, D=3)
grain <- 0.001

sprintf(
  "((C:%f,D:%f)B:%f)A;",
  edges["C"],
  edges["D"],
  edges["B"]
) %>%
  read.tree(text=.) -> ex1

trait <- list()

seq(0, edges["B"], grain) %>%
  cbind(
    x = .,
    y = OUEvolve(., from=0, opt=opt["B"], alpha=alpha["B"], sigma=sigma["B"])
  ) -> trait[["B"]]

seq(0, edges["C"], grain) %>%
  {. + tail(trait[["B"]][,1L],1L)} %>%
  cbind(
    x = .,
    y = OUEvolve(., from=tail(trait[["B"]][,"y"],1L), opt=opt["C"],
                 alpha=alpha["C"], sigma=sigma["C"])
  ) -> trait[["C"]]

seq(0, edges["D"], grain) %>%
  {. + tail(trait[["B"]][,1L],1L)} %>%
  data.frame(
    x = .,
    y = OUEvolve(., from=tail(trait[["B"]][,"y"],1L), opt=opt["D"],
                 alpha=alpha["D"], sigma=sigma["D"])
  ) -> trait[["D"]]

## ylim <- range(trait[["B"]][,"y"],trait[["C"]][,"y"],trait[["D"]][,"y"], opt)
ylim=c(-5,5)

par(mfrow=c(2L,1L), mar=c(1,4.1,1,2.1))

plot(ex1, show.node.label = TRUE, edge.color = c("black","blue","red")) -> pl1

par(mar=c(4,4.1,1,2.1))

trait[["B"]] %>%
  {plot(x=.[,"x"], y=.[,"y"], type="l", xlim=pl1$x.lim, ylim=ylim, las=1L,
        xlab="Time", ylab="Trait value")}
if(alpha["B"])
  lines(x=range(trait[["B"]][,"x"]), y=rep(opt["B"],2L), lty=3L)

trait[["C"]] %>%
  {lines(x=.[,"x"], y=.[,"y"], col="blue")}
if(alpha["C"])
  lines(x=range(trait[["C"]][,"x"]), y=rep(opt["C"],2L), lty=3L, col="blue")

trait[["D"]] %>%
  {lines(x=.[,"x"], y=.[,"y"], col="red")}
if(alpha["D"])
  lines(x=range(trait[["D"]][,"x"]), y=rep(opt["D"],2L), lty=3L, col="red")









matrix(
  data = c(0.7,0.2,0.2,0.2,0.7,0.1,0.1,0.1,0.7),
  nrow = length(opt),
  ncol = length(opt),
  dimnames = list(from=opt, to=opt)
) -> transit

nsp <- 25L

rtree(
  n = nsp,
  tip.label = paste("Species", 1L:nsp, sep="")
) -> tree2

paste("N", 1L:tree2$Nnode, sep="") -> tree2$node.label

plot(tree2, show.node.label = TRUE)

EvolveOptimMarkovTree(
  tp = tree2,
  tw = transit,
  p = 10,
  anc = 2
) -> reg

TraitOUsimTree(
  tp = tree2,
  a = 0,
  sigma = 1,
  opt = opt[reg[,1L]],
  p = 10L
) -> y1

TraitOUsimTree(
  tp = tree2,
  a = 1,
  sigma = 1,
  opt = opt[reg[,1L]],
  p = 10L
) -> y2

TraitOUsimTree(
  tp = tree2,
  a = 10,
  sigma = 1,
  opt = opt[reg[,1L]],
  p = 10L
) -> y3

displayOUprocess <- function(tp, trait, regime, mvalue) {
  layout(matrix(1:2,1,2))
  n <- length(tp$tip.label)
  ape::plot.phylo(
    x=tp, show.tip.label=TRUE, show.node.label=TRUE, root.edge=FALSE,
    direction="rightwards", adj=0,
    edge.color=rainbow(length(trait))[regime[tp$edge[,2L]]]
  )
  plot(
    y=1:n, x=mvalue[1L:n], type="b", xlim=c(-5,5), ylab="", xlab="Trait value",
    yaxt="n", bg=rainbow(length(trait))[regime[1L:n]], pch=21L
  ) 
  text(
    trait[regime[1L:n]], y=1:n, x=5, col=rainbow(length(trait))[regime[1L:n]]
  )
  abline(v=0)
  invisible(NULL)
}

displayOUprocess(tree2, opt, reg[,1L], y1[,1L])
displayOUprocess(tree2, opt, reg[,1L], y2[,1L])
displayOUprocess(tree2, opt, reg[,1L], y3[,1L])

grain <- 0.1
vertices <- c(A=2.7, B=3.4, C=3.1)
edges <- list(A="A000:0.1", B="B000:0.1", C="(%s,%s)C000:0.1")

## i="A"
for(i in names(edges)) {
  ## j=1L
  for(j in 1L:(round(vertices[i]/grain) - 1L))
    edges[[i]] <- sprintf("(%s)%s%03d:%s",edges[[i]], i, j, grain)
}
edges$C %>%
  substr(
    start = nchar(.) - 6L,
    stop = nchar(.) - 4L
  ) %>%
  as.integer %>%
  {. + 1L} %>%
  sprintf("(%%s)C%03d;",.) %>%
  sprintf(sprintf(edges$C, edges$A, edges$B)) %>%
  read.tree(text=.) -> ex1

par(mfrow=c(3L,1L))
plot(ex1)#, show.node.label = TRUE)

ex1 %>%
  Phylo2DirectedGraph -> dg1
## str(dg1)
## attr(dg1,"vlabel")
reg <- integer(attr(dg1,"ev")[2L])
reg[substr(attr(dg1,"vlabel"),1L,1L) == "C"] <- 2L
reg[substr(attr(dg1,"vlabel"),1L,1L) == "A"] <- 1L
reg[substr(attr(dg1,"vlabel"),1L,1L) == "B"] <- 3L

TraitOUsimTree(
  tp = ex1,
  a = 10,
  sigma = 1,
  opt = opt[reg],
  p = 10L
) -> y4

rownames(y4) %>%
  {data.frame(
    edge = substr(.,1L,1L),
    order = as.integer(substr(.,2L,4L))
  )} -> eid

c(
  y4[eid$edge=="C",1L][order(eid$order[eid$edge=="C"], decreasing = TRUE)],
  y4[eid$edge=="B",1L][order(eid$order[eid$edge=="B"], decreasing = TRUE)]
  ) -> y41

c(
  y4[eid$edge=="C",1L][order(eid$order[eid$edge=="C"], decreasing = TRUE)],
  y4[eid$edge=="A",1L][order(eid$order[eid$edge=="A"], decreasing = TRUE)]
  ) -> y42

xlim <- c(0,max(grain*length(y41),grain*length(y42)))

plot(x=seq(0, grain*length(y41), length.out = length(y41)), y=y41, xlim=xlim)
plot(x=seq(0, grain*length(y42), length.out = length(y42)), y=y42, xlim=xlim)




ex1$edge


dg1$vertex$species
dg1$edge$distance


EvolveOptimMarkovTree(
  tp = ex1,
  tw = transit,
  p = 1,
  anc = 2
) %>% length


ex1$edge
ex1$edge.length
ex1$Nnode
ex1$node.label
ex1$tip.label

ex1 %>%
  Phylo2DirectedGraph -> dgr1
## str(dgr1)
dgr1$vertex$species[] <- TRUE

dgr1 %>%
  PEMInfluence -> infl1




