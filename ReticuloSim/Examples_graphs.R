
## Example graphs.

pop.graph(
  n = 10,
  vertex = list(
    species = rep(TRUE,10),
    x = c(2,3,2,4,3,4,2,1,1,0),
    y = c(-2,1,2,0,-0.5,-2,0,-1,1,0)
  ),
  label = sprintf("V%d",1:10)
) %>%
  add.edge(
    from = c(10,10,9,9,8,8,3,7,7,10,2,2,5,1,4,5),
    to = c(9,8,3,7,7,1,2,2,5,2,1,4,4,4,6,6),
    edge = list(distance=c(1,1,1,1,1,1,1,1,1,4,2,1,1,3,1,1)),
    label = sprintf("E%d",1:16)
  ) -> gr_ex
gr_ex

## Setting the RNG seed to obtain a consistent set of examples:
set.seed(248274957)

## The simplest of all graph: a chain of vertices.
randomGraph(
  NV = 100,
  NC = function(...) 1,
  NP = function(...) 1,
  timestep = function(...) 1,
  maxDist = function(...) NULL
) -> gr_linear

## graphDist(gr_linear)
## 

## A chain of randomly-spaced vertices.
randomGraph(
  NV = 100,
  NC = function(...) 1,
  NP = function(...) 1,
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(...) NULL,
  ts_min = 2,
  ts_max = 5
) -> gr_linear_r

## graphDist(gr_linear_r)

## A dichotomic tree:
randomGraph(
  NV = 100,
  NC = function(...) 2L,        ## A maximum of two children per vertex
  NP = function(...) 1L,        ## Trees have single parents per vertex
  timestep = function(...) 1,
  maxDist = function(...) NULL
) -> gr_dich

## graphDist(gr_dich)

## A dichotomic tree with random branch lengths:
randomGraph(
  NV = 100,
  NC = function(...) 2L,        ## A maximum of two children per vertex
  NP = function(...) 1L,        ## Trees have single parents per vertex
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(...) NULL,
  ts_min = 2,
  ts_max = 5
) -> gr_dich_r

## graphDist(gr_dich_r)

## An other random tree with random number of children per vertex:
randomGraph(
  NV = 100,
  NC = function(nc_mean, ...) 1 + rpois(1L, nc_mean),
  NP = function(...) 1L,        ## Trees have single parents per vertex
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(...) NULL,
  nc_mean = 3,
  ts_min = 2,
  ts_max = 5
) -> gr_dich_r2

## graphDist(gr_dich_r2)

## A dichotomic network with a maximum of two parents per vertex:
randomGraph(
  NV = 100,
  NC = function(...) 2,
  NP = function(...) 2,
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(...) 3,
  ts_min = 2,
  ts_max = 5
) -> gr_net

## graphDist(gr_net)

## A bigger network with more complex embedded functions:
randomGraph(
  NV = 100L,
  NC = function(nc_lambda = 2, ...) 1L + rpois(1L, nc_lambda),
  NP = function(np_lambda = 1, ...) 1L + rpois(1L, np_lambda),
  timestep = function(ts_mu = 0, ts_sigma = 1, ...)
    exp(rnorm(1L, ts_mu, ts_sigma)),
  maxDist = function(md_mu = 2, md_sigma = 1, ...)
    exp(rnorm(1L, md_mu, md_sigma)),
  verbose = FALSE,
  nc_lambda = 4,
  np_lambda = 3,
  ts_mu = 1
) -> gr_big

N <- 100
coords <- cbind(x=runif(N,-1,1), y=runif(N,-1,1))
rownames(coords) <- sprintf("N%d",1:N)
dst <- dist(coords)
gr_dst <- dstGraph(d=dst, th=0.35, origin=15)
rm(N,coords,dst)

ls() %>%
  .[substr(.,1L,3L) == "gr_"] %>%
  paste(collapse=",") %>%
  sprintf("save(%s,file=\"gr_all.rds\")",.) %>%
  parse(text=.) %>%
  eval
