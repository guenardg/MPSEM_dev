
set.seed(1234567L)

randomGraph(
  NV = 50L,
  NC = function(lambda_child, ...) 1 + rpois(1, lambda_child),
  NP = function(lambda_parent, ...) 1 + rpois(1, lambda_parent),
  timestep = function(ts_min, ts_max, ...) runif(1, ts_min, ts_max),
  maxDist = function(max_anc, ...) runif(1, 0, max_anc),
  mean = "harmonic",
  a = 1,
  lambda_child = 3.0,
  lambda_parent = 3.0,
  ts_min = 5,
  ts_max = 12,
  max_anc = 13,
  verbose = FALSE
) -> gr_ex

## gr_ex

## plot(gr_ex, cex.min=3, cex.lab=0.5)

vtx <- c("V2","V5")

mm <- isUnderVertex(gr_ex, vtx)

if(FALSE) {
  
  plotv <- function(x, mm, v) {
    
    tmp <- integer(nrow(x))
    tmp[edge(x)[[1L]][as.logical(mm[,v])]] <- 1
    tmp[edge(x)[[2L]][as.logical(mm[,v])]] <- 1
    
    plot(x, y = tmp, bg=c("red","blue"))
    
    invisible(NULL)
  }
  
  ## plotv(gr_ex, mm, vtx[1L])
  ## plotv(gr_ex, mm, vtx[2L])
  
  tmp <- factor(rep("V1", nrow(gr_ex)), levels=c("V1",vtx))
  
  for(v in vtx) {
    tmp[edge(gr_ex)[[1L]][as.logical(mm[,v])]] <- v
    tmp[edge(gr_ex)[[2L]][as.logical(mm[,v])]] <- v
  }
  
  plot(gr_ex, y=as.numeric(tmp), bg=rainbow(7L))
  rm(plotv, tmp, v)
}

b_alpha <- c(-4,2,-4.2)
alpha <- exp(as.numeric(cbind(1, mm) %*% b_alpha))

b_sigma <- c(-0.54,-0.21)
sigma <- exp(as.numeric(mm %*% b_sigma))

reg <- unique(cbind(alpha, sigma))

edge(gr_ex)$tem <- NA

# i=1L
for(i in 1L:nrow(reg))
  edge(gr_ex)$tem[which((alpha == reg[i,1L]) & (sigma == reg[i,2L]))] <- i

tem <- list()

# i=1L
for(i in 1L:nrow(reg))
  traitEvolSim(
    name = sprintf("Trait%d", i),
    sigma = reg[i,"sigma"],
    alpha = reg[i,"alpha"],
    optima = c(10,20,30),
    transition = matrix(c(NA,0.1,0.2,0.1,NA,0.1,0.2,0.1,NA), 3L, 3L)
  ) -> tem[[i]]

## all(reg[edge(gr_ex)$tem,1L] == alpha)
## all(reg[edge(gr_ex)$tem,2L] == sigma)

tr_ex <- simulateTrait(gr_ex, tem, state=2L, value=20, a=1)

## plot(gr_ex, y=tr_ex$state, bg=c("red","green","blue"))
## plot(gr_ex, y=tr_ex$value, bg=head(rainbow(1200),1000))

## range(sigma)
## range(alpha)

## data.frame(optim=c(10,20,30)[tr_ex$state], value=tr_ex$value)

## x <- PEM.fitSimple(y = tr_ex$value[gr_ex$species], w=gr_ex)
## gsc <- getAncGraphLocations(gr_ex, gr_ex)
## loc <- Locations2PEMscores(x, gsc)

## y <- matrix(NA, length(tr_ex$state), ncol(x$u))
## y[gr_ex$species,] <- x$u
## y[!gr_ex$species,] <- loc$scores

## plot(gr_ex, y=y[,1L], bg=head(rainbow(1200),1000))

PEM2(
  gr_ex,
  a =   c(-4.00,  2.00, -4.20),
  psi = c(-0.54, -0.21       ),
  mm_a = cbind(CTE=1, mm),
  mm_psi = mm
) -> pem1

## pem1$d
## pem1$sp
## pem1$tol
## pem1$mm
## pem1$nsp
## pem1$B

## pem1$print()
## pem1

## pem1$graph()
## edge(pem1$graph())$a
## edge(pem1$graph())$psi

## pem1$pem()
## as.matrix(pem1)
## as.data.frame(pem1)

## pem1$par()

pem1$update(
  a =   c(-4.00,  2.10, -4.20),
  psi = c(-0.54, -0.21       ),
)

## pem1$par()

## edge(pem1$graph())$a
## edge(pem1$graph())$psi

## pem1$var()
## pem1$inv_var()
## pem1$logdet()

tr_ex2 <- simulateTrait(gr_ex, tem, state=1L, value=10, a=1)

y <- cbind(tr_ex$value[gr_ex$species],tr_ex2$value[gr_ex$species])
x <- matrix(rnorm(2*nrow(y),0,1), ncol=2L)

## pem1$dev(y[,1L,drop=FALSE], NULL)
## pem1$dev(y[,1L,drop=FALSE], x[,1L,drop=FALSE])
## pem1$dev(y[,1L,drop=FALSE], x[,2L,drop=FALSE])
## pem1$dev(y[,1L,drop=FALSE], x)

## pem1$dev(y, NULL)
## pem1$dev(y, x[,1L,drop=FALSE])
## pem1$dev(y, x[,2L,drop=FALSE])
## pem1$dev(y, x)

## pem1$S2(y[,1L,drop=FALSE], NULL)
## pem1$S2(y[,1L,drop=FALSE], x[,1L,drop=FALSE])
## pem1$S2(y[,1L,drop=FALSE], x[,2L,drop=FALSE])
## pem1$S2(y[,1L,drop=FALSE], x)

## pem1$S2(y, NULL)
## pem1$S2(y, x[,1L,drop=FALSE])
## pem1$S2(y, x[,2L,drop=FALSE])
## pem1$S2(y, x)

PEM2(
  gr_ex,
  a =   c(0, 0, 0),
  psi = c(0      ),
  mm_a = cbind(CTE=1, mm),
  mm_psi = mm[,1L,drop=FALSE]
) -> pem2

opt <- evolution.model.PEM2(pem2, y=y, x=x)

## pem2$par()$a == opt$par$a
## pem2$par()$psi == opt$par$psi
## pem2$dev(y, x) == opt$value

data.frame(
  species = as.logical(c(1,0,1,0,1,0,0,0,1,1,1,1,1,1,1)),
  type = c(2,2,3,1,2,2,2,2,2,2,3,3,3,3,3),
  x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5,5,5),
  y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5,1,0.5),
  row.names = sprintf("V%d",1:15)
) %>%
  st_as_sf(
    coords=c("x","y"),
    crs = NA
  ) %>%
  graph %>%
  add.edge(
    from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,3,6 ,9 ,10,10,3 ,3 ),
    to =   c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,6,13,10,12,11,14,15),
    data = data.frame(
      distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,
                   4.8,3.2,3.5,4.4,2.5,3.4,4.3,3.1,2.2,2.1,0.9),
      row.names = sprintf("E%d",1:21)
    )
  ) -> gr_ex2

saveRDS(gr_ex2, file="gr_ex2.rds")

tr2_ex <- simulateTrait(gr_ex2, tem[[1L]], state=2L, value=20, a=1)
tr2_ex2 <- simulateTrait(gr_ex2, tem[[1L]], state=1L, value=10, a=1)

y2 <- cbind(tr2_ex$value[gr_ex2$species],tr2_ex2$value[gr_ex2$species])
x2 <- matrix(rnorm(2*nrow(y2),0,1), ncol=2L)

PEM2(
  gr_ex2,
  a =   c(0),
  psi = c()
) -> pem3

opt2 <- evolution.model.PEM2(pem3, y=y2, x=x2)





