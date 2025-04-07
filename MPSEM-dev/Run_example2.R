
set.seed(1234567L)

randomGraph(
  NV = 200L,
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
## plot(gr_ex)

gr_ex %<>% purge.terminal
gr_ex %<>% purge.median

## gr_ex

## plot(gr_ex, cex.min=3, cex.lab=0.5)

vtx <- c("V2","V5")

mm <- isUnderVertex(gr_ex, vtx)

b_alpha <- c(-4,2,-4.2)
alpha <- exp(as.numeric(cbind(1, mm) %*% b_alpha))

b_sigma <- c(-0.54,-0.21)
sigma <- exp(as.numeric(mm %*% b_sigma))

reg <- unique(cbind(alpha, sigma))

edge(gr_ex)$tem <- NA

for(i in 1L:nrow(reg))
  edge(gr_ex)$tem[which((alpha == reg[i,1L]) & (sigma == reg[i,2L]))] <- i

tem <- list()

for(i in 1L:nrow(reg))
  traitEvolSim(
    name = sprintf("Trait%d", i),
    sigma = reg[i,"sigma"],
    alpha = reg[i,"alpha"],
    optima = c(10,30),
    transition = matrix(c(NA,0.025,0.05,NA), 2L, 2L)
  ) -> tem[[i]]

tr_ex <- simulateTrait(gr_ex, tem, state=1L, value=20, a=1)
