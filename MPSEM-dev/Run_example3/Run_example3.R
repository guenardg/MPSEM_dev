
if(FALSE) {
  
  library(MPSEM)
  
  set.seed(1234567L)
  
  randomGraph(
    NV = 100L,
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
  
  b_alpha <- c(0.3,-1.1,0.5)
  alpha <- exp(as.numeric(cbind(1, mm) %*% b_alpha))
  
  b_sigma <- c(0,0)
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
  
  tr_ex <- list(state=matrix(NA,100L,12L), value=matrix(NA,100L,12L))
  
  for(i in 1L:12L) {
    tmp <- simulateTrait(gr_ex, tem, state=1L, value=20, a=1)
    tr_ex$state[,i] <- c(10,30)[tmp$state]
    tr_ex$value[,i] <- tmp$value
    s <- sqrt(mean(tr_ex$value[,i]^2))
    tr_ex$value[,i] %<>% {./s}
    tr_ex$state[,i] %<>% {./s}
  }
  rm(i,tmp,s)
  
  ## tr_ex$state
  ## tr_ex$value
  
  data.frame(
    habitat = tr_ex$value[,12L] %>%
      cut(
        breaks = c(-Inf,quantile(.,c(1/4,1/2,3/4)),Inf),
        labels = c("Mountain","Prairie","Forest","Wetland")
      ),
    crest_color = tr_ex$value[,11L] %>%
      cut(
        breaks = c(-Inf,quantile(.,c(1/3,2/3)),Inf),
        labels = c("brown","orange","red")
      ),
    gape_width = 3*tr_ex$value[,10L],
    leg_length = 2 + 3*tr_ex$value[,9L],
    nb_spike = ceiling( 10 + 3*tr_ex$value[,8L])
  ) -> indep
  
  mdat <- model.data(~., data=indep)$x
  
  matrix(
    rnorm(ncol(mdat)*7L, 0, 1),
    ncol(mdat), 7L,
    dimnames = list(colnames(mdat), NULL)
  ) -> B
  
  mdat %*% B %>%
    {. %*% diag(1/sqrt(colMeans(.^2)))} -> eff
  
  ## colMeans(tr_ex$value[,1L:7L]^2)
  ## colMeans(eff^2)
  f <- function(x,y,p) p*x + (1 - p)*y
  data.frame(
    log_morpho1 = -0.3 + 2*f(tr_ex$value[,1L], eff[,1L], 0.55),
    log_morpho2 = -0.1 + 2*f(tr_ex$value[,2L], eff[,2L], 0.50),
    log_morpho3 = -1 + 3*f(tr_ex$value[,3L], eff[,3L], 0.55),
    log_morpho4 = 2 - 1.9*f(tr_ex$value[,4L], eff[,4L], 0.60),
    log_morpho5 = 0 + 2.1*f(tr_ex$value[,5L], eff[,5L], 0.55),
    log_morpho6 = 2.3 - 3.2*f(tr_ex$value[,6L], eff[,6L], 0.50),
    log_morpho7 = -1 + 2*f(tr_ex$value[,7L], eff[,7L], 0.45)
  ) -> dep
  
  df1 <- cbind(dep, indep)
  
  lm(
    cbind(log_morpho1,log_morpho2,log_morpho3,log_morpho4,log_morpho5,
          log_morpho6,log_morpho7)~.,
    data = df1
  ) -> lm1
  
  lm1$rank
  summary(lm1)
  anova(lm1)
  coef(lm1)
  
  PEM(
    gr_ex,
    a = c(0,0,0),
    mm_a = cbind(1, mm)#,
    # psi = c(0,0),
    # mm_psi = mm
  ) -> PEM1
  
  model.data(
    cbind(log_morpho1,log_morpho2,log_morpho3,log_morpho4,log_morpho5,
          log_morpho6,log_morpho7)~.,
    data = df1
  ) -> mdat
  
  PEM1 %>%
    evolution.model(
      y = mdat$y,
      x = mdat$x
    ) -> evmod1
  
  
  ## PEM1$par()
  ## PEM1$S2()
  ## edge(PEM1$graph())
  ## as.data.frame(PEM1)
  
  save(gr_ex, vtx, mm, df1, mdat, PEM1, evmod1, file="../Example3.rda")
  ## save(gr_ex, vtx, mm, df1, mdat, PEM1, evmod1, file="Example3.rda")
  rm(list=ls())
  load(file="../Example3.rda")
} else
  load(file="Example3.rda")
