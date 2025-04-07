##
### Development script.
##
## rm(list=ls())
compile <- function() {
  library(MPSEM)
  try(dyn.unload("workfile.so"), silent=TRUE)
  system("R CMD SHLIB -o workfile.so workfile.c")
  dyn.load("workfile.so")
  source("fun1.R")
  source("Run_example3/Run_example3.R")
}

compile()

tg <- rownames(gr_ex)[which(gr_ex$species)]
ngrp <- 10L
grpi <- sample(rep(1L:ngrp, length.out = length(tg)))
list(
  obs = df1[gr_ex$species,1L:7L],
  aux = matrix(NA, sum(gr_ex$species), 7L,
               dimnames=list(tg,colnames(df1)[1L:7L])),
  prd = array(NA, c(sum(gr_ex$species),7L,3L),
              list(tg,colnames(df1)[1L:7L],c("fit","lwr","upr")))
) -> pred
rownames(pred$obs) <- tg

## locate(gr_ex, c("V55","V66"))
## locate(gr_ex, c("V58","V78"))

## i=1L
for(i in 1L:ngrp) {
  
  cat("Fold:",i,"\n")
  
  loc <- locate(gr_ex, tg[grpi == i])
  
  PEM(
    loc$x,
    a = PEM1$par()$a,
    mm_a = PEM1$mm$a[!attr(loc$x,"removedEdge"),]
  ) -> PEM_tmp
  
  ## PEM_tmp$S2(
  ##   y = mdat$y[!attr(loc$x,"removedVertex"),][loc$x$species,1L],
  ##   x = mdat$x[!attr(loc$x,"removedVertex"),-1L][loc$x$species,]
  ## )
  
  PEM_tmp$S2(
    y = mdat$y[!attr(loc$x,"removedVertex"),][loc$x$species,],
    x = mdat$x[!attr(loc$x,"removedVertex"),-1L][loc$x$species,]
  )
  
  ## PEM_tmp$S2()
  ## edge(PEM_tmp$graph())
  ## predict(PEM_tmp, loc$location)
  
  ## pemlm(
  ##   log_morpho1~.,
  ##   data = df1[!attr(loc$x,"removedVertex"),][loc$x$species,-(2L:7L)],
  ##   pem = PEM_tmp
  ## ) -> mod
  
  pemlm(
    cbind(log_morpho1,log_morpho2,log_morpho3,log_morpho4,log_morpho5,
          log_morpho6,log_morpho7)~.,
    data = df1[!attr(loc$x,"removedVertex"),][loc$x$species,],
    pem = PEM_tmp
  ) -> mod
  
  ## summary(mod$auxModel)
  ## anova(mod$auxModel)
  ## mod$getIncluded()
  ## mod$getCandidate()
  ## mod$resetPEM()
  mod$forward()
  ## mod$getIncluded()
  ## mod$getCandidate()
  ## summary(mod$pemModel())
  ## anova(mod$pemModel())
  
  ## predict(
  ##   object = mod,
  ##   newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
  ##   newloc = loc$location
  ## )
  
  ## predict(
  ##   object = mod,
  ##   newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
  ##   newloc = loc$location,
  ##   se.fit = TRUE
  ## )
  
  ## predict(
  ##   object = mod,
  ##   newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
  ##   newloc = loc$location,
  ##   interval = "c"
  ## )
  
  predict(
    object = mod$auxModel,
    newdata = df1[match(rownames(loc$location), rownames(gr_ex)),]
  ) -> pred$aux[grpi == i,]
  
  predict(
    object = mod,
    newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
    newloc = loc$location,
    interval = "p"
  ) -> tmp
  
  ## loc$location
  
  ## predict(
  ##   object = mod,
  ##   newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
  ##   newloc = loc$location,
  ##   interval = "c",
  ##   se.fit = TRUE
  ## )
  
  ## predict(
  ##   object = mod,
  ##   newdata = df1[match(rownames(loc$location), rownames(gr_ex)),],
  ##   newloc = loc$location,
  ##   interval = "p",
  ##   se.fit = TRUE
  ## )
  
  pred$prd[grpi == i,,"fit"] <- tmp[,,"fit"]
  pred$prd[grpi == i,,"lwr"] <- tmp[,,"lwr"]
  pred$prd[grpi == i,,"upr"] <- tmp[,,"upr"]
}

rng <- range(pred$obs, pred$prd)
col <- c("red","orange","yellow","green","blue","purple","violet")
psq <- list(aux=numeric(7L), phy=numeric(7L))
par(mar=c(5,5,2,2))
plot(NA, xlim=rng, ylim=rng, xlab = "obs", ylab="pred", asp=1)
abline(0,1)
## i=1L
for(i in 1L:7L) {
  arrows(x0 = pred$obs[,i], x1=pred$obs[,i], y0=pred$prd[,i,2L],
         y1=pred$prd[,i,3L], angle=90, length=0.05, col=col[i], code=3L)
  points(x = pred$obs[,i], y = pred$prd[,i,1L], pch=21L, bg=col[i])
  psq$aux[i] <- Psquare(pred$obs[,i], pred$aux[,i])
  psq$phy[i] <- Psquare(pred$obs[,i], pred$prd[,i,1L])
}
## Psquare(as.matrix(pred$obs), as.matrix(pred$aux))
## Psquare(as.matrix(pred$obs), as.matrix(pred$prd[,,1L]))
legend(x=-3.8, y=4.3, legend=sprintf("%s (%0.2f)", colnames(pred$obs), psq$phy),
       pch=21L, pt.bg=col)
psq











data.frame(
  species = as.logical(c(1,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1)),
  x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5,5,5,2.33),
  y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5,1,0.5,-1),
  row.names = sprintf("V%d",1:16)
) %>%
  st_as_sf(
    coords=c("x","y"),
    crs = NA
  ) %>%
  graph %>%
  add.edge(
    from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,3,6 ,9 ,10,10,3 ,3 ,7 ,9 ),
    to =   c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,6,13,10,12,11,14,15,16,16),
    data = data.frame(
      distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,4.8,3.2,3.5,
                   4.4,2.5,3.4,4.3,3.1,2.2,2.1,0.9,1.0,2.1),
      row.names = sprintf("E%d",1:23)
    )
  ) -> gr1

plot(gr1)

## 1. All the targets must be $species == TRUE, otherwise stop.
## 2. Begin with the targets on terminal vertices.
## 3. Then continue with the targets in median vertices.
## 4. Purge the non-species marked terminal vertices while taking care of
##    following the references of any previous target references when necessary.
## 5. Purge the non-species marked median vertices while taking care of
##    following the references of any previous target references when necessary.
## 6. Remove the species status of any remaining non-removable vertices.

## source("fun1.R")




gr_ex_purged <- purge.terminal(gr_ex)
attr(gr_ex_purged,"removedVertex")
attr(gr_ex_purged,"removedEdge")
gr_ex_purged %>% {.$species[getTerminal(.)]}

gr_ex_purged %<>% purge.median
attr(gr_ex_purged,"removedVertex")
attr(gr_ex_purged,"removedEdge")
gr_ex_purged %>% {.$species[getMedian(.)]}

plot(gr_ex)
plot(gr_ex_purged)

par(mfrow=c(2L,1L), mar=c(0,0,0,0))
tmp <- fun1(gr_ex_purged, getTerminal(gr_ex_purged))
plot(gr_ex_purged)
tmp2 <- tmp$x
rownames(tmp2) <- sprintf("VIDX_%d", 1L:nrow(tmp2))

plot(tmp2)
tmp$target

tr1 <- rtree(25L)
gr1 <- as.graph(tr1)
## attr(purge.terminal(gr1),"removedVertex")
## attr(purge.median(gr1),"removedVertex")

par(mfrow=c(2L,1L), mar=c(0,0,0,0))
for(i in names(getTerminal(gr1))) {
  tmp <- fun1(gr1, i)
  plot(gr1)
  tmp2 <- tmp$x
  plot(tmp2)
  tmp3 <- tmp$target
  if(is.na(tmp3$dist)) {
    tmp3$ref <- rownames(tmp2)[tmp3$ref]
  } else 
    tmp3$ref <- edgenames(tmp2)[tmp3$ref]
  print(tmp3)
  if(is.null(locator(1L)))
    break
}

## Pour la purgation d'un réseau:

## Est retirable durant la purgation d'un réseau:
## 
## 1. un sommet non-espèce et terminale, peu importe le nombre de ses arêtes
##    entrantes.
## 2. un sommet non-espèce et médian. Si la nouvelle arête formée par la
##    consolidation de l'arête entrente avec l'arête sortante joint les mêmes
##    sommets qu'une des arêtes préexistentes, les arêtes sont elle-mème
##    consolidés en une arête. Le calcul de la distance représentée par cette
##    arête doit rendre une distance plus petite que la plus petite des
##    distance. La fonction utilise actuellement l'inverse de la somme de
##    l'inverse des distances comme solution de calcul.
##
## Est retirable durant le calcul des cibles de modélisation:
## 1. un sommet terminal ayant une seule arête entrante.
## 2. un sommet median si la nouvelle arête formée par la consolidation de
##    l'arête entrante avec l'arête sortante ne joint pas les mêmes sommets
##    qu'une des arêtes préexistentes  (règle de non-duplication).
##
## Est retirable suivant le retrait des cibles de modélisation:
## 1. un sommet non-espèce terminal ayant une seule arête entrante.
## 2. un sommet non-espèce median si la nouvelle arête formée par la
##    consolidation de l'arête entrante avec l'arête sortante ne joint pas les
##    mêmes sommets qu'une des arêtes préexistentes (règle de non-duplication).
##





plot(object)

rownames(object)[ttab$ref]
target
ttab


## 



cbind(erm,mask)



vrm


## ttab[!trm,]





## Target function
## Input: a graph + a target list
## Output: a list with a subgraph + a table (id, distance (NA: vertex), LA)

object <- gr_ex2
target <- "V13"










## The new getGraphLocations
##
## A vertex 'v' can be removed if:
## - it is not and extant species ($species[v] == FALSE)
## - it has no more than one incoming edge.
## - it has no more than two outgoing edge.
##
## Any remaining vertex that is non-extant and has only a single incoming edge
## and a single outgoing edge can be purged.
## 
## For the vertices that cannot be removed, the species flag has to be turned
## off.
##


