##
## rm(list=ls())
compile <- function() {
  library(MPSEM)
  ## try(dyn.unload("workfile.so"), silent=TRUE)
  ## system("R CMD SHLIB -o workfile.so workfile.c")
  ## dyn.load("workfile.so")
  ## source("workfile.R")
  ## source("fun1.R")
}

compile()

data(perissodactyla,package="caper")
knitr::kable(perissodactyla.data)

par(mar=c(2,2,2,2))
plot(perissodactyla.tree)
par(mar=c(5,4,4,2))

knitr::kable(perissodactyla.data)

match(
  perissodactyla.tree$tip.label,
  perissodactyla.data[,1L]
) -> spmatch

drop.tip(
  perissodactyla.tree,
  perissodactyla.tree$tip.label[is.na(spmatch)]
) -> perissodactyla.tree

cbind(perissodactyla.tree$tip.label, perissodactyla.data[,1L])

match(
  perissodactyla.tree$tip.label,
  perissodactyla.data[,1L]
) -> spmatch

perissodactyla.data[spmatch,] -> perissodactyla.data

all(perissodactyla.tree$tip.label == perissodactyla.data[,1L])

perissodactyla.data[,1L] -> rownames(perissodactyla.data)
perissodactyla.data[,-1L] -> perissodactyla.data

knitr::kable(perissodactyla.data)

par(mar=c(4.5,4.5,1,7) + 0.1)
d <- seq(0, 2, length.out=1000)
a <- c(0,0.33,0.67,1,0.25,0.75,0)
psi <- c(1,1,1,1,0.65,0.65,0.4)
cc <- c(1,1,1,1,1,1,1)
ll <- c(1,2,2,2,3,3,3)
trial <- cbind(a, psi)
colnames(trial) <- c("a","psi")
ntrials <- nrow(trial)
nd <- length(d)
matrix(
  NA,
  ntrials,
  nd,
  dimnames=list(paste("a=", trial[,"a"], ", psi=", trial[,"psi"], sep=""),
                paste("d=", round(d,4), sep=""))
) -> w
for(i in 1:ntrials)
  w[i,] <- PEMweights(d, trial[i,"a"], trial[i,"psi"])
plot(NA, xlim=c(0,2), ylim=c(0,1.6), ylab="Weight", xlab="Distance", axes=FALSE)
axis(1, at=seq(0,2,0.5), label=seq(0,2,0.5))
axis(2, las=1)
text(expression(paste(~~~a~~~~~~~psi)),x=2.2,y=1.57,xpd=TRUE,adj=0)
for(i in 1:ntrials) {
  lines(x=d, y=w[i,], col=cc[i], lty=ll[i])
  text(paste(sprintf("%.2f", trial[i,1]), sprintf("%.2f",trial[i,2]), sep="  "),
       x=rep(2.2,1), y=w[i,1000], xpd=TRUE, adj=0)
}
rm(d,a,psi,cc,ll,trial,ntrials,nd,w,i)

gr1 <- as.graph(perissodactyla.tree)

str(gr1)

perissodactyla.PEM <- PEM(gr1)
perissodactyla.PEM

perissodactyla.tree -> tree
sprintf("N%s",1L:tree$Nnode) -> tree$node.label

par(mar=c(2,2,2,2))
plot(tree, show.node.label=TRUE)

edgelabels(
  sprintf("E%d",1L:nrow(tree$edge)),
  edge=1L:nrow(tree$edge),
  bg="white",
  cex=0.75
)

cbind(
  isUnderVertex(gr1, "N9"),
  isUnderEdge(gr1,"E11")
) -> mm

## Calculate the PEM:
PEM(
  x = gr1,
  a = c(-6,5,8),            ## steepness sub model parameters
  psi = c(-1,1),            ## evolution rate sub model parameters
  mm_a = cbind(1, mm),      ## model matrix: steepness sub model
  mm_psi = mm               ## model matrix: evolution rate sub model
) -> perissodactyla.PEM

## Show the PEM object:
perissodactyla.PEM

## Extract the graph from within the PEM-class object:
gr <- perissodactyla.PEM$graph()

## Show the steepness and evolution rate at the edges:
round(edge(gr)[,c("a","psi")],3)

tmp <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(1.1,1.1,2.6,0.1))

## Singular vectors are extracted using the as.matrix method:
perissodactyla.U <- as.matrix(perissodactyla.PEM)

plot(perissodactyla.tree, x.lim=60, cex=1.0)

par(mar=c(1.1,0.1,2.6,1.1))
plot(NA, xlim=c(1,ncol(perissodactyla.U)), ylim=c(1,nrow(perissodactyla.U)),
     ylab="", xlab="", axes=FALSE, cex=1.5)

axis(3, at=1:ncol(perissodactyla.U), tick=FALSE, cex.axis=1.1,
     label=parse(text=sprintf("bold(u)[%d]",1:ncol(perissodactyla.U))))

absrng <- max(abs(perissodactyla.U))

for(i in 1:ncol(perissodactyla.U))
  points(
    x = rep(i,nrow(perissodactyla.U)),
    y = 1:nrow(perissodactyla.U),
    cex = 4*abs(perissodactyla.U[,i])/absrng,
    bg = grey(c(0,1))[1.5 + 0.5*sign(perissodactyla.U[,i])],
    pch=22
  )

par(tmp)

perissodactyla.PEM <- PEM(gr1, a = 0)   ## A simpler, single-steepness, model

## perissodactyla.PEM$S2()

evolution.model(
  object = perissodactyla.PEM,
  y = perissodactyla.data[,"log.neonatal.wt"]
) -> opt

opt

perissodactyla.PEM_aux <- PEM(gr1, a = 0)

## Show the usage of model.data
model.data(
  formula = log.neonatal.wt~log.female.wt+Territoriality,
  data = perissodactyla.data
) -> mdat

evolution.model(
  object = perissodactyla.PEM_aux,
  y = mdat$y,
  x = mdat$x
) -> opt_aux

opt_aux

perissodactyla.PEM_aux$S2()

predict(
  object = perissodactyla.PEM,
  newdata = data.frame(
    row.names = sprintf("target_%d",1:6),
    ref = c(1,3,5,2,3,4),
    dist = c(NA,0.1,NA,NA,0.4,1.1),
    lca = c(0,2,2.3,0,2,1.1)
  )
)

predict(
  object = perissodactyla.PEM_aux,
  newdata = data.frame(
    row.names = sprintf("target_%d",1:6),
    ref = c(1,3,5,2,3,4),
    dist = c(NA,0.1,NA,NA,0.4,1.1),
    lca = c(0,2,2.3,0,2,1.1)
  )
)

## A first pemlm model without an auxiliary trait:
pemlm(
  log.neonatal.wt~1,
  perissodactyla.data,
  perissodactyla.PEM
) -> lm1

## Summary of the first model:
summary(lm1)

## Anova of the first model:
anova(lm1)

## A second pemlm model with two auxiliary traits:
pemlm(
  log.neonatal.wt~log.female.wt+Territoriality,
  perissodactyla.data,
  perissodactyla.PEM_aux
) -> lm2

## Summary of the second model:
summary(lm2)

## Anova of the second model:
anova(lm2)

## Adding eigenfunctions to the first model until the smallest AICc is reached
lm1$forward()
lm1

## Calculating the summary and anova 
summary(lm1)
anova(lm1)

## Adding eigenfunctions to the second model until the smallest AICc is reached
lm2$forward()
lm2

summary(lm2)
anova(lm2)

### Making predictions:

sp <- "Tapirus pinchaque"     ## A selected species.
train <- locate(gr1, sp)      ## The locate method.
train

## Initial PEM (single a = 0.5, single psi = 1) built from the residual graph:
pem.train <- PEM(train$x, a = 0)

## Estimate the evolution model:
evolution.model(
  object = pem.train,
  y = mdat$y[names(mdat$y) != sp],
  x = mdat$x[rownames(mdat$x) != sp,]
) -> opt

## Build a linear model:
pemlm(
  formula = log.gestation.length~log.female.wt+log.neonatal.wt+Territoriality,
  data = perissodactyla.data %>% .[rownames(.) != sp,],
  pem = pem.train
) -> lm3

## Summary and anova for the auxiliary trait model:
summary(lm3)
anova(lm3)

## Add the phylogenetic eigenfunctions:
lm3$forward()

## Summary and anova for the final model:
summary(lm3)
anova(lm3)

## Make the prediction:
predict(
  object = lm3,
  newdata = perissodactyla.data[sp,],
  newloc = train$location
) -> prd
prd

## Substitute the missing value for the estimated gestation time:
perissodactyla.data[sp,"log.gestation.length"] <- prd[,"fit"]

## Cross-valisation:

## Obtaining the updated model data:
model.data(
  formula = log.neonatal.wt~log.female.wt+log.gestation.length+Territoriality,
  data = perissodactyla.data
) -> mdat

## Table storing the results:
data.frame(
  observed = mdat$y,
  auxiliary = NA,
  predictions = NA,
  lower = NA,
  upper = NA,
  row.names = names(mdat$y)
) -> perissodactyla.pred

## For each species i:
for(i in 1L:nrow(perissodactyla.pred)) {
  
  ## Calculate the residual graph and location:
  train <- locate(gr1, rownames(perissodactyla.pred)[i])
  
  ## Calculate a PEM:
  pem.train <- PEM(train$x, a = 0)
  
  ## Estimate the evolution model:
  evolution.model(
    object = pem.train,
    y = mdat$y[-i],
    x = mdat$x[-i,]
  ) -> opt
  
  ## Build an empty (auxiliary trait only) pemlm model:
  pemlm(
    formula = log.neonatal.wt~log.female.wt+log.gestation.length+Territoriality,
    data = perissodactyla.data[-i,],
    pem = pem.train
  ) -> lm_cv
  
  ## Make prediction using the empty model:
  predict(
    lm_cv$auxModel,
    perissodactyla.data[i,]
  ) -> perissodactyla.pred[i,2L]
  
  ## Add the PEM eigenfunction(s) on the basis of the AICc:
  lm_cv$forward()
  
  ## Make the prediction using the PEM-based model, including the limits of the
  ## 95% prediction interval:
  predict(
    object = lm_cv,
    newdata = perissodactyla.data[i,],
    newloc = train$location,
    interval = "prediction"
  ) -> perissodactyla.pred[i,3L:5L]
}

## Prediction coefficient using the auxiliary traits alone:
Psquare(perissodactyla.pred$observed, perissodactyla.pred$auxiliary)

## Prediction coefficient using the PEM eigenfunctions:
Psquare(perissodactyla.pred$observed, perissodactyla.pred$predictions)

## Plotting:

## Calculate the range of the whole data (predictions and interval limits):
rng <- range(perissodactyla.pred)

## Save the graphical parameters:
p <- par(no.readonly = TRUE)

## Generates an empty plot:
par(mar=c(5,5,2,2))
plot(NA, xlim=rng, ylim=rng, xlab="Observed (log) neonatal",
     ylab="predicted (log) neonatal", las=1, asp=1)

## Show the predictions without the PEM:
points(x=perissodactyla.pred$observed, y=perissodactyla.pred$auxiliary, pch=21,
       bg="blue")

## Show the predictions and their prediction intervals with the PEM:
arrows(x0=perissodactyla.pred$observed, x1=perissodactyla.pred$observed,
       y0=perissodactyla.pred$lower, y1=perissodactyla.pred$upper, code=3,
       angle=90, length=0.05)
points(x=perissodactyla.pred$observed, y=perissodactyla.pred$predictions,
       pch=21, bg="black")

## 1:1 line:
abline(0,1)

## Restore the graphical parameters:
par(p)

InflMat(gr1) -> res
res

## Update the PEM object, setting a to 0.9:
update(perissodactyla.PEM, a = log(0.9/(1 - 0.9)), psi=NULL)

## Access the parameter values
perissodactyla.PEM$par()

## The parameter values are also available from the edge properties:
edge(perissodactyla.PEM$graph())

## Old trait variance value:
perissodactyla.PEM$S2()

## Recalculating trait variance:
perissodactyla.PEM$S2(y = perissodactyla.data[,"log.neonatal.wt"])

## The new trait variance value is available thereafter:
perissodactyla.PEM$S2()
