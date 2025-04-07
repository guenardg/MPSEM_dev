
lmforAICc <- function(y, x, u, ...) {
  
  df1 <- as.data.frame(u)
  nm <- colnames(df1)
  
  if(missing(x)) {
    AuxTrait <- FALSE
    p1 <- "y~"
  } else {
    AuxTrait <- TRUE
    x <- cbind(aux = x)
    df1 <- data.frame(x, u)
    p1 <- paste("y~", paste(colnames(x), collapse = "+"), sep = "")
  }
  
  included <- numeric(0L)
  candidates <- 1L:length(nm)
  
  repeat {
    
    p2 <- paste(if(length(included)) {
      
      paste(
        if(AuxTrait) "+",
        paste(nm[included], collapse = "+")
      )
    } else
      if(AuxTrait) "" else "1", sep = "")
    
    lm1 <- lm(as.formula(paste(p1, p2, sep = "")), data = df1)
    
    k1 <- length(lm1$coef)
    
    AICc1 <- AIC(lm1) + (2 * k1 * (k1 + 1)/(length(y) - k1 - 1))
    
    AICc2 <- rep(NA, length(candidates))
    
    for(i in 1L:length(candidates)) {
      
      lm2 <- lm(as.formula(paste(p1, p2, " + ", nm[candidates[i]], sep = "")),
                data = df1)
      k2 <- length(lm2$coef)
      
      AICc2[i] <- AIC(lm2) + (2 * k2 * (k2 + 1)/(length(y) - k2 - 1))
    }
    
    if(min(AICc2, na.rm = TRUE) < AICc1) {
      
      i <- candidates[which.min(AICc2)]
      
      cat("Including:",nm[candidates[which.min(AICc2)]],"\n")
      
      included <- c(included, i)
      
      candidates <- candidates[-which.min(AICc2)]
      
    } else {
      
      lm1$AICc <- AICc1
      
      lm1$call <- match.call()
      
      return(lm1)
    }
  }
}

lmforSidak <- function(y, x, u, ..., alpha = 0.05) {
  
  df1 <- as.data.frame(u)
  nm <- colnames(df1)
  
  if(missing(x)) {
    AuxTrait <- FALSE
    p1 <- "y~"
  } else {
    AuxTrait <- TRUE
    x <- cbind(aux = x)
    df1 <- data.frame(x, u)
    p1 <- paste("y~", paste(colnames(x), collapse = "+"), sep = "")
  }
  
  included <- numeric(0L)
  candidates <- 1L:length(nm)
  
  repeat {
    
    p2 <- paste(if(length(included)) {
      paste(
        if(AuxTrait) "+",
        paste(nm[included], collapse = "+")
      )
    } else
      if(AuxTrait) "" else "1", sep = "")
    
    lm1 <- lm(as.formula(paste(p1, p2, sep = "")), data = df1)
    
    pval <- rep(NA, length(candidates))
    
    for(i in 1L:length(candidates)) {
      
      lm2 <- lm(as.formula(paste(p1, p2, " + ", nm[candidates[i]], sep = "")),
                data = df1)
      
      aovcomp <- anova(lm1, lm2)
      
      pow <- length(candidates) - length(included)
      
      pval[i] <- 1 - (1 - aovcomp[["Pr(>F)"]][2L])^pow
      
    }
    
    if(min(pval, na.rm = TRUE) < alpha) {
      
      i <- candidates[which.min(pval)]
      
      cat("Including:",nm[candidates[which.min(pval)]],"\n")
      
      included <- c(included, i)
      
      candidates <- candidates[-which.min(pval)]
      
    } else {
      
      aovlm1 <- anova(lm1)
      lm1[["Familiwise"]] <- aovlm1[["Pr(>F)"]]
      nc <- ncol(u$pem()$u)
      
      match(
        sprintf("U_%d", (1L:nc)[included]),
        rownames(aovlm1)
      ) -> idx
      
      pow <- nc:(nc - length(included) + 1L)
      lm1[["Familiwise"]][idx] <- 1 - (1 - aovlm1[["Pr(>F)"]][idx])^pow
      
      lm1$call <- match.call()
      
      return(lm1)
    }
  }
}
