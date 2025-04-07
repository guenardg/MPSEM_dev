
predint <- function(object, lm, location, newx, ..., interval, level = 0.95)
  UseMethod("predint")

predint.PEM2 <- function(object, lm, location, newx, ...,
                         interval = c("confidence","prediction"),
                         level = 0.95) {
  
  interval <- match.arg(interval)
  
  if(!inherits(lm,"lm"))
    stop("Argument 'lm' must be a lm-class object.")
  
  ## if(inherits(lm,"mlm"))
  ##   stop("")
  
  if(!missing(newx)) {
    
    if(!is.matrix(newx))
      newx <- cbind(aux = newx)
    
    if(nrow(location) != nrow(newx))
      stop("The number of target locations (argument 'location': ",
           nrow(location), ") does not correspond to the number of auxiliary ",
           "trait values (argument 'newx': ", nrow(newx), ")")
  }
  
  s <- predict(object, location)#, ...)
  
  rvar <- diag(t(lm$residuals) %*% lm$residuals)/lm$df.residual
  
  Xh <- if(!missing(newx)) cbind(newx, s) else s
  Xh <- cbind(1, Xh[,attr(lm$terms, "term.labels"),drop=FALSE])
  
  pred <- Xh %*% lm$coefficients
  
  if(interval == "none") {
    
    pred
    
  } else {
    
    R <- qr.R(lm$qr)
    invXtX <- solve(t(R) %*% R)
    XhinvXtXtXh <- diag(Xh %*% invXtX %*% t(Xh))
    
    if(interval == "confidence")
      sqrt(
        t(
          attr(s,"vf")/object$nsp +
            matrix(rvar, length(rvar), length(XhinvXtXtXh)) *
            XhinvXtXtXh
        )
      ) -> S
    
    if(interval == "prediction")
      sqrt(
        t(attr(s,"vf") +
            matrix(rvar, length(rvar), length(XhinvXtXtXh)) *
            (1 + XhinvXtXtXh)
        )
      ) -> S
    
    list(
      value = pred,
      lower = pred + S * qt(0.5 * (1 - level), lm$df.residual),
      upper = pred + S * qt(0.5 * (1 - level), lm$df, lower.tail = FALSE)
    )
  }

  ## p <- predict(object=lm, newdata=newdata, se.fit=TRUE, ...)
  ## 
  ## tl <- 0.5*(1 - level)
  ## qt <- qt(p = c(1 - tl,tl), df = p$df, lower.tail = FALSE)
  ## se <- sqrt(p$se.fit^2 + p$residual.scale*attr(s,"vf"))
  ## 
  ## out <- cbind(p$fit + se*qt[1L], p$fit, p$fit + se*qt[2L])
  ## 
  ## c(paste(round(100*tl,1L),"%"),
  ##   "mean",
  ##   paste(round(100*(1 - tl),1L),"%")) -> colnames(out)
  ## 
  ## out
}
