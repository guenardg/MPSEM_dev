
locateAlt <- function(x, target, ...) {
  
  if(is.character(target)) {
    
    idx <- match(target, rownames(x))
    if(any(is.na(idx)))
      stop("Unknown target(s): ", paste(target[is.na(idx)], collapse=", "))
    
    data.frame(
      row.names = target,
      ref = idx,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
    
    target <- idx
    
  } else {
    
    outrange <- which((target < 0) | (target > nrow(x)))
    if(length(outrange))
      stop("Unknown target(s): ", paste(target[outrange], collapse=", "))
    
    data.frame(
      row.names = names(target),
      ref = target,
      dist = as.double(rep(NA, length(target))),
      ladist = double(length(target))
    ) -> ttab
  }
  
  if(is.null(x$species))
    stop("'x' has no vertex property called 'species'")
  
  if(!all(x$species[target]))
    stop("Non-species Target(s): ",
         paste(target[!x$species[target]], collapse=", "))
  
  ## Sanity check: is there any terminal vertex not marked as species?
  if(!all(x$species[getTerminal(x)]))
    stop("Sanity check failed: the graph has one or more terminal vertices ",
         "not marked as species.\nFunction purge.terminal() can be used to ",
         "discard them.")
  
  ## Sanity check: is there any median vertex not marked as a species?
  if(!all(x$species[getMedian(x)]))
    stop("Sanity check failed: the graph has one or more median vertices not ",
         "marked as species.\nFunction purge.median() can be used to discard ",
         "them automatically.")
  
  edge <- edge(x)
  
  if(is.null(edge$distance))
    stop("'x' has no edge property called 'distance'")
  
  vrm <- logical(nrow(x))
  erm <- logical(nrow(edge))
  trm <- logical(length(target))
  
  ## Step 1: removing any target that is terminal vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=1L
    for(i in which(!trm))
      if(!any(!erm & (edge[[1L]] == target[i]))) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          erm[up] <- TRUE
          vrm[target[i]] <- TRUE
          trm[i] <- TRUE
          
          ## Problem here...
          ## ttab$ref[i] <- edge[up,1L]
          ## ttab$ladist[i] <- edge$distance[up]
          ##
          ## Should look more like it:
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          ##
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 2: removing any target that is a median vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=1L
    for(i in which(!trm)) {
      
      down <- which(!erm & (edge[[1L]] == target[i]))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == target[i]))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[target[i]] <- TRUE
            
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            
            erm[down] <- TRUE
            trm[i] <- TRUE
            
            s <- which(ttab$ref == edge[down,1L])
            ttab$ref[s] <- up
            ttab$dist[s] <- ifelse(is.na(ttab$dist[s]), dup, ttab$dist[s] + dup)
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Step 3: Purging any new terminal vertex (ie., a previously non-terminal
  ## vertex that have become a terminal vertex following the removal of terminal
  ## target species) that is not marked as a species.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=10L
    for(i in which(!(vrm | x$species)))
      if(!any(!erm & (edge[[1L]] == i))) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          s <- which(is.na(ttab$dist) & ttab$ref == edge[up,2L])
          ttab$ref[s] <- edge[up,1L]
          ttab$ladist[s] <- ttab$ladist[s] + edge$distance[up]
          
          erm[up] <- TRUE
          vrm[edge[up,2L]] <- TRUE
          
          end <- FALSE
          break
        }
      }
  }
  
  ## Step 4: purging any new median vertex (ie., a previously non-median vertex
  ## that have become a median vertex following the previous removal of terminal
  ## vertice; either target species at step 1 or vertices not marked as a
  ## species at step 2-3) that is not-marked as a species vertex.
  end <- FALSE
  while(!end) {
    end <- TRUE
    
    ## i=3L
    for(i in which(!(vrm | x$species))) {
      
      down <- which(!erm & (edge[[1L]] == i))
      
      if(length(down) == 1L) {
        
        up <- which(!erm & (edge[[2L]] == i))
        
        if(length(up) == 1L) {
          
          if(!any(!erm & (edge[,1L] == edge[up,1L]) &
                  (edge[,2L] == edge[down,2L]))) {
            
            vrm[i] <- TRUE
            edge[up,2L] <- edge[down,2L]
            dup <- edge$distance[up]
            edge$distance[up] <- sum(edge$distance[c(up,down)])
            erm[down] <- TRUE
            
            ## Treat vertices differently from edges
            ## Vertices:
            s <- which(is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$dist[s] <- dup
            
            ## Edges (if any):
            s <- which(!is.na(ttab$dist) & (ttab$ref == edge[down,1L]))
            ttab$ref[s] <- up
            ttab$ladist[s] <- ttab$ladist[s] + ttab$dist[s]
            ttab$dist[s] <- dup
            
            end <- FALSE
            break
          }
        }
      }
    }
  }
  
  ## Vertices that cannot be removed simply lose their species status.
  x$species[target[!trm]] <- FALSE
  
  ## Recalculating the new edge indices:
  mask <- rep(NA, nrow(edge))
  mask[!erm] <- 1L:sum(!erm)
  
  ## Changing the edge indices from the target table:
  ttab$ref[!is.na(ttab$dist)] <- mask[ttab$ref[!is.na(ttab$dist)]]
  
  ## Removing the edges marked for removal:
  edge <- edge[!erm,]
  
  ## Recalculating the new vertex indices:
  mask <- rep(NA, nrow(x))
  mask[!vrm] <- 1L:sum(!vrm)
  
  ## Changing the vertex indices from the target table:
  ttab$ref[is.na(ttab$dist)] <- mask[ttab$ref[is.na(ttab$dist)]]
  
  ## Reindexing the vertices:
  edge[[1L]] <- mask[edge[[1L]]]
  edge[[2L]] <- mask[edge[[2L]]]
  
  ## Removing the vertices marked for removal and reassign the edges:
  out <- x[!vrm,,drop=FALSE]
  edge(out) <- edge
  class(out) <- c("graph",class(out)[-which(class(out) == "graph")])
  
  if(!is.null(attr(out,"processOrder")))
    attr(out,"processOrder") <- MPSEM:::getProcessOrder(out)
  
  if(!is.null(attr(out,"dist")))
    attr(out,"dist") <- graphDist(out)
  
  attr(out,"removedVertex") <- which(vrm)
  attr(out,"removedEdge") <- which(erm)
  
  list(
    x = out,
    location = ttab
  )
}
