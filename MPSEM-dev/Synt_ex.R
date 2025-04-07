
## Synthetic example
data.frame(
  species = as.logical(c(0,0,1,0,0,0,0,0,0,0,1,1,1)),
  type = c(2,2,3,1,2,2,2,2,2,2,3,3,3),
  x = c(1,3,4,0,1.67,4,1,1.33,2.33,3.33,4.33,4,5),
  y = c(1,1,1,0,0.5,0,-1,0,0,-0.5,-1,-0.5,-0.5),
  row.names = sprintf("V%d",1:13)
) %>%
  st_as_sf(
    coords=c("x","y"),
    crs = NA
  ) %>%
  graph %>%
  add.edge(
    from = c(1,2,1,5,4,4,5,9,4,8,9,4,7,7,6,6,9,10,10),
    to = c(2,3,5,2,1,5,9,2,8,9,6,7,8,9,3,13,10,12,11),
    data = data.frame(
      distance = c(4.2,4.7,3.9,3.0,3.6,2.7,4.4,3.4,3.6,3.3,
                   4.8,3.2,3.5,4.4,2.5,3.4,4.3,3.1,2.2),
      row.names = sprintf("E%d",1:19)
    )
  ) -> gr_ex

## Plot the graph:
plot(gr_ex, cex.min=3, cex.lab=0.6)

## Show the edges of the graph:
edge(gr_ex)

## Identify the edges that are under or at the edges E7 or E17 using a binary
## contrast matrix:
isUnderEdge(gr_ex, c("E7","E17"))

## Identify the edges that are under vertices V5 or V9 using a binary
## contrast matrix:
isUnderVertex(gr_ex, c("V5","V9"))
