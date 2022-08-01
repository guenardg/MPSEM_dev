##
### Development script.
##
## rm(list=ls())
compile <- function() {
  try(dyn.unload("../MPSEM/src/MPSEM.so"),silent=TRUE)
  system("R CMD SHLIB -o ../MPSEM/src/MPSEM.so ../MPSEM/src/*.c")
  dyn.load("../MPSEM/src/MPSEM.so")
  for(i in list.files("../MPSEM/R","*.R"))
    source(file.path("../MPSEM/R",i))
}
## compile()
library(MPSEM)
##
