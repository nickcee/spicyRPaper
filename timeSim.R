# Load packages
library(spicyR)
library(spatstat)
library(BiocParallel)
library(tidyverse)
library(plotROC)




## INITIALISE

seed = 51773
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))
nPatients <- 40
nIm <- 1
nSim <- 500
nCores <- 50
counts <- 10^seq(1,5,0.25)


time1 <- time2 <- NULL


for(count in counts){
  
  
  
  
  sCount1 <- count/2
  sCount2 <- count/2
  a <- rpoispp(sCount1/1000^2, win = window)
  b <- rpoispp(sCount1/1000^2, win = window)
  
  x <- c(a$x, b$x)
  y <- c(a$y, b$y)
  
  cellType <- c(rep("A", a$n), rep("B", b$n))
  imageID <- c(rep(paste(1,1,sep = "_"), a$n+b$n))
  
  imageID <- factor(imageID)
  
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  
  cellExp <- spicyR::SegmentedCells(cellExp, verbose = FALSE)
  
  s1 <- Sys.time()
  assoc <- getPairwise(cellExp)
  s2 <- Sys.time()
  cat(s2-s1, "\n")
  
  time1[as.character(count)] <- s1
  time2[as.character(count)] <- s2
  
}


pdf("Output/timeSim.pdf", height = 4, width = 5)
df <- data.frame(time = time2-time1, cellNumber = counts)
ggplot(df, aes(cellNumber, time)) + geom_point()  + scale_x_log10() + scale_y_log10() + theme_classic() + labs(y = "Seconds", x= "Number of cells in image")
dev.off()


