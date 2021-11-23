# Load packages
library(spatstat)
library(BiocParallel)
library(tidyverse)
library(plotROC)
library(spicyR)



## INITIALISE
s1 <- Sys.time()


seed = 123
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))
nPatients <- 40
nIm <- 3
nSim <- 500
nCores <- 50
counts <- seq(from = 20, to = 400, by = 10)
Rs <- seq(from = 10, to = 100, by = 10)

lambda = 40 


  ## SIGNAL
  s1 <- Sys.time()
  
  sim  <- function(i, counts, nPatients, nIm, window, lambda){
    
    set.seed(i)
    
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda + lambda/3)
    adjustSigma = c(g1,g2)+1
    
    x <- c()
    y <- c()
    cellType <- c()
    imageID <- c()
    
    for (p in 1:nPatients) {
      for (j in 1:nIm) {
        sCount1 <- sample(counts,1)
        sCount2 <- sample(counts,1)
        a <- rpoispp(sCount1/1000^2, win = window)
        aDens <- density(a, sigma = adjustSigma[p], kernel = "disc")
        aDens$v <- pmax(aDens$v,0)*sCount2/sCount1
        b <- rpoispp(aDens)
        
        x <- c(x, a$x, b$x)
        y <- c(y, a$y, b$y)
        
        cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
        imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
      }
    }
    
    imageID <- factor(imageID)
    
    cellExp <- data.frame(
      x = x,
      y = y,
      cellType = factor(cellType),
      imageID = imageID
    )
    cellExp <- spicyR::SegmentedCells(cellExp, verbose = FALSE)
    
    phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    spicyR::imagePheno(cellExp) <- phenoData
    
    
    test.Weights <- spicyR::spicy(cellExp,
                                  condition = "condition",
                                  subject = "subject",
                                  from = "B",
                                  to = "A",
                                  fast = TRUE, Rs = Rs,
                                  edgeCorrect = FALSE,
                                  weights = TRUE,
                                  verbose = FALSE)
    
    test.NoWeights <- spicyR::spicy(cellExp,
                                    condition = "condition",
                                    subject = "subject",
                                    from = "B",
                                    to = "A",
                                    fast = TRUE, Rs = Rs,
                                    edgeCorrect = FALSE,
                                    weights = FALSE,
                                    verbose = FALSE)
    
    c(Weights = test.Weights$p.value[1,"conditionGroup2"], NoWeights = test.NoWeights$p.value[1,"conditionGroup2"])
    
    
  }
  
  
  res <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, BPPARAM=MulticoreParam(workers=nCores))
  
  res <- do.call('rbind', res)
  
  s2 <- Sys.time()
  
  cat("\n",s2-s1,"\n")
  
  res <- cbind(data.frame(res), simulation = "Signal")
  
  
  resultsTW <- res
  






## NO SIGNAL

  
  s1 <- Sys.time()
  
  
  sim  <- function(i, counts, nPatients, nIm, window, lambda){
    
    set.seed(i)
    
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda )
    adjustSigma = c(g1,g2)+1
    
    x <- c()
    y <- c()
    cellType <- c()
    imageID <- c()
    
    for (p in 1:nPatients) {
      for (j in 1:nIm) {
        sCount1 <- sample(counts,1)
        sCount2 <- sample(counts,1)
        a <- rpoispp(sCount1/1000^2, win = window)
        aDens <- density(a, sigma = adjustSigma[p], kernel = "disc")
        aDens$v <- pmax(aDens$v,0)*sCount2/sCount1
        b <- rpoispp(aDens)
        
        x <- c(x, a$x, b$x)
        y <- c(y, a$y, b$y)
        
        cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
        imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
      }
    }
    
    imageID <- factor(imageID)
    
    cellExp <- data.frame(
      x = x,
      y = y,
      cellType = factor(cellType),
      imageID = imageID
    )
    cellExp <- spicyR::SegmentedCells(cellExp, verbose = FALSE)
    
    phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    spicyR::imagePheno(cellExp) <- phenoData
    
    
    test.Weights <- spicyR::spicy(cellExp,
                                  condition = "condition",
                                  subject = "subject",
                                  from = "B",
                                  to = "A",
                                  fast = TRUE, Rs = Rs,
                                  edgeCorrect = FALSE,
                                  weights = TRUE,
                                  verbose = FALSE)
    
    test.NoWeights <- spicyR::spicy(cellExp,
                                    condition = "condition",
                                    subject = "subject",
                                    from = "B",
                                    to = "A",
                                    fast = TRUE, Rs = Rs,
                                    edgeCorrect = FALSE,
                                    weights = FALSE,
                                    verbose = FALSE)
    
    c(Weights = test.Weights$p.value[1,"conditionGroup2"], NoWeights = test.NoWeights$p.value[1,"conditionGroup2"])
    
    
  }
  
  
  res <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, BPPARAM=MulticoreParam(workers=nCores))
  
  res <- do.call('rbind', res)
  
  s2 <- Sys.time()
  cat("\n",s2-s1,"\n")
  
  res <- cbind(data.frame(res), simulation = "NoSignal")
  
  
  resultsFW <- res
  



results <- rbind(resultsTW,resultsFW)



save(results, file = "Output/results.baseSim.RData")



results <- pivot_longer(results,cols = contains("Weights"), names_to = "method", values_to = "p.value")


g1 <-ggplot(results, aes(m = p.value, d = simulation, colour = factor(method))) + geom_roc(n.cuts = 0, increasing = FALSE) + theme_classic()

pdf("Output/baseSimROC.pdf", height = 6, width = 6)
g1 + labs(x = "False positive rate", y = "True positive rate", colour = "Method")
dev.off()
# 


# 
# 
# 

g2 <- results %>%
  mutate(method = factor(method, levels = unique(method))) %>%
  group_by(simulation, method) %>%
  summarise(pos = mean(p.value < 0.05)) %>%
 # filter(simulation != "noSignal") %>%
  ggplot(aes(x = method, y = pos, fill = simulation)) + geom_bar(stat = "identity",position = "dodge")+ theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=1)) + geom_hline(yintercept = 0.05, linetype = 2, size = 2)

pdf("Output/baseSimPvalue.pdf", height = 4, width = 4)
g2 + labs(x = "Method", y = "Proportion of p-values < 0.05", fill = "Simulation")
dev.off()






#ggplot(results, aes(m = p.value, d = simulation, colour = method)) + geom_roc(n.cuts = 0, increasing = FALSE)

# df <- rbind(cbind(resultsOrig, weights = "orig"), cbind(resultsNew, weights = "new"))
# ggplot(df, aes(m = p.value, d = simulation, colour = method, linetype = weights)) + geom_roc(n.cuts = 0, increasing = FALSE) + theme_classic()
# 
# 
# df %>%
#   group_by(simulation, method, weights) %>%
#   summarise(mean(p.value < 0.05))



#ggplot(mutate(results, p.value = p.value + 0.0000001), aes(m = p.value, d = simulation, colour = method)) + geom_roc( increasing = FALSE) + theme_classic()


