# Load packages
library(spicyR)
library(spatstat)
library(BiocParallel)
library(tidyverse)
library(plotROC)




## INITIALISE
s1 <- Sys.time()


seed = 51773
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))
nPatients <- 40
nIm <- 1
nSim <- 500
nCores <- 50
counts <- seq(from = 20, to = 400, by = 10)
Rs <- seq(from = 10, to = 100, by = 10)

lambda = 40 
multi = c(1,5,10,20)

resultsTW = NULL




for(m in multi){
  
  ## SIGNAL
  s1 <- Sys.time()

  cat("\n",m,"\n")
  
  sim  <- function(i, counts, nPatients, nIm, window, lambda, m){
    
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
        sCount1 <- sample(counts,1)*m
        sCount2 <- sample(counts,1)*m
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
                                  #subject = "subject",
                                  from = "B",
                                  to = "A",
                                  fast = TRUE, Rs = Rs,
                                  edgeCorrect = FALSE,
                                  weights = TRUE,
                                  verbose = FALSE)
    
    test.NoWeights <- spicyR::spicy(cellExp,
                                    condition = "condition",
                                    #subject = "subject",
                                    from = "B",
                                    to = "A",
                                    fast = TRUE, Rs = Rs,
                                    edgeCorrect = FALSE,
                                    weights = FALSE,
                                    verbose = FALSE)
    
    c(Weights = test.Weights$p.value[1,"conditionGroup2"], NoWeights = test.NoWeights$p.value[1,"conditionGroup2"])
    
    
  }
  
  
  res <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, m = m, BPPARAM=MulticoreParam(workers=nCores))
  
  res <- do.call('rbind', res)
  
  s2 <- Sys.time()
  
  cat("\n",s2-s1,"\n")
  
  res <- cbind(data.frame(res), simulation = "Signal", multi = m)
  
  
  resultsTW <- rbind(resultsTW,res)
  
}



resultsFW = NULL

## NO SIGNAL
for(m in multi){
  
  cat("\n",m,"\n")
  s1 <- Sys.time()

  
  sim  <- function(i, counts, nPatients, nIm, window, lambda, m){
    
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
        sCount1 <- sample(counts,1)*m
        sCount2 <- sample(counts,1)*m
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
                          #subject = "subject",
                          from = "B",
                          to = "A",
                          fast = TRUE, Rs = Rs,
                          edgeCorrect = FALSE,
                          weights = TRUE,
                          verbose = FALSE)
    
    test.NoWeights <- spicyR::spicy(cellExp,
                                  condition = "condition",
                                  #subject = "subject",
                                  from = "B",
                                  to = "A",
                                  fast = TRUE, Rs = Rs,
                                  edgeCorrect = FALSE,
                                  weights = FALSE,
                                  verbose = FALSE)
    
    c(Weights = test.Weights$p.value[1,"conditionGroup2"], NoWeights = test.NoWeights$p.value[1,"conditionGroup2"])
    
    
  }
  
  
  res <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, m = m, BPPARAM=MulticoreParam(workers=nCores))
  
  res <- do.call('rbind', res)
  
  s2 <- Sys.time()
  cat("\n",s2-s1,"\n")
  
  res <- cbind(data.frame(res), simulation = "NoSignal", multi = m)
  
  
  resultsFW <- rbind(resultsFW,res)
  
}

# 
# RESFW = resultsFW
# 
# resFW <- NULL
# for(i in Rs){
#   r <- RESFW
#   r$lambda <- i
#   resFW <- rbind(r, resFW)
# }
# 
# 
# 
 results <- rbind(resultsTW,resultsFW)
#results <- results[,-2]


save(results, file = "Output/results.multiSim.RData")



results <- pivot_longer(results,cols = contains("Weights"), names_to = "method", values_to = "p.value")


g1 <-ggplot(results, aes(m = p.value, d = simulation, colour = factor(method))) + geom_roc(n.cuts = 0, increasing = FALSE) + theme_classic() + facet_wrap(~(multi*mean(counts)))

pdf("Output/multiSimROC.pdf", height = 6, width = 6)
g1 + labs(x = "False positive rate", y = "True positive rate", colour = "Method")
dev.off()
# 

auc <- calc_auc(g1)
auc$PANEL <- factor(auc$PANEL, levels = levels(auc$PANEL), labels = unique(results$multi))
auc$group <- factor(auc$group, levels = unique(auc$group), labels = levels(factor(results$method)))
auc$PANEL <- as.numeric(as.character(auc$PANEL))
auc$multi <- auc$PANEL

pdf("Output/multiSimAUC.pdf", height = 5, width = 5)
ggplot(auc, aes(x = multi*mean(counts), y = AUC, colour = factor(group))) + geom_point() + theme_classic() + geom_smooth(method = "loess", span = 2, se = FALSE) + labs(x = "Average number of cells per image", colour = "Method")
dev.off()


# 
# 
# 
# results %>%
#   mutate(method = factor(method, levels = unique(method))) %>%
#   group_by(simulation, method, multi) %>%
#   summarise(pos = mean(p.value < 0.05)) %>%
#  # filter(simulation != "noSignal") %>%
#   ggplot(aes(x = method, y = pos)) + geom_bar(stat = "identity") + facet_wrap(~ multi)+theme(axis.text.x=element_text(angle=90,hjust=1))
# 







#ggplot(results, aes(m = p.value, d = simulation, colour = method)) + geom_roc(n.cuts = 0, increasing = FALSE)

# df <- rbind(cbind(resultsOrig, weights = "orig"), cbind(resultsNew, weights = "new"))
# ggplot(df, aes(m = p.value, d = simulation, colour = method, linetype = weights)) + geom_roc(n.cuts = 0, increasing = FALSE) + theme_classic()
# 
# 
# df %>%
#   group_by(simulation, method, weights) %>%
#   summarise(mean(p.value < 0.05))



#ggplot(mutate(results, p.value = p.value + 0.0000001), aes(m = p.value, d = simulation, colour = method)) + geom_roc( increasing = FALSE) + theme_classic()

