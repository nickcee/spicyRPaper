## Load packages
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
nsimBoot <- 100
nCores <- 40
counts <- seq(from = 20, to = 400, by = 10)
Rs <- seq(from = 10, to = 100, by = 10)

s1 <- Sys.time()


resultsTW = NULL

for(lam in Rs){
  
  lambda = lam  
  
  ## SIGNAL
  
  sim  <- function(i, counts, nPatients, nIm, window, lambda){
    
    set.seed(i)
    
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda + lambda/2)
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
    
    res <- NULL
    for(lamb in Rs){
      test <- spicyR::spicy(cellExp,
                    condition = "condition",
                    from = "B",
                    to = "A",
                    fast = TRUE, Rs = lamb,
                    edgeCorrect = FALSE,
                    weights = FALSE,
                    verbose = FALSE)
      
      res = c(res,test$p.value[1,"conditionGroup2"])
    }
    
    test <- spicyR::spicy(cellExp,
                  condition = "condition",
                  from = "B",
                  to = "A",
                  fast = TRUE, Rs = Rs,
                  edgeCorrect = FALSE,
                  weights = FALSE,
                  verbose = FALSE)
    
    res = c(res,test$p.value[1,"conditionGroup2"])
    
    names(res) <- paste0("Lam",c(Rs,"All"))
    res
    
  }
  
  
  res <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, BPPARAM=MulticoreParam(workers=nCores))
  
  res <- do.call('rbind', res)
  
  res <- cbind(data.frame(res), simulation = "Signal", lambda = lam)
  
  
  resultsTW <- rbind(resultsTW,res)
  
}


s2 <- Sys.time()
s2-s1


## NO SIGNAL

s1 <- Sys.time()

set.seed(seed)
sim  <- function(i, counts, nPatients, nIm, window){
  set.seed(i)
  x <- c()
  y <- c()
  cellType <- c()
  imageID <- c()
  for (p in 1:nPatients) {
    for (j in 1:nIm) {
      sCount1 <- sample(counts,1)
      sCount2 <- sample(counts,1)
      a <- rpoispp(sCount1/1000^2, win = window)
      # aDens <- density(a, adjust = 1000)
      # aDens <- aDens#/sum(aDens)
      b <- rpoispp(sCount2/1000^2, win = window)
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
  
  res <- NULL
  for(lamb in Rs){
    test <- spicyR::spicy(cellExp,
                     condition = "condition",
                     from = "B",
                     to = "A",
                     fast = TRUE, Rs = lamb,
                     edgeCorrect = FALSE,
                     weights = FALSE,
                     verbose = FALSE)
    
    res = c(res,test$p.value[1,"conditionGroup2"])
  }
  
  test <- spicyR::spicy(cellExp,
                   condition = "condition",
                   from = "B",
                   to = "A",
                   fast = TRUE, Rs = Rs,
                   edgeCorrect = FALSE,
                   weights = FALSE,
                   verbose = FALSE)
  
  res = c(res,test$p.value[1,"conditionGroup2"])
  
  names(res) <- paste0("Lam",c(Rs,"All"))
  res
  
}

results <- bplapply(as.list(seq_len(nSim)+seed),sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, BPPARAM=MulticoreParam(workers=nCores))

results <- do.call('rbind', results)

s2 <- Sys.time()
s2-s1


#save(results, file = "resultsFW.RData")
resultsFW <- results

resultsFW <- cbind(data.frame(resultsFW), simulation = "noSignal", lambda = lam)



RESFW = resultsFW

resFW <- NULL
for(i in Rs){
  r <- RESFW
  r$lambda <- i
  resFW <- rbind(r, resFW)
}



results <- rbind(resFW,resultsTW)



save(results, file = "Output/manyRsim.RData")



results <- pivot_longer(results,cols = starts_with("Lam", ignore.case = FALSE), names_to = "method", values_to = "p.value")

results$lambda <- as.factor(results$lambda)
results$method <- factor(results$method, levels = unique(results$method))


g1 <-ggplot(results, aes(m = p.value, d = simulation, colour = method)) + geom_roc(n.cuts = 0, increasing = FALSE) + theme_classic() + facet_wrap(~lambda)
# g1

auc <- calc_auc(g1)
auc$PANEL <- factor(auc$PANEL, levels = levels(auc$PANEL), labels = levels(results$lambda))
auc$group <- factor(auc$group, levels = unique(auc$group), labels = gsub("Lam","",levels(results$method)))

pdf("Output/manyRAvgAUC.pdf", height = 4, width = 8)
ggplot(auc, aes(x = group, y = AUC)) + geom_bar(stat = "identity") + 
  geom_hline(data = filter(auc, group == "All"), aes(yintercept = AUC), colour = "red", linetype = 2) + facet_wrap(~factor(PANEL),nrow = 2) + 
  theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs(x = "Radius")
dev.off()
# 
# df <- auc %>%
#   group_by(PANEL) %>%
#   mutate(rank = rank(AUC))
# 
# ggplot(df, aes(x = group, y = rank)) + geom_bar(stat = "identity") + 
#   geom_hline(data = filter(df, group == "LamAll"), aes(yintercept = rank), colour = "red", linetype = 2) + facet_wrap(~factor(PANEL)) + 
#   theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1))


pdf("Output/manyRAvgRelMaxMin.pdf", height = 5, width = 5)
auc %>%
  group_by(PANEL) %>%
  mutate(relMax = AUC/max(AUC), relMin = AUC/min(AUC))  %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(aveRelMax = mean(relMax), aveRelMin = mean(relMin)) %>%
  ggplot(aes(x = aveRelMax, y = aveRelMin)) + theme_classic() + geom_text(aes(label = gsub("Lam", "", group))) + labs(x = "Average AUC relative to max AUC", y = "Average AUC relative to min AUC")
dev.off()



# 
# results %>%
#   mutate(method = factor(method, levels = unique(method))) %>%
#   group_by(simulation, method, lambda) %>%
#   summarise(pos = mean(p.value < 0.05)) %>%
#   filter(simulation != "noSignal") %>%
#   ggplot(aes(x = method, y = pos)) + geom_bar(stat = "identity") + facet_wrap(~ lambda)+theme(axis.text.x=element_text(angle=90,hjust=1))
# 

# results %>%
#   mutate(method = factor(method, levels = unique(method))) %>%
#   group_by(simulation, method, lambda) %>%
#   summarise(pos = mean(p.value < 0.05)) %>%
#   filter(simulation == "noSignal") %>%
#   ggplot(aes(x = method, y = pos)) + geom_bar(stat = "identity") + facet_wrap(~ lambda)+theme(axis.text.x=element_text(angle=90,hjust=1))
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

