## INITIALISE
s1 <- Sys.time()

library(spicyR)
library(spatstat)

set.seed(51773)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))
nPatients <- 20
nIm <- 10
nSim <- 10000
nsimBoot <- 1000
nCores <- 75
counts <- seq(from = 50, to = 500, by = 10)
adjP <- seq(from = 0.75, to = 1.5, length.out = nPatients)

## SIGNAL
s1 <- Sys.time()
set.seed(51773)


library(parallel)

sim  <- function(i){
    print(i)
    x <- c()
    y <- c()
    cellType <- c()
    imageID <- c()
    
    for (p in 1:nPatients) {
        for (j in 1:nIm) {
            sCount <- sample(counts,1)
            a <- rpoispp(sCount/1000^2, win = window)
            aDens <- density(a, adjust = adjP[p])
            aDens <- aDens#/sum(aDens)
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
    cellExp <- SegmentedCells(cellExp)
    
    phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    imagePheno(cellExp) <- phenoData
    
    
    
    test.MEM.Weights <- spicy(cellExp,
                              condition = "condition",
                              subject = "subject",
                              from = "B",
                              to = "A",
                              fast = TRUE, Rs = c(20,50,100))
    
    
    test.MEM.NoWeights <- spicy(cellExp,
                                condition = "condition",
                                subject = "subject",
                                from = "B",
                                to = "A",
                                weights = FALSE,
                                fast = TRUE, Rs = c(20,50,100))
    
    test.MEM.Weights.Boot <- spicy(cellExp,
                                   condition = "condition",
                                   subject = "subject",
                                   from = "B",
                                   to = "A",
                                   nsim = nsimBoot,
                                   fast = TRUE, Rs = c(20,50,100))
    
    
    test.MEM.NoWeights.Boot <- spicy(cellExp,
                                     condition = "condition",
                                     subject = "subject",
                                     from = "B",
                                     to = "A",
                                     weights = FALSE,
                                     nsim = nsimBoot,
                                     fast = TRUE, Rs = c(20,50,100))
    
    
    c(MEM.Weights = test.MEM.Weights$p.value[1,"conditionGroup2"], MEM.NoWeights = test.MEM.NoWeights$p.value[1,"conditionGroup2"],
      MEM.Weights = test.MEM.Weights.Boot$p.value[1,"conditionGroup2"], MEM.NoWeights = test.MEM.NoWeights.Boot$p.value[1,"conditionGroup2"])
}


results <- mclapply(as.list(seq_len(nSim)),sim, mc.cores = nCores)

results <- do.call('rbind', results)

s2 <- Sys.time()
s2-s1


colMeans(results<0.05)
save(results, file = "resultsTW.RData")
resultsTW <- results

## NO SIGNAL
s1 <- Sys.time()
set.seed(51773)
sim  <- function(i){
    print(i)
    x <- c()
    y <- c()
    cellType <- c()
    imageID <- c()
    for (p in 1:nPatients) {
        for (j in 1:nIm) {
            sCount <- sample(counts,1)
            a <- rpoispp(sCount/1000^2, win = window)
            aDens <- density(a, adjust = 1)
            aDens <- aDens#/sum(aDens)
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
    cellExp <- SegmentedCells(cellExp)
    phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    imagePheno(cellExp) <- phenoData
    test.MEM.Weights <- spicy(cellExp,
                              condition = "condition",
                              subject = "subject",
                              from = "B",
                              to = "A",
                              fast = TRUE, Rs = c(20,50,100))
    test.MEM.NoWeights <- spicy(cellExp,
                                condition = "condition",
                                subject = "subject",
                                from = "B",
                                to = "A",
                                weights = FALSE,
                                fast = TRUE, Rs = c(20,50,100))
    test.MEM.Weights.Boot <- spicy(cellExp,
                                   condition = "condition",
                                   subject = "subject",
                                   from = "B",
                                   to = "A",
                                   nsim = nsimBoot,
                                   fast = TRUE, Rs = c(20,50,100))
    test.MEM.NoWeights.Boot <- spicy(cellExp,
                                     condition = "condition",
                                     subject = "subject",
                                     from = "B",
                                     to = "A",
                                     weights = FALSE,
                                     nsim = nsimBoot,
                                     fast = TRUE, Rs = c(20,50,100))
    c(MEM.Weights = test.MEM.Weights$p.value[1,"conditionGroup2"], MEM.NoWeights = test.MEM.NoWeights$p.value[1,"conditionGroup2"],
      MEM.Weights = test.MEM.Weights.Boot$p.value[1,"conditionGroup2"], MEM.NoWeights = test.MEM.NoWeights.Boot$p.value[1,"conditionGroup2"])
}

results <- mclapply(as.list(seq_len(nSim)),sim, mc.cores = nCores)

results <- do.call('rbind', results)

s2 <- Sys.time()
s2-s1

colMeans(results<0.05)
save(results, file = "resultsFW.RData")
resultsFW <- results
