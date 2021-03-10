data <- load("data.rda")

set.seed(57713)
spicyTest <- spicy(data, 
                   condition = "stage", 
                   subject = "case",
                   weights = T,
                   nsim = 20000,
                   fast = F)

signifPlot(spicyTest, 
           breaks=c(-3, 3, 0.5), 
           fdr = F, 
           marksToPlot = c("alpha", "beta", "gamma", "delta", "B", 
                           "naiveTc", "Th", "Tc", "neutrophil", "macrophage"))
