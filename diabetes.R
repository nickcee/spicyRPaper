

library(inline)
openblas.set.num.threads <- cfunction( signature(ipt="integer"),
                                       body = 'openblas_set_num_threads(*ipt);',
                                       otherdefs = c ('extern void openblas_set_num_threads(int);'),
                                       libargs = c ('/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so'),
                                       language = "C",
                                       convention = ".C"
)
openblas.set.num.threads(1)




load("/dskh/nobackup/biostat/datasets/spatial/IMC_Diabetes_Damond2019/cellExp.RData")



ph <- imagePheno(cellExp)
ph$case <- factor(ph$case)

imagePheno(cellExp) <- ph



t1 <- Sys.time()


spicyTest <- spicy(cellExp, 
                   condition = "stage", 
                   subject = "case",
                   BPPARAM = BiocParallel::MulticoreParam(50)
                  )

t2 <- Sys.time()

t2 - t1


top <- topPairs(spicyTest, n = 20, coef = 2)
top

signifPlot(spicyTest)



case <- imagePheno(cellExp)$case
names(case) <- imagePheno(cellExp)$imageID

tab <- table( case[imageID(cellExp)], cellType(cellExp))

tab <- table( imageID(cellExp), cellType(cellExp))


boxplot( spicyTest$pairwiseAssoc$beta__Tc ~ imagePheno(cellExp)$case + imagePheno(cellExp)$stage, las = 2)
boxplot( tab[imagePheno(cellExp)$imageID,"beta"] ~ imagePheno(cellExp)$case + imagePheno(cellExp)$stage, las = 2)

boxplot( tab[imagePheno(cellExp)$imageID,"Tc"] ~ imagePheno(cellExp)$case + imagePheno(cellExp)$stage, las = 2)






df <- data.frame(beta__Tc = spicyTest$pairwiseAssoc$beta__Tc, beta = tab[imagePheno(cellExp)$imageID,"beta"], Tc = tab[imagePheno(cellExp)$imageID,"Tc"], case = imagePheno(cellExp)$case, stage = imagePheno(cellExp)$stage)

df[is.na(df)] <- -250

ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, beta__Tc, colour = stage)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 


p1 <- ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(log10(beta+1), beta__Tc, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 
p2 <- ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")), aes(log10(Tc + 1), beta__Tc, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm")

library(patchwork)

p1 + p2





cellExp2 <- cellExp[imagePheno(cellExp)$stage%in%c("Non-diabetic", "Onset"),]

ph <- droplevels(imagePheno(cellExp2))
imagePheno(cellExp2) <- ph

cS <- cellSummary(cellExp2)
cS <- cS[!cS$cellType%in%c("unknown"),]
cS$cellType <- droplevels(cS$cellType)
cellExp2 <- SegmentedCells(cS)
imagePheno(cellExp2) <- ph


case <- imagePheno(cellExp2)$case
names(case) <- imagePheno(cellExp2)$imageID
tab <- table( case[imageID(cellExp2)], cellType(cellExp2))


set.seed(51773)
use <- 1:nrow(cellExp2)#sample(1:nrow(cellExp2),50)

t1 <- Sys.time()


spicyTest2 <- spicy(cellExp2[use,], 
                   condition = "stage", 
                   subject = "case",
                   BPPARAM = BiocParallel::MulticoreParam(50),
                   weights = TRUE,
                   #sigma = 20,
                   Rs = seq(10,100,10),#c(20,50,100),
                   weightsByPair = TRUE,
                   weightFactor = 1,
                   fast = TRUE,
                   includeZeroCells = TRUE
)

t2 <- Sys.time()

t2 - t1


top2 <- topPairs(spicyTest2, n = 20)
top2


signifPlot(spicyTest2, breaks = c(-2,2,0.1))


pval <- spicyTest2$p.value[grep("beta",rownames(spicyTest2$p.value)),]
pval[order(pval[,2]),]

tab <- table( imageID(cellExp2), cellType(cellExp2))


cell1 <- "naiveTc"
cell2 <- "macrophage"
lab <- paste(cell1,cell2, sep = "__")

df <- data.frame(assoc = spicyTest2$pairwiseAssoc[[lab]], cell1 = tab[imagePheno(cellExp2)$imageID, cell1], cell2 = tab[imagePheno(cellExp2)$imageID, cell2], imageID = imagePheno(cellExp2)$imageID, case = imagePheno(cellExp2)$case, stage = imagePheno(cellExp2)$stage, weights = spicyTest2$weights[[lab]])

#df[is.na(df)] <- -250

ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, assoc, colour = stage, alpha = weights)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 


ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, assoc, colour = stage, alpha = weights)) + geom_boxplot() 

p1 <- ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(log10(cell1 + 1), assoc, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 
p2 <- ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")), aes(log10(cell2 + 1), assoc, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm")

library(patchwork)

p1 + p2


ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, cell1, colour = stage, alpha = weights)) + geom_point() +geom_boxplot()

ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, cell2, colour = stage, alpha = weights)) + geom_point() +geom_boxplot()




d <- cellSummary(cellExp2["M18",]) %>% as.data.frame()
ggplot(d,aes(x,y)) + geom_point(colour = "grey60") + geom_point(data = filter(d, cellType %in% c("alpha", "beta", "delta", "gamma", "macrophage")), aes(colour = cellType))



tab <- table( imageID(cellExp2), cellType(cellExp2[use,]))


cell1 <- "naiveTc"
cell2 <- "beta"
lab <- paste(cell1,cell2, sep = "__")

df <- data.frame(assoc = spicyTest2$pairwiseAssoc[[lab]], cell1 = tab[imagePheno(cellExp2[use,])$imageID, cell1], cell2 = tab[imagePheno(cellExp2[use,])$imageID, cell2], case = imagePheno(cellExp2[use,])$case, stage = imagePheno(cellExp2[use,])$stage, weights = spicyTest2$weights[[lab]])

#df[is.na(df)] <- -250

ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, assoc, colour = stage, alpha = weights)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 


ggplot(dplyr::filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, assoc, colour = stage, alpha = weights)) + geom_boxplot() 












df <- data.frame(Th__B = spicyTest2$pairwiseAssoc$Th__B, Th = tab[imagePheno(cellExp2)$imageID,"Th"], B = tab[imagePheno(cellExp2)$imageID,"B"], case = imagePheno(cellExp2)$case, stage = imagePheno(cellExp2)$stage)

#df[is.na(df)] <- -250

ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(case, Th__B, colour = stage)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 


p1 <- ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")),  aes(log10(B+1), Th__B, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm") 
p2 <- ggplot(filter(df, stage %in% c("Non-diabetic", "Onset")), aes(log10(Th + 1), Th__B, colour = stage, shape = case)) + geom_point() + geom_smooth(se = FALSE, method = "lm")

library(patchwork)

p1 + p2





t1 <- Sys.time()


spicyTest <- spicy(cellExp, 
                   condition = "stage", 
                   subject = "case",
                   BPPARAM = BiocParallel::MulticoreParam(20),
                   weights = TRUE,
                   #sigma = 20,
                   Rs = c(20, 50, 100),
                   weightsByPair = TRUE
)

t2 <- Sys.time()

t2 - t1










library(imcdatasets)
sce <- DamondPancreas2019Data(data_type = "sce")

cellExpDam <- SegmentedCells(colData(sce), imageIDString = "ImageName", spatialCoords = c("Pos_X", "Pos_Y"), cellTypeString = "CellType")

ph <- unique(colData(sce)[, c("ImageName", "case", "stage", "slide")])
ph$imageID <- ph$ImageName
rownames(ph) <- ph$imageID
ph$case <- as.factor(ph$case)
ph$stage <- factor(ph$stage, levels = c("Non-diabetic", "Onset", "Long-duration"))

imagePheno(cellExpDam) <- ph



t1 <- Sys.time()

spicyTestDam <- spicy(cellExpDam, 
                   condition = "stage", 
                   subject = "case",
                   BPPARAM = BiocParallel::MulticoreParam(50),
                   weights = FALSE,
                   sigma = 20,
                   Rs = c(20,50,100),
                   weightFactor = 1
                  )

t2 <- Sys.time()

t2 - t1


top <- topPairs(spicyTestDam, n = 20)
top


signifPlot(spicyTestDam)


library(pheatmap)
fdr = FALSE
breaks = c(-3, 3, 0.1)
colors = c("blue", "white", "red")
marksToPlot = NULL
results = spicyTestDam
pVal <- results$p.value[, 2]
marks <- unique(results$comparisons$from)
if (is.null(marksToPlot)) 
  marksToPlot <- marks
if (min(pVal, na.rm = TRUE) == 0) {
  pVal[pVal == 0] <- pVal[pVal == 0] + 10^floor(log10(min(pVal[pVal > 
                                                                 0], na.rm = TRUE)))
}
if (fdr) {
  pVal <- p.adjust(pVal, method = "fdr")
}
isGreater <- which(results$coefficient[, 2] > 0)
pVal <- log10(pVal)
pVal[isGreater] <- abs(pVal[isGreater])
pVal <- matrix(pVal, nrow = length(marks))
colnames(pVal) <- marks
rownames(pVal) <- marks
breaks <- seq(from = breaks[1], to = breaks[2], by = breaks[3])
pal <- colorRampPalette(colors)(length(breaks))
heatmap <- pheatmap(pVal[marksToPlot, marksToPlot], color = pal, 
                    breaks = breaks, cluster_rows = FALSE, cluster_cols = FALSE)
heatmap



w <- spicyTestDam$weights
nCells <- spicyTestDam$nCells
m1 <- unlist(lapply(strsplit(names(w), "__"), function(x)x[1]))
m2 <- unlist(lapply(strsplit(names(w), "__"), function(x)x[2]))
counts1 <- as.vector(nCells[,m1])
counts2 <- as.vector(nCells[,m2])
df <- data.frame(counts1, counts2, weights = unlist(w))
ggplot(df, aes(sqrt(counts1),sqrt(counts2), colour = (weights)^(1/1))) + geom_point() + theme_classic() + scale_color_continuous(type = "viridis")#+ scale_colour_gradient()




ggplot(df, aes(sqrt(counts1),sqrt(counts2), colour = weights)) + geom_point() + theme_classic() + scale_color_continuous(type = "viridis")#+ scale_colour_gradient()




assoc <- spicyTest2$pairwiseAssoc[["gamma__macrophage"]]





w <- spicyTest2$weights
nCells <- spicyTest2$nCells
m1 <- unlist(lapply(strsplit(names(w), "__"), function(x)x[1]))
m2 <- unlist(lapply(strsplit(names(w), "__"), function(x)x[2]))
counts1 <- as.vector(nCells[,m1])
counts2 <- as.vector(nCells[,m2])
df <- data.frame(counts1, counts2, weights = unlist(w))
ggplot(df, aes(sqrt(counts1),sqrt(counts2), colour = log10(weights)^(1/1))) + geom_point() + 
  theme_classic() + scale_color_continuous(type = "viridis") +
  xlim(0,65) + ylim(0,65)+
  labs(x = "sqrt ( cell type A )", y = "sqrt ( cell type B )", colour = "log10 ( weights )")





df <- data.frame(count1ToWeight, count2ToWeight, resSqToWeight)
ggplot(df, aes(sqrt(count1ToWeight),sqrt(count2ToWeight), colour = log10(resSqToWeight+1))) + geom_point() + 
  theme_classic() + scale_color_continuous(type = "viridis") +
  xlim(0,65) + ylim(0,65)+
  labs(x = "sqrt ( cell type A )", y = "sqrt ( cell type B )", colour = "log10 ( weights )")




log10(resSqToWeight+1)~ s(log10(count1ToWeight+1),bs="mpd")+s(log10(count2ToWeight+1)

