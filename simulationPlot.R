library(ggplot2)
library(pROC)

load("resultsFW.Rdata")
resultsFW <- results
load("resultsTW.Rdata")
resultsTW <- results

results <- rbind(resultsFW, resultsTW)

# ROC Curves
rocobj1 <- roc(rep(c("Null", "Signal"), each = nrow(resultsFW)), results[,1])
rocobj2 <- roc(rep(c("Null", "Signal"), each = nrow(resultsFW)), results[,2])
rocobj3 <- roc(rep(c("Null", "Signal"), each = nrow(resultsFW)), results[,3])
rocobj4 <- roc(rep(c("Null", "Signal"), each = nrow(resultsFW)), results[,4])

df <- data.frame(Specificity = c(rocobj1$specificities,
                                 rocobj2$specificities,
                                 rocobj3$specificities,
                                 rocobj4$specificities),
                 Sensitivity = c(rocobj1$sensitivities,
                                 rocobj2$sensitivities,
                                 rocobj3$sensitivities,
                                 rocobj4$sensitivities),
                 Test = rep(c("Weights", "No Weights", "Weights\n+Boot", "No Weights\n+Boot"),
                            c(length(rocobj1$sensitivities),
                              length(rocobj2$sensitivities),
                              length(rocobj3$sensitivities),
                              length(rocobj4$sensitivities))),
                 Group = rep(c("No Boot", "No Boot", "Boot", "Boot"),
                             c(length(rocobj1$sensitivities),
                               length(rocobj2$sensitivities),
                               length(rocobj3$sensitivities),
                               length(rocobj4$sensitivities))))

df$Test <- factor(df$Test, levels = c("No Weights\n+Boot", "No Weights", "Weights\n+Boot", "Weights"))
df$Group <- factor(df$Group, levels = c("Boot", "No Boot"))

ggplot(data=df, aes(x=Specificity, y=Sensitivity, colour=Test)) +
    geom_line(aes(linetype = Group)) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    theme_classic() +
    scale_colour_manual(values = c("steelblue4", "orange", "firebrick2", "springgreen3")) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_x_reverse() +
    coord_fixed()

## BARPLOT
pValsNull <- colMeans(results2[,1:4]<0.05)*100
pValsSignal <- colMeans(results2[,5:8]<0.05)*100

df <- data.frame(pVals = c(pValsNull, pValsSignal),
                 test = rep(c("Null", "Signal"), each = 4),
                 simulation = rep(c("Weights", "No Weights", "Weights\n+Boot", "No Weights\n+Boot"), times = 2))
df$simulation = factor(df$simulation, levels = c("No Weights", "No Weights\n+Boot", "Weights", "Weights\n+Boot"))


ggplot(df, aes(x=simulation, y=pVals, fill=test)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 5,
             color = "black",
             size = 2,
             alpha = 0.5,
             linetype = 2) +
  theme_classic() +
  labs(x = "",
       y = "% of Significant Simulations",
       fill = "Simulation Type")

