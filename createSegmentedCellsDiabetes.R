# Load libraries
library(tidyverse)
library(spicyR)

# Data is from https://data.mendeley.com/datasets/cydmwsfztj/2

# Read in datasets.
allCells <- read.csv("All_Cells.csv")
celltype <- read.csv("CellTypes.csv")
meta <- read.csv('Metadata.csv')

# Clean celltype
celltype <- celltype %>%
  mutate(imageID = core, ImageNumber = as.numeric(as.factor(celltype$core)), imageCellID = id, ObjectNumber = as.numeric(lapply(strsplit(id,"_"),function(x)x[2])), cellType = factor(CellType))

cells <- allCells %>%
  mutate(x = AreaShape_Center_X, y = AreaShape_Center_Y) %>%
  select(x,y,ImageNumber, ObjectNumber, starts_with("Intensity_MeanIntensity_CleanStack_"), starts_with("AreaShape_")) %>%
  inner_join(celltype, by = c("ImageNumber", "ObjectNumber"))

cellExp <- SegmentedCells(cells, intensityString = "Intensity_MeanIntensity_CleanStack_", morphologyString = "AreaShape_")

meta$stage <- factor(as.character(meta$stage),levels = c("Non-diabetic", "Onset", "Long-duration"))
meta <- meta %>% mutate(imageID = image) %>% select(-image)

imagePheno(cellExp) <- meta

save(cellExp, file = "cellExpDiabetes.RData")






