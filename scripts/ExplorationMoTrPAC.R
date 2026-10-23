# This script aims to explore metabolomic in humans

#libraries----
library(dplyr) #to use pipe
library(ggplot2) #to graph
library(readr) #to read csv
library(tidyr)  # to use drop-na()
library(ggpubr) 
library(purrr)
library(Hmisc)
library(lme4)
library(lmerTest)
library(emmeans)
library(pracma)
library(lubridate)
library(stringr)
library(forcats)
library(patchwork)
library(ggpattern)
library(car)
library(broom) 
library(rstatix)
library(ARTool)
library(ggrepel)


my_data <- read_tsv("Documents/GitHub/data/data/motrpac-data-download-132a11fcaaf79154498d3d0974728863/analysis/human-precovid-sed-adu/c2.0/transcriptomics/metadata/human-precovid-sed-adu_t11-adipose_transcript-rna-seq_metadata_features_v2.0.txt")
my_data2 <- read_tsv("Documents/GitHub/data/data/motrpac-data-download-132a11fcaaf79154498d3d0974728863/132a11fcaaf79154498d3d0974728863.list.manifest.json")
