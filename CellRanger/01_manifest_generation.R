## In case you have excel sheet and want to safe it as txt file

library(tidyverse)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Fastq manifest
samplesInfo <- readxl::read_xlsx("")

## Cellranger output manifest
samplesInfo  %>%
  write.table("samples_manifest.txt", row.names = F, col.names = F, quote = F, sep = "\t")
