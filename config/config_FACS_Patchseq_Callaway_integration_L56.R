#################   Setup   #######################
knitr::opts_chunk$set(echo = TRUE)
#Extra memory
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 4000 * 1024^2)
options(stringsAsFactors=FALSE)


SeuratIntegration <- TRUE
SeuratLabelTransfer <- FALSE 
filtSubclasses <- c("L5 IT", "L6 IT")
filtCallawayProjection <- c("AL","PM")
filtCallawayType <- NULL #c("L56_AL","L56_PM")
filtCallawayMappedTType <- c('L5 IT VISp Hsd11b1 Endou', 'L5 IT VISp Whrn Tox2', 'L5 IT VISp Batf3', 
                             'L5 IT VISp Col6a1 Fezf2','L5 IT VISp Col27a1','L6 IT VISp Penk Col27a1', 
                             'L6 IT VISp Penk Fst','L6 IT VISp Col18a1', 'L6 IT VISp Car3','L6 IT VISp Col23a1 Adamts2')
layers <- "L56" # for filenames
onlyVISp <- TRUE
nGenes  <- 2000
useCuratedGeneList <- TRUE
subSampleFACS <- TRUE
maxSamples <- 100

#################   Directories   #######################

## Folder location for reference feather files (default to FACS)
FACSFolder <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/Mm_VISp_AIT2.3.0_20047_202005/"
PatchseqFolder <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20210818_collapsed40_cpm/"
dendFile <- file.path(FACSFolder, "dend.RData")

callawayFolder <-"/allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data/" #(Location of Callaway data)
callawaySuppFolder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data/Kim2020_supplement/Tables-S4-10/"
callawayDataFile <- "GSE133230_rawCount_945filtered_cells.tsv"
callawayAnnoFile <- "TableS4-Single nuclei sequencing sample metadata.xlsx"
callawayMappingFile <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/EXC_patchseq_paper_2020/Callaway_results_on_v1_and_v1_alm_taxonomies.csv"

outputFolder <- paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/derived_data/ALPM_projections/",layers,"/")
curatedGenesFileName <-  "//allen/programs/celltypes/workgroups/rnaseqanalysis/agatab/exc_mouse_patchseq_R/data/L56_gene_list_for_integration.csv"

outputDataFileName <- paste0("input_data_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.RData")
outputIntegratedDataFileName    <- paste0("SeuratIntegration_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.RData")
outputCSVFileName <- paste0("integrated_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv")
predictionsCallawayCSVFileName <- paste0("Seurat_predictions_",layers,"_Callaway_",layers,"IT_VISp_FACS_Patchseq.csv")
predictionsPatchseqCSVFileName <- paste0("Seurat_predictions_",layers,"_Patchseq_",layers,"IT_VISp_FACS_Patchseq.csv")