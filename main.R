install.packages('pcr')
# Needs data in table format
devtools::install_github("ewallace/tidyqpcr", build_vignettes = TRUE)
BiocManager::install("ReadqPCR")
devtools::install_github("jimrperkins/NormqPCR")


BiocManager::install("EasyqpcR") # Corrupted dependencies

source("https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R")



library(magrittr)

### plateDesignFile needs 3 columns: well, group and sample names
### complexGeneDesignFilePresent - chodzi o to, że na jednej płytce masz różne geny. wtedy musisz załączyć do każdego pliku z raw data plik z projektem genów na płytce. Jeśli F, to po prostu nazwa genu dla tej płytki. Jeśli T to plik.
### skipRows - used when there is some shit rows before actual colnames. just set it to skip all shit rows, but not colnames
### sampleRepsIn - tells if we chould look for replication information in plateDesign "plate" or in files descriptions "files"
### namesForReplicates - lista nazw dla kolejnych replikatów
### format must be either tsv or excelSheets. In excelSheets format each "file" is separate sheet
### geneCol is not in the data, its produced based on other information
qpcr <- list(
  "opts" = list(
    "sampleRepsIn" = "files",
    "comparisons" = list(c("DOR", "VENT")),
    "outputDir" = "qpcr_output",
    "plateDesign" = list(
      "file" = "plateDesign_gjt_extraction.tsv",
      "wellColname" = "well", ### This needs to have the same values as position column in runs
      "sampleColname" = "sample",
      "groupColname" = "group",
      "RepsColName" = NA),
    "runs" = list(
      "format" = "excelSheets",
      "excelFile" = "sekcja hipp_dane.xlsx",
      "colTypes" = "ccccnccc", ### "ccccncccccnncccccnncccn" for LC96, 
      "wellCol" = "Position",
      "cqCol" = "Cq",
      "skipRows" = 0, ### 0 for LC96
      "complexGeneDesignFilePresent" = F,
      "geneCol" = "gene",
      "emptyName" = NA, ### this cannot be changed, as !is.na is used further, do use this value to change empty names to NA
      "blankName" = "H2O",
      "ladder" = list("x2" = 1250, "x10" = 250, "x50" = 50, "x250" = 10, "x1250" = 2), ### !!! beware of non-ordered ladder names
      "refGenes" = c("gapdh", "hmbs", "tbp", "ywhaz"),
      "namesForReplicates" = list(1, 2, 3) ### !!! co jeśli nazwy replicatów zostały pdoane w plate file? Need to write some checking function
    )
  ))

dir.create(qpcr$opts$outputDir)

qpcr$files = list(
  "lct_1" = list(
    "data" = "LCT", 
    "genesDesign" = "lct", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "lct_2" = list(
    "data" = "LCT (2)", 
    "genesDesign" = "lct", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "lct_3" = list(
    "data" = "LCT (3)", 
    "genesDesign" = "lct", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "ttr_1" = list(
    "data" = "TTR", 
    "genesDesign" = "ttr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "ttr_2" = list(
    "data" = "TTR (2)", 
    "genesDesign" = "ttr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "ttr_3" = list(
    "data" = "TTR (3)", 
    "genesDesign" = "ttr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "trhr_1" = list(
    "data" = "TRHR", 
    "genesDesign" = "trhr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "trhr_2" = list(
    "data" = "TRHR (2)", 
    "genesDesign" = "trhr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "trhr_3" = list(
    "data" = "TRHR (3)", 
    "genesDesign" = "trhr", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "gapdh_1" = list(
    "data" = "GAPDH", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "gapdh_2" = list(
    "data" = "GAPDH (2)", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "gapdh_3" = list(
    "data" = "GAPDH (3)", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "hmbs_1" = list(
    "data" = "HMBS", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "hmbs_2" = list(
    "data" = "HMBS (2)", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "hmbs_3" = list(
    "data" = "HMBS (3)", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "tbp_1" = list(
    "data" = "TBP", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "tbp_2" = list(
    "data" = "TBP (2)", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "tbp_3" = list(
    "data" = "TBP (3)", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "ywhaz_1" = list(
    "data" = "YWHAZ", 
    "genesDesign" = "ywhaz", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "ywhaz_2" = list(
    "data" = "YWHAZ (2)", 
    "genesDesign" = "ywhaz", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "ywhaz_3" = list(
    "data" = "YWHAZ (3)", 
    "genesDesign" = "ywhaz", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]])
)



### Here additional step could be added where we check if columns with given names exist
qpcr$plateDesign <- readr::read_tsv(qpcr$opts$plateDesign$file, col_names = T)



qpcr$input <- purrr::map(
  .x = qpcr$files, 
  .f = function(files){
    
    if (qpcr$opts$runs$format == "tsv") {
      prepareInput(
        rawFileName = files$data, 
        genesDesign = files$genesDesign, 
        repsNameFiles = files$Rep)
    } else if (qpcr$opts$runs$format == "excelSheets") {
      prepareInput(
        rawFileName = qpcr$opts$runs$excelFile, 
        genesDesign = files$genesDesign, 
        repsNameFiles = files$Rep,
        sheetName = files$data)
    }
  })



qpcr$inputMerged <- rlist::list.rbind(qpcr$input)



library(Hmisc)
qpcr$perGene <- purrr::map(
  .x = levels(as.factor(qpcr$inputMerged$gene)), 
  .f = function(gene_){
    prepareDataFromInput(gene_)
  })
names(qpcr$perGene) <- levels(as.factor(qpcr$inputMerged$gene))




### MANUAL STEP - REMOVE CORRUPTED VALUES, BAD LADDER CONCENTRATIONS AND RECALCULATE VALUES ###
### are all samples within ladder? Are blanks reasonably far from last ladder?
qpcr$perGeneQa <- qpcr$perGene

# rows/columns
qpcr$perGeneQa$lct$valuesToNaCoordinates <- list(
  "gene" = list( c(14,6), c(15,6), c(16,6),c(17,7), c(22,6), c(23,6), c(24,6)),
  "ladder" = list( c(5,5), c(5,6), c(5,7) )
)

qpcr$perGeneQa$hmbs$valuesToNaCoordinates <- list(
  "ladder" = list( c(5,5), c(5,6), c(5,7) ))

qpcr$perGeneQa$tbp$valuesToNaCoordinates <- list(
  "ladder" = list( c(5,5), c(5,6), c(5,7) ))

qpcr$perGeneQa$trhr$valuesToNaCoordinates <- list(
  "ladder" = list( c(5,6) ))

qpcr$perGeneQa$ywhaz$valuesToNaCoordinates <- list(
  "ladder" = list( c(5,6) ))

qpcr$perGeneQaPassed <- purrr::map(
  .x = qpcr$perGeneQa, 
  .f = function(gene){
    
    RepairCorruptedValues(gene)
  })

qpcr$perGeneQa <- NULL

### MANUAL STEP - REMOVE CORRUPTED VALUES, BAD LADDER CONCENTRATIONS AND RECALCULATE VALUES ###


# set lowestDilutionToUse and highestDilutionToUse with pmap
qpcr$ladderData <- purrr::map(
  .x = qpcr$perGeneQaPassed, 
  .f = function(gene_){
    calculateLadderStats(gene_$ladder$data)
  })



CheckIfAllReferenceGenesHaveData()



qpcr$refGenesToUse <- purrr::map(
  .x = list("geNorm" = "geNorm", "NormFinder" = "NormFinder" ), 
  .f = function(method){
    GetReferenceValues(
      method, 
      perGene = qpcr$perGeneQaPassed,
      perGeneNames = names(qpcr$perGeneQaPassed),
      refGenesToRemove = "ywhaz")
  })



qpcr$deltaCt <- purrr::map2(
  .x = qpcr$perGeneQaPassed,
  .y = qpcr$ladderData,
  .f = function(gene, ladder){
    CalculateDeltaCt(gene$data, ladder$efficiency)
  })



qpcr$ratio <- rlist::list.clean(
  purrr::map(
    .x = qpcr$deltaCt,
    .f = function(gene){
      CalculateRatio(gene, qpcr$deltaCt$hmbs, "hmbs")
      }) )







qpcr$normalityHomogenity <- purrr::map2(
  .x = qpcr$ratio,
  .y = names(qpcr$ratio),
  .f = function(gene, geneName){
    TidyTestNormalityAndHomogenityOfVarianceForOneVariable(gene, geneName)
  })





qpcr$pairwiseStat <- purrr::map2(
  .x = qpcr$ratio,
  .y = qpcr$normalityHomogenity,
  .f = function(gene, normalityHomogenityStat){
    TidyWilcoxonOrTTestForOneVariable(
      tidyData = gene, 
      tidyTestNormalityAndHomogenityOfVarianceOutput = normalityHomogenityStat)
  })




qpcr$anovas <- purrr::map(
  .x = qpcr$ratio,
  .f = function(gene){
    TidyAnovas(tidyData = gene)
  })





wb <- openxlsx::createWorkbook()
purrr::walk2(
  .x = qpcr$anovas,
  .y = names(qpcr$anovas),
  .f = function(gene, geneName){
    PrintTidyStatisticsToExcel(tidyStatisticSet = gene, variableName = geneName, openxlsxWBookToWriteInto = wb)
  })
openxlsx::saveWorkbook(wb, "test.xlsx", overwrite = TRUE)

qpcr$ratiosBind <- rlist::list.rbind(qpcr$ratio)



TidyDrawGroupBoxplotsForEachVariable(tidyData = qpcr$ratiosBind)



qpcr$correlations <- TidyDrawCorrelationBetweenTwoVariables(tidyData = qpcr$ratiosBind)







###########
### PCA ###

purrr::walk2(
  .x = metadata_analysis$corr, 
  .y = names(metadata_analysis$corr),
  .f = function(dataset, dataset_name){
    
    ggsave(
      filename = paste0(dataset$dir, dataset_name, '_pca_metadata.png'), 
      plot = factoextra::fviz_pca_ind(dataset$pca_scaled, repel = T), 
      dpi = 250)
    
    if (dataset_name != 'no_nas_pca_fil') {
      
      ggsave(
        filename = paste0(dataset$dir, dataset_name, '_pca_metadata_samples.png'), 
        plot = factoextra::fviz_pca_biplot(
          dataset$pca_scaled_samples, 
          repel = T, 
          habillage = dataset$habillage),
        dpi = 250)
    }
  })

### PCA ###
###########











# qpcr <- qpcr_removed_bad_measurements
# 
# 
# qpcr$for_methods$one <- purrr::map2_df(
#   .x = qpcr[["perGeneQaPassed"]],
#   .y = names(qpcr[["perGeneQaPassed"]]),
#   .f = function(analysis, analysis_name){
#     
#     
#     purrr::map2_df(
#       .x = analysis,
#       .y = names(analysis),
#       .f = function(gene, gene_name){
#         
#         data.frame(
#           "ladder_dilutions" = paste0(gene$ladder$data$sample[!is.na(gene$ladder$data$mean)], collapse = ", "),
#           "mean_Cq_of_the_highest_dilution_of_dilution_curve" = max(gene$ladder$data$mean),
#           "mean_Cq_of_NTC" = mean(gene$blanks$data$mean),
#           "analysis_name" = paste0(gene_name, "_", analysis_name))
#       })
#   })
# 
# qpcr$for_methods$ladder <- purrr::map2_df(
#   .x = qpcr$ladderData,
#   .y = names(qpcr$ladderData),
#   .f = function(analysis, analysis_name){
#     
#     
#     purrr::map2_df(
#       .x = analysis,
#       .y = names(analysis),
#       .f = function(gene, gene_name){
#         
#         gene$analysis_name <- paste0(gene_name, "_", analysis_name)
#         
#         return(gene)
#       })
#   })
# 
# qpcr$for_methods$qpcr_qa <- merge(x = qpcr$for_methods$one, y = qpcr$for_methods$ladder, by = "analysis_name")
# 
# openxlsx::write.xlsx(x = qpcr$for_methods$qpcr_qa, file = "qpcr_qa.xlsx")
# 
# qpcr$for_methods$refgenes <- purrr::map2_dfc(
#   .x = qpcr$refGenesToUse,
#   .y = names(qpcr$refGenesToUse),
#   .f = function(analysis, analysis_name){
#     
#     analysis_ <- data.frame(qpcr[["refGenesToUse"]][["h"]][["advOutout"]])
#     
#     colnames(analysis_) <- analysis_name
#     
#     return(analysis_)
#   })
# 
# openxlsx::write.xlsx(x = qpcr$for_methods$refgenes, file = "refgenes.xlsx")

































#only if format == excelSheet
prepareInput <- function(
  rawFileName,
  genesDesign,
  repsNameFiles = NA,
  format = qpcr$opts$runs$format,
  sheetName = NA, 
  sampleRepsIn = qpcr$opts$sampleRepsIn,
  repsNamePlate = qpcr$opts$plateDesign$RepsColName,
  colTypes = qpcr$opts$runs$colTypes, 
  skipRows = qpcr$opts$runs$skipRows, 
  wellCol = qpcr$opts$runs$wellCol,
  cqCol = qpcr$opts$runs$cqCol,
  geneCol = qpcr$opts$runs$geneCol,
  plateDesign = qpcr$plateDesign,
  wellColname = qpcr$opts$plateDesign$wellColname,
  sampleColname = qpcr$opts$plateDesign$sampleColname,
  groupColname = qpcr$opts$plateDesign$groupColname)
{
  library(Hmisc)
  
  if(format == "tsv") {
    
    colNames_ <- readr::read_tsv(rawFileName, col_types = colTypes, col_names = T, skip = skipRows, n_max = 2)
    
    rawFile <- readr::read_tsv(rawFileName, col_types = colTypes, col_names = F, skip = skipRows + 1)
    
    colnames(rawFile) <- colnames(colNames_)
    
  } else if (format == "excelSheets") {
    
    rawFile <- openxlsx::read.xlsx(xlsxFile = rawFileName, sheet = sheetName, colNames = T)
    
    rawFile[[cqCol]] <- as.numeric(rawFile[[cqCol]])
    
  } else if (format != "tsv" && format != "excelSheets") {
    
    stop("Format needs to be either tsv or excelSheets")
  }
  
  rawFile2 <-subset(rawFile, select = c(wellCol, cqCol))

  
  
### !!! Jeśli byśmy chcieli to rozwijać, to warto sprawdzić, czy rawFile nie zawiera kolumny o nazwie geneCol, bo jeśli zawiera to po prostu ją trzeba użyć jako geneDesign
  if (geneCol %in% colnames(rawFile))
  {
    rawFile2[[geneCol]] <- rawFile[[geneCol]]
  } else rawFile2[[geneCol]] <- genesDesign
  
  
  
  if (sampleRepsIn == "files")
  {
    if (!is.na(repsNameFiles))
    {
      rawFile2[["Rep"]] <- repsNameFiles
    } else stop("Reps names for files not provided")
    
  } else if (sampleRepsIn == "plate")
  {
    if (!is.na(repsNamePlate))
    {
      rawFile2[["Rep"]] <- repsNamePlate
    } else { stop("Reps names for plates not provided") }
    
  } else
  {
      stop("qpcr$opts$sampleRepsIn need to be either  plate or flies")
  }
    
  
    
  rawFile2[[sampleColname]] <- recode_values_based_on_key(
    to_recode_chrvec = rawFile2[[wellCol]],
    replace_this_chrvec = plateDesign[[wellColname]],
    with_this_chrvec = plateDesign[[sampleColname]])
  
  rawFile2[[groupColname]] <- recode_values_based_on_key(
    to_recode_chrvec = rawFile2[[wellCol]],
    replace_this_chrvec = plateDesign[[wellColname]],
    with_this_chrvec = plateDesign[[groupColname]])
  
  return(rawFile2)
}








prepareDataFromInput <- function(
    geneName,
    repColName = NA,
    plateDesign = qpcr$plateDesign,
    inputMerged = qpcr$inputMerged,
    sampleColumn = qpcr$opts$plateDesign$sampleColname,
    geneColumn = qpcr$opts$runs$geneCol,
    blankName = qpcr$opts$runs$blankName,
    ladderNames = names(qpcr$opts$runs$ladder),
    cqCol = qpcr$opts$runs$cqCol,
    RepCol = "Rep")
{
  library(Hmisc)
  
  if ( !is.na(repColName)  ) {
    inputMerged <- inputMerged[order(inputMerged[[repColName]]),]
  }
  
  preAnalyzedData <- inputMerged[inputMerged[[geneColumn]] == geneName & !is.na(inputMerged[[sampleColumn]]),]
  
  repLevels <- levels(as.factor(preAnalyzedData[[RepCol]]))
  
  preAnalyzedData <- tidyr::pivot_wider(
    data = preAnalyzedData, names_from = RepCol, values_from = cqCol)
  
  preAnalyzedData <- preAnalyzedData[match(plateDesign[[sampleColumn]], preAnalyzedData[[sampleColumn]]),]
  
  if ( any(plateDesign[[sampleColumn]] != preAnalyzedData[[sampleColumn]]) ) {
    stop("ladderNames != ladderData[[sampleColumn]]")
  }
  
  preAnalyzedData <- AddMeanSdFlagToDF(preAnalyzedData, repLevels)
  
  preAnalyzedData <-preAnalyzedData[order(preAnalyzedData[[sampleColumn]]),]
  
  ladderData <- preAnalyzedData[preAnalyzedData[[sampleColumn]] %in% ladderNames,]
  
  ladderData <- ladderData[match(ladderNames, ladderData[[sampleColumn]]),]
  if ( any(ladderNames != ladderData[[sampleColumn]]) ) {
    stop("ladderNames != ladderData[[sampleColumn]]")
  }
  
  preAnalyzedDataList <- list(
    "data" = preAnalyzedData[preAnalyzedData[[sampleColumn]] %nin%  c(blankName, ladderNames),],
    "ladder" = list(
      "data" = ladderData),
    "blanks" = list(
      "data" = preAnalyzedData[preAnalyzedData[[sampleColumn]] %in% blankName,]))
  
  return(preAnalyzedDataList)
}





AddMeanSdFlagToDF <- function(
  fullDfWithOneColPerSampleReplicate,
  namesForColumnsWithTechnicalReplicates
)
{
  columnsWithTechnicalReplicates <- fullDfWithOneColPerSampleReplicate[, namesForColumnsWithTechnicalReplicates]
  
  fullDfWithOneColPerSampleReplicate$mean <- apply(X = columnsWithTechnicalReplicates, MARGIN = 1, FUN = mean, na.rm=TRUE)
  
  fullDfWithOneColPerSampleReplicate$mean[is.nan(fullDfWithOneColPerSampleReplicate$mean)] <- NA
  
  fullDfWithOneColPerSampleReplicate$sd <- apply(X = columnsWithTechnicalReplicates, MARGIN = 1, FUN = sd, na.rm=TRUE)
  
  fullDfWithOneColPerSampleReplicate$flag <- dplyr::case_when(
    fullDfWithOneColPerSampleReplicate$sd <= 0.5 ~ "OK",
    fullDfWithOneColPerSampleReplicate$sd < 1 & fullDfWithOneColPerSampleReplicate$sd > 0.5 ~ "DANGER",
    fullDfWithOneColPerSampleReplicate$sd >= 1 ~ "ERROR")
  
  return(fullDfWithOneColPerSampleReplicate)
}






calculateLadderStats <- function(
  ladderDataMeans,
  lowestDilutionToUse = 1,
  ladderConcentrations = qpcr$opts$runs$ladder,
  sampleColname = qpcr$opts$plateDesign$sampleColname)
{
  ladderDataMeans <- ladderDataMeans[match(names(ladderConcentrations), ladderDataMeans[[sampleColname]]),]
  
  ladderMeans <- ladderDataMeans[["mean"]]
  
  highestDilutionToUse <- length(ladderMeans[!is.na(ladderMeans)])
  
  fit <- lm(ladderMeans[lowestDilutionToUse:highestDilutionToUse] ~ log10(as.numeric(ladderConcentrations)[lowestDilutionToUse:highestDilutionToUse]) )
  
  slope <- coef(fit)[2]
  
  efficiency <- (10 ^ (-1/slope))
  
  efficiencyPerc <- (efficiency - 1) * 100
  
  r2 <- summary(fit)$r.squared
  
  return_ <- data.frame(
    "slope" = slope, 
    "efficiency" = efficiency, 
    "efficiencyPerc" = efficiencyPerc, 
    "r2" = r2)
  
  rownames(return_) <- "data"
  
  return(return_)
}





GetReferenceValues <- function(
  method,
  perGene,
  perGeneNames,
  refGenesToRemove = NA,
  refGenes = qpcr$opts$runs$refGenes,
  namesForReplicates = qpcr$opts$runs$namesForReplicates,
  groupColname = qpcr$opts$plateDesign$groupColname,
  sampleColname = qpcr$opts$plateDesign$sampleColname
)
{
  if ( !is.na(refGenesToRemove) ) {
    library(Hmisc)
    refGenes <- refGenes[refGenes %nin% refGenesToRemove]
  }
  
  
  groupNamesTemp <- perGene[[1]][["data"]]
  
  groupNames <- groupNamesTemp[order(groupNamesTemp[[sampleColname]]),]
  
  groupNames <- groupNames[groupColname]
  
  
  refDf <- purrr::pmap_dfc(
    .l = list(perGene, perGeneNames), 
    .f = function(gene, geneName)
    {
      
      if (tolower(geneName) %in% tolower(refGenes))
      {
        gene_data <- gene[["data"]]
        gene_data <- gene_data[order(gene_data[[sampleColname]]),]
        gene_data <- gene_data["mean"]
        colnames(gene_data) <- c(geneName)
        return(gene_data)
      }
    })
  
  normRefDf <- NormqPCR::selectHKs(as.matrix(refDf), Symbols = colnames(refDf), group = groupNames[[groupColname]], method = method, minNrHKs = 2, log = TRUE, trace = TRUE, na.rm = TRUE)
  
  advOutout <- capture.output(NormqPCR::selectHKs(as.matrix(refDf), Symbols = colnames(refDf), group = groupNames[[groupColname]], method = method, minNrHKs = 2, log = TRUE, trace = TRUE, na.rm = TRUE), type = c("output", "message"))
  
  print(paste0("FINISHED ANALYSIS USING ", method))
  
  return(list("output" = normRefDf, "advOutout" = advOutout))
}






CheckIfAllReferenceGenesHaveData <- function(
  perGeneDatasets = qpcr$perGene,
  refGenesNames = qpcr$opts$runs$refGenes
)
{
  for (gene in refGenesNames)
  {
    if (gene %in% names(perGeneDatasets)) {
      print(paste0(gene, " is present in datasets"))
    } else stop(paste0(gene, " IS ABSENT IN DATASETS. BEWARE OF CAPITALIZED LETTERS."))
  }
}







CalculateDeltaCt <- function(
  geneData,
  eneEfficiency,
  meanFromAllSamplesOrControlGroup = "all", ### "all" or "control"
  controlGroupName = NA, ### only for meanFromAllSamplesOrControlGroup = "control"
  groupName = qpcr$opts$plateDesign$groupColname,
  cqCol = qpcr$opts$runs$cqCol,
  deltaCtColName = "deltaCt"
)
{
  if (meanFromAllSamplesOrControlGroup == "all")
  {
    meanOfMeans <- mean(geneData[["mean"]], na.rm = TRUE)
    
  } else if (meanFromAllSamplesOrControlGroup == "control")
  {
    meanOfMeans <- mean( 
      subset(
        x = geneData,
        subset = geneData[[groupName]] == controlGroupName,
        select = groupName ),
      na.rm = TRUE)
    
    assertthat::not_empty(meanOfMeans)
    
  } else { stop("meanFromAllSamplesOrControlGroup need to be either all or control") }
  
  geneData[[deltaCtColName]] <- NA
  
  for (sampleNb in seq_along(geneData[["mean"]]))
  {
    geneData[[deltaCtColName]][[sampleNb]] <- eneEfficiency ^ (meanOfMeans - geneData[["mean"]][[sampleNb]])
  }
  
  return(geneData)
}








CalculateRatio <- function(
  geneExperimentalData,
  geneReferenceData,
  geneReferenceName,
  geneCol = qpcr$opts$runs$geneCol,
  CqCol = qpcr$opts$runs$cqCol,
  deltaCtColName = "deltaCt",
  ratioColName = "ratio_vs_"
)
{
  if (geneExperimentalData[[geneCol]][[1]] != geneReferenceName)
  {
    ddCtName <- paste0(ratioColName, geneReferenceName)
    
    geneExperimentalData[[ddCtName]] <- NA
    
    geneExperimentalData[[ddCtName]] <- geneExperimentalData[[deltaCtColName]]/geneReferenceData[[deltaCtColName]]
    
    return(geneExperimentalData)
    
  } else return(NULL)

}







TidyTestNormalityAndHomogenityOfVarianceForOneVariable <- function(
  tidyData,
  variableName,
  valuesColumnRegex = "ratio_vs_.*",
  groupColName = qpcr$opts$plateDesign$groupColname,
  outputDir = qpcr$opts$outputDir,
  histogramsDir = "histograms"
)
{
  library(ggplot2)
  
  dirName <- paste0(outputDir, '/', histogramsDir)
  
  dir.create(dirName)
  
  valuesColnameBool <- stringr::str_detect(colnames(tidyData), pattern = valuesColumnRegex)
  
  valuesColname <- colnames(tidyData)[valuesColnameBool]
  
  groups <- unique(tidyData[[groupColName]])
  
  shapiro <- purrr::map_dfr(
    .x = groups, 
    .f = function(group){
      
      histogramName <- paste0(dirName, '/', variableName, "_", group, "_hist.png")
      dotplotName <- paste0(dirName, '/', variableName, "_", group, "_dot.png")
      
      groupSubset <- tidyData[tidyData[[groupColName]] == group,]
      
      histogram_ <- ggplot2::ggplot(
        data = groupSubset, 
        mapping = ggplot2::aes(x = !!as.symbol(valuesColname))) +
        ggplot2::geom_density()
      ggplot2::ggsave(filename = histogramName, plot = histogram_, dpi = 250)
      
      dotplot_ <- ggplot2::ggplot(
        data = groupSubset, 
        mapping = ggplot2::aes(x = !!as.symbol(valuesColname))) +
        ggplot2::geom_dotplot()
      ggplot2::ggsave(filename = dotplotName, plot = dotplot_, dpi = 250)
      
      
      
      shapiroForGroup <- broom::tidy(shapiro.test(groupSubset[[valuesColname]]))
      
      shapiroForGroup[["group"]] <- group
      
      return(shapiroForGroup)
    })
  
  formula_ <- as.formula(paste0(valuesColname, "~", groupColName))
  
  levene <- car::leveneTest(formula_, data = tidyData, center = "median")# Levene’s test: A robust alternative to the Bartlett’s test that is less sensitive to departures from normality.
  
  fligner <- broom::tidy(fligner.test(formula_, data = tidyData))# Fligner-Killeen’s test: a non-parametric test which is very robust against departures from normality.
  
  return(list("shapiro" = shapiro, "levene" = levene, "fligner" = fligner))
}







TidyWilcoxonOrTTestForOneVariable <- function(
  tidyData,
  tidyTestNormalityAndHomogenityOfVarianceOutput,
  comparisons = qpcr$opts$comparisons, ### I think we should just use list with charvectors with control and exp factor
  significanceLevel = 0.05,
  groupColname = qpcr$opts$plateDesign$groupColname,
  geneCol = qpcr$opts$runs$geneCol,
  valueColnameRegex = "ratio_vs_.*" ### !!! redo the name to valueColnameRegex
)
{
  valueColnameBool <- stringr::str_detect(colnames(tidyData), pattern = valueColnameRegex)
  
  valueColname <- colnames(tidyData)[valueColnameBool]
  
  formula_ <- as.formula(paste0(valueColname, " ~ ", groupColname))
  
  
  
  stats <- purrr::map(
    .x = comparisons, 
    .f = function(comparison){
      
      compSubset <- tidyData[tidyData[[groupColname]] %in% comparison,]
      
      areGroupsNormallyDistributed <- tidyTestNormalityAndHomogenityOfVarianceOutput$shapiro[tidyTestNormalityAndHomogenityOfVarianceOutput$shapiro[[groupColname]] %in% comparison,]$p.value
      
      
      if (tidyTestNormalityAndHomogenityOfVarianceOutput$levene[["Pr(>F)"]][[1]] > significanceLevel && areGroupsNormallyDistributed[[1]]  > significanceLevel && areGroupsNormallyDistributed[[2]]  > significanceLevel)
      {
        stat_ <- broom::tidy( t.test(formula = formula_, data = compSubset, na.action = na.omit, var.equal = T) )
        
      } else if (tidyTestNormalityAndHomogenityOfVarianceOutput$levene[["Pr(>F)"]][[1]] <= significanceLevel && areGroupsNormallyDistributed[[1]]  > significanceLevel && areGroupsNormallyDistributed[[2]]  > significanceLevel)
      {
        stat_ <- broom::tidy( t.test(formula = formula_, data = compSubset, na.action = na.omit, var.equal = F) )
        
      } else if (areGroupsNormallyDistributed[[1]]  <= significanceLevel || areGroupsNormallyDistributed[[2]]  <= significanceLevel)
      {
        stat_ <- broom::tidy( wilcox.test(formula_, data = compSubset) )
      }
      
      stat_$comparison <- paste0(comparison[[1]], "_vs_", comparison[[2]])
      stat_$gene <- compSubset[[geneCol]][[1]]
      
      return(stat_)
    })
  
  names(stats) = purrr::map(
    .x = comparisons, 
    .f = function(comp){
      paste0(comp[[1]], "_vs_", comp[[2]])
    } )
  
  return(stats)
}












TidyDrawGroupBoxplotsForEachVariable <- function(
  tidyData,
  generateIndividualPlotsBasedOnThisColumn = "gene",
  valueGroupsForSinglePlot = "group",
  valuesColumnRegex = "ratio_vs_.*", ### !!! change into regex
  outputDir = qpcr$opts$outputDir,
  boxplotDir = "boxplots_"
)
{
  dirName <- paste0(outputDir, '/', boxplotDir, generateIndividualPlotsBasedOnThisColumn)
  
  dir.create(dirName)
  
  valuesColnameBool <- stringr::str_detect(colnames(tidyData), pattern = valuesColumnRegex)
  
  valuesColumn <- colnames(tidyData)[valuesColnameBool]
  
  purrr::walk(
    .x = unique(tidyData[[generateIndividualPlotsBasedOnThisColumn]]), 
    .f = function(singlePlot){
      
      singleVariable <- tidyData[tidyData[[generateIndividualPlotsBasedOnThisColumn]] == singlePlot,]
      
      png_name <- paste0(dirName, '/', generateIndividualPlotsBasedOnThisColumn, '_', singlePlot, '.png')
      x <- ggplot(
        data = singleVariable, 
        mapping = aes(x = !!as.symbol(valueGroupsForSinglePlot), y = !!as.symbol(valuesColumn), col = !!as.symbol(valueGroupsForSinglePlot))) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(filename = png_name, plot = x)
    })
}








TidyAnovas <- function(
  tidyData,
  groupColname = qpcr$opts$plateDesign$groupColname,
  controlGroupForPostHoc = NA,
  valueColnameRegex = "ratio_vs_.*", 
  significanceLevel = 0.05,
  multipleCompMethod = NA ### !!! implement
)
{
  valueColnameBool <- stringr::str_detect(colnames(tidyData), pattern = valueColnameRegex)
  
  valueColname <- colnames(tidyData)[valueColnameBool]
  
  formula_ <- as.formula(paste0(valueColname, " ~ ", groupColname))
  
  tidyData[[groupColname]] <- as.factor(tidyData[[groupColname]])
  

  
  return_ <- list("anova" = broom::tidy(aov(formula_, data = tidyData)),
                  "welchAnova" = broom::tidy(oneway.test(formula_, data = tidyData, var.equal = F))
  )
  ### !!! Brown-Forstyhe ANOVA?
  
  if (is.na(controlGroupForPostHoc) && return_$anova$p.value[[1]] < significanceLevel)
  {
    return_[["anova_TukeyHSD"]] <- broom::tidy(TukeyHSD(aov(formula_, data = tidyData), groupColname))
    
  } else if (!is.na(controlGroupForPostHoc) && return_$anova$p.value[[1]] < significanceLevel)
  {
    return_[["anova_dunnet"]] <- as.data.frame(DescTools::DunnettTest(formula_, data = tidyData, control = controlGroupForPostHoc)[[1]])
    return_[["anova_dunnet"]]$comparison <- rownames(return_[["anova_dunnet"]])
  }
  
  if (is.na(controlGroupForPostHoc) && return_$welchAnova$p.value[[1]] < significanceLevel)
  {
    return_[["welshAnova_GH"]] <- rstatix::games_howell_test(tidyData, formula_)
    
  } else if (!is.na(controlGroupForPostHoc) && return_$welchAnova$p.value[[1]] < significanceLevel)
  {
    ### !!! to implement 
  }
  
  
  
  return_[["kruskalWallis"]] <- broom::tidy(kruskal.test(formula_, data = tidyData))
  
  if (is.na(controlGroupForPostHoc) && return_$kruskalWallis$p.value[[1]] < significanceLevel)
  {
    
    return_[["kruskalWallis_dunn"]] <- as.data.frame(DescTools::DunnTest(formula_, data = tidyData)[[1]])
    return_[["kruskalWallis_dunn"]]$comparison <- rownames(return_[["kruskalWallis_dunn"]])
    
  } else if (!is.na(controlGroupForPostHoc) && return_$kruskalWallis$p.value[[1]] < significanceLevel)
  {
    ### !!! to implement Mann-Whitney test (for KW) with multiple comparisons
  }
  
  
  return(return_)
}








PrintTidyStatisticsToExcel <- function(
  tidyStatisticSet,
  variableName,
  openxlsxWBookToWriteInto
)
{
  openxlsx::addWorksheet(openxlsxWBookToWriteInto, sheetName = variableName)
  
  currentRow <- 1
  
  for (statNb in seq_along(tidyStatisticSet)) {
    
    nbOfRowsPlusColnames <- length(tidyStatisticSet[[statNb]][[1]]) + 1
    
    openxlsx::writeData(
      openxlsxWBookToWriteInto, 
      sheet = variableName, 
      x = names(tidyStatisticSet)[[statNb]], 
      startRow = currentRow)
    
    openxlsx::writeDataTable(
      openxlsxWBookToWriteInto,
      sheet = variableName,
      x = tidyStatisticSet[[statNb]],
      colNames = TRUE,
      startRow = currentRow + 1,
      tableName = paste0(variableName, "__", names(tidyStatisticSet)[[statNb]]) )
    
    currentRow <- currentRow + nbOfRowsPlusColnames + 3
  }
}










TidyDrawCorrelationBetweenTwoVariables <- function(
  tidyData,
  samplesToCorrelateColumn = "gene",
  correlateGroupsFromThisColumn = "sample",
  valuesColumnRegex = "ratio_vs_.*", ### !!! change into regex
  outputDir = qpcr$opts$outputDir,
  correlationsDir = "correlations_"
)
{
  dirName <- paste0(outputDir, '/', correlationsDir)
  
  dir.create(dirName)
  
  valuesColnameBool <- stringr::str_detect(colnames(tidyData), pattern = valuesColumnRegex)
  
  valuesColumn <- colnames(tidyData)[valuesColnameBool]
  
  png_name <- paste0(dirName, '/', 'correlation_', samplesToCorrelateColumn, "_", correlateGroupsFromThisColumn, '.png')
  
  matrixData <- tidyData[c(correlateGroupsFromThisColumn, samplesToCorrelateColumn, valuesColumn)]
  
  matrixData <- tidyr::pivot_wider(
    data = matrixData, 
    values_from = !!as.symbol(valuesColumn),
    names_from = !!as.symbol(samplesToCorrelateColumn) )
  
  rownames(matrixData) <- matrixData[[correlateGroupsFromThisColumn]]
  matrixData[[correlateGroupsFromThisColumn]] <- NULL
  
  cor <- Hmisc::rcorr(as.matrix(matrixData))
  
  cor_plot <- ggcorrplot::ggcorrplot(cor$r,
                                     hc.order = TRUE,
                                     type = "lower",
                                     lab = TRUE,
                                     p.mat = ggcorrplot::cor_pmat(cor$r))
  
  ggsave(filename = png_name, plot = cor_plot, dpi = 250)
  
  return(cor)
}






RepairCorruptedValues <- function(
  geneQaPassed_,
  whereInGQPAreGeneData = "data",
  whereInGQPAreLadderData = c("ladder", "data"),
  whereInGQPAreNACoordinates = "valuesToNaCoordinates",
  whereInNACoordAreGeneCoords = "gene",
  whereInNACoordAreLadderCoords = "ladder",
  repLevels = as.character(qpcr$opts$runs$namesForReplicates),
  rowCoord = 1,
  colCoord = 2
)
{
  
  if (!is.null(
    geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreGeneCoords]])
  ) {
    
    for (value in seq_along(
      geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreGeneCoords]])) {
      
      geneQaPassed_[[ whereInGQPAreGeneData ]][[ geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreGeneCoords]][[value]][[rowCoord]], geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreGeneCoords]][[value]][[colCoord]] ]] <- NA
      
      geneQaPassed_[[ whereInGQPAreGeneData ]] <- AddMeanSdFlagToDF(geneQaPassed_[[ whereInGQPAreGeneData ]], repLevels)
    }
  }
  
  
  
  if (!is.null(geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreLadderCoords]])
  ) {
    
    for (value in seq_along(
      geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreLadderCoords]])) {
      
      geneQaPassed_[[ whereInGQPAreLadderData[[1]] ]][[ whereInGQPAreLadderData[[2]] ]][[geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreLadderCoords]][[value]][[rowCoord]], geneQaPassed_[[whereInGQPAreNACoordinates]][[whereInNACoordAreLadderCoords]][[value]][[colCoord]] ]] <- NA
      
      geneQaPassed_[[ whereInGQPAreLadderData[[1]] ]][[ whereInGQPAreLadderData[[2]] ]] <- AddMeanSdFlagToDF(geneQaPassed_[[ whereInGQPAreLadderData[[1]] ]][[ whereInGQPAreLadderData[[2]] ]], repLevels)
    }
  }
  
  return(geneQaPassed_)
}














# TestNormalityAndHomogenityOfVariance <- function(
#   geneData,
#   ratioColnameRegex = "ratio_vs_.*",
#   groupColName = qpcr$opts$plateDesign$groupColname
#   
# )
# {
#   groups <- unique(geneData[[groupColName]])
#   
#   ratioColnameBool <- stringr::str_detect(colnames(geneData), pattern = ratioColnameRegex)
#   
#   ratioColname <- colnames(geneData)[ratioColnameBool]
#   
#   
#   
#   shapiro <- purrr::map_dfr(
#     .x = groups, 
#     .f = function(group){
#       
#       groupSubset <- geneData[geneData[[groupColName]] == group,]
#       
#       shapiroForGroup <- broom::tidy(shapiro.test(groupSubset[[ratioColname]]))
#       
#       shapiroForGroup[["group"]] <- group
#       
#       return(shapiroForGroup)
#     })
#   
#   formula_ <- as.formula(paste0(ratioColname, "~", groupColName))
#   
#   levene <- car::leveneTest(formula_, data = geneData)# Levene’s test: A robust alternative to the Bartlett’s test that is less sensitive to departures from normality.
#   
#   fligner <- broom::tidy(fligner.test(formula_, data = geneData))# Fligner-Killeen’s test: a non-parametric test which is very robust against departures from normality.
#   
#   return(list("shapiro" = shapiro, "levene" = levene, "fligner" = fligner))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# WilcoxonOrTTest <- function(
#   ratioData,
#   normalityHomogenityDataList,
#   comparisons = qpcr$opts$comparisons, ### I think we should just use list with charvectors with control and exp factor
#   significanceLevel = 0.05,
#   groupColname = qpcr$opts$plateDesign$groupColname,
#   geneCol = qpcr$opts$runs$geneCol,
#   ratioColnameRegex = "ratio_vs_.*"
# )
# {
#   ratioColnameBool <- stringr::str_detect(colnames(ratioData), pattern = ratioColnameRegex)
#   
#   ratioColname <- colnames(ratioData)[ratioColnameBool]
#   
#   formula_ <- as.formula(paste0(ratioColname, " ~ ", groupColname))
#   
#   
#   
#   stats <- purrr::map(
#     .x = comparisons, 
#     .f = function(comparison){
#       
#       compSubset <- ratioData[ratioData[[groupColname]] %in% comparison,]
#       
#       areGroupsNormallyDistributed <- normalityHomogenityDataList$shapiro[normalityHomogenityDataList$shapiro[[groupColname]] %in% comparison,]$p.value
#       
#       
#       if (normalityHomogenityDataList$levene[["Pr(>F)"]][[1]] > significanceLevel && areGroupsNormallyDistributed[[1]]  > significanceLevel && areGroupsNormallyDistributed[[2]]  > significanceLevel)
#       {
#         stat_ <- broom::tidy( t.test(formula = formula_, data = compSubset, na.action = na.omit, var.equal = T) )
#         
#       } else if (normalityHomogenityDataList$levene[["Pr(>F)"]][[1]] <= significanceLevel && areGroupsNormallyDistributed[[1]]  > significanceLevel && areGroupsNormallyDistributed[[2]]  > significanceLevel)
#       {
#         stat_ <- broom::tidy( t.test(formula = formula_, data = compSubset, na.action = na.omit, var.equal = F) )
#         
#       } else if (areGroupsNormallyDistributed[[1]]  <= significanceLevel || areGroupsNormallyDistributed[[2]]  <= significanceLevel)
#       {
#         stat_ <- broom::tidy( wilcox.test(formula_, data = compSubset) )
#       }
#       
#       stat_$comparison <- paste0(comparison[[1]], "_vs_", comparison[[2]])
#       stat_$gene <- compSubset[[geneCol]][[1]]
#       
#       return(stat_)
#     })
#   names(stats) = purrr::map(
#     .x = comparisons, 
#     .f = function(comp){
#       paste0(comp[[1]], "_vs_", comp[[2]])
#     } )
#   
#   return(stats)
# }
