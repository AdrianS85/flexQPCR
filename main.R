install.packages('pcr')
# Needs data in table format
devtools::install_github("ewallace/tidyqpcr", build_vignettes = TRUE)
BiocManager::install("ReadqPCR")



BiocManager::install("EasyqpcR") # Corrupted dependencies

source("https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R")



library(magrittr)

### plateDesignFile needs 3 columns: well, group and sample names
### complexGeneDesignFilePresent - chodzi o to, że na jednej płytce masz różne geny. wtedy musisz załączyć do każdego pliku z raw data plik z projektem genów na płytce. Jeśli F, to po prostu nazwa genu dla tej płytki. Jeśli T to plik.
### skipRows - used when there is some shit rows before actual colnames. just set it to skip all shit rows, but not colnames
### sampleRepsIn - tells if we chould look for replication information in plateDesign "plate" or in files descriptions "files"
### namesForReplicates - lista nazw dla kolejnych replikatów
qpcr <- list(
  "opts" = list(
    "sampleRepsIn" = "files",
    "plateDesign" = list(
      "file" = "plateDesign.tsv",
      "wellColname" = "well", ### This needs to have the same values as position column in runs
      "sampleColname" = "sample",
      "groupColname" = "group",
      "RepsColName" = NA),
    "runs" = list(
      "colTypes" = "ccccnccc", ### "ccccncccccnncccccnncccn" for LC96, 
      "wellCol" = "Pos",
      "cqCol" = "Cp",
      "skipRows" = 1, ### 0 for LC96
      "complexGeneDesignFilePresent" = F,
      "geneCol" = "gene",
      "emptyName" = NA, ### this cannot be changed, as !is.na is used further, do use this value to change empty names to NA
      "blankName" = "BLANK", ### !!! There may be some issues if more than one blank is present
      "ladder" = list("x10" = 10000, "x50" = 2000, "x250" = 400, "x1250" = 80, "x6250" = 16),
      "refGenes" = c("gapdh", "hmbs", "tbp"),
      "namesForReplicates" = list(1, 2, 3) ### !!! co jeśli nazwy replicatów zostały pdoane w plate file? Need to write some checking function
    )
  ))



qpcr$files = list(
  "aqp1_1" = list(
    "data" = "STQPCR3EXPAQP1.txt", 
    "genesDesign" = "aqp1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "aqp1_2" = list(
    "data" = "STQPCR3EXPAQP1(2).txt", 
    "genesDesign" = "aqp1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "aqp1_3" = list(
    "data" = "STQPCR3EXPAQP1(3).txt", 
    "genesDesign" = "aqp1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "hbb-b1_1" = list(
    "data" = "STQPCR3EXPHBB-B1.txt", 
    "genesDesign" = "hbb-b1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "hbb-b1_2" = list(
    "data" = "STQPCR3EXPHBB-B1[2].txt", 
    "genesDesign" = "hbb-b1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "hbb-b1_3" = list(
    "data" = "STQPCR3EXPHBB-B1[3].txt", 
    "genesDesign" = "hbb-b1", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "gapdh_1" = list(
    "data" = "STQPCR3EXPGAPDH.txt", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "gapdh_2" = list(
    "data" = "STQPCR3EXPGAPDH[2].txt", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "gapdh_3" = list(
    "data" = "STQPCR3EXPGAPDH[3].txt", 
    "genesDesign" = "gapdh", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "hmbs_1" = list(
    "data" = "STQPCR3EXPHMBS.txt", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "hmbs_2" = list(
    "data" = "STQPCR3EXPHMBS[2].txt", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "hmbs_3" = list(
    "data" = "STQPCR3EXPHMBS[3].txt", 
    "genesDesign" = "hmbs", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]]),
  "tbp_1" = list(
    "data" = "STQPCR3EXPTBP.txt", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[1]]),
  "tbp_2" = list(
    "data" = "STQPCR3EXPTBP[2].txt", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[2]]),
  "tbp_3" = list(
    "data" = "STQPCR3EXPTBP[3].txt", 
    "genesDesign" = "tbp", 
    "Rep" = qpcr$opts$runs$namesForReplicates[[3]])
)


### Here additional step could be added where we check if columns with given names exist
qpcr$plateDesign <- readr::read_tsv(qpcr$opts$plateDesign$file, col_names = T)



qpcr$input <- purrr::map(
  .x = qpcr$files, 
  .f = function(files){
    prepareInput(
      rawFileName = files$data, 
      genesDesign = files$genesDesign, 
      repsNameFiles = files$Rep)
  })



qpcr$inputMerged <- rlist::list.rbind(qpcr$input)



library(Hmisc)
qpcr$perGene <- purrr::map(
  .x = levels(as.factor(qpcr$inputMerged$gene)), 
  .f = function(gene_){
    prepareDataFromInput(gene_)
  })
names(qpcr$perGene) <- levels(as.factor(qpcr$inputMerged$gene))



# set lowestDilutionToUse and highestDilutionToUse with pmap
qpcr$ladderData <- purrr::map(
  .x = qpcr$perGene, 
  .f = function(gene_){
    calculateLadderStats(gene_$ladder$data$mean)
  })


CheckIfAllReferenceGenesHaveData()




test <- GetReferenceValues()

GetReferenceValues <- function(
  perGene = qpcr$perGene,
  perGeneNames = qpcr$perGene,
  refGenes = qpcr$opts$runs$refGenes,
)
{
  purrr::pmap(
    .l = list(perGene, perGeneNames), 
    .f = function(gene, geneName)
    {
      
      if (tolower(geneName) %in% tolower(refGenes))
      {
        geneMatrix <- as.matrix(gene[c("1", "2", "3")])
      }
      
    })
  

  
}



NormqPCR::selectHKs(qPCRBatch, group, method = "geNorm", minNrHKs = 2, log = TRUE, Symbols, trace = TRUE, na.rm = TRUE)


# Select reference genes

# BiocManager::install("NormqPCR"), BiocManager::install("SLqPCR") : normPCR, selectHKgenes

# calculate reference ct value from either all groups or control group. Eff POTEGA(Ctmean-Ctpróby). cdd method only if I am going to develop this

# Calculate pfaffl values compared to selected reference genes

# Visualizations

# Check normality and variance

# go to stats










prepareInput <- function(
  rawFileName,
  genesDesign,
  repsNameFiles = NA,
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
  
  colNames_ <- readr::read_tsv(rawFileName, col_types = colTypes, col_names = T, skip = skipRows, n_max = 2)
  
  rawFile <- readr::read_tsv(rawFileName, col_types = colTypes, col_names = F, skip = skipRows + 1)
  
  colnames(rawFile) <- colnames(colNames_)
  
  
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
  inputMerged = qpcr$inputMerged,
  sampleColumn = qpcr$opts$plateDesign$sampleColname,
  geneColumn = qpcr$opts$runs$geneCol,
  blankName = qpcr$opts$runs$blankName,
  ladderNames = names(qpcr$opts$runs$ladder),
  cqCol = qpcr$opts$runs$cqCol,
  RepCol = "Rep")
{
  preAnalyzedData <- inputMerged[inputMerged[[geneColumn]] == geneName & !is.na(qpcr$inputMerged[[sampleColumn]]),]
  
  repLevels <- levels(as.factor(preAnalyzedData[[RepCol]]))
  
  preAnalyzedData <- tidyr::pivot_wider(
    data = preAnalyzedData, names_from = RepCol, values_from = cqCol)
  
  cqValues <- subset(x = preAnalyzedData, select = repLevels)
  
  preAnalyzedData$mean <- rowMeans(cqValues, na.rm=TRUE)
  
  preAnalyzedData$sd <- apply(X = cqValues, MARGIN = 1, FUN = sd, na.rm=TRUE)
  
  preAnalyzedData$flag <- dplyr::case_when(
    preAnalyzedData$sd <= 0.5 ~ "OK",
    preAnalyzedData$sd < 1 & preAnalyzedData$sd > 0.5 ~ "DANGER",
    preAnalyzedData$sd >= 1 ~ "ERROR"
  )
  
  
  preAnalyzedDataList <- list(
    "data" = preAnalyzedData[preAnalyzedData[[sampleColumn]] %nin%  c(blankName, ladderNames),],
    "ladder" = list(
      "data" = preAnalyzedData[preAnalyzedData[[sampleColumn]] %in% ladderNames,]),
    "blanks" = list(
      "data" = preAnalyzedData[preAnalyzedData[[sampleColumn]] == blankName,]))
  
  return(preAnalyzedDataList)
}




calculateLadderStats <- function(
  ladderMeans,
  lowestDilutionToUse = 1,
  highestDilutionToUse = length(qpcr$opts$runs$ladder),
  ladderConcentrations = qpcr$opts$runs$ladder)
{
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
