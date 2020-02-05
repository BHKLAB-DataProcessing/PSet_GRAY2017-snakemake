library(PharmacoGx)
library(readxl)
library(openxlsx)
library(tximport)
library(Biobase)

getGRAYP <-
  function (
    verbose=FALSE,
    nthread=1){
    
    options(stringsAsFactors=FALSE)
    
    
 myDirPrefix <- "/pfs/"
args = commandArgs(trailingOnly=TRUE)
rnaseq_select <- args
print(rnaseq_select)
rnaseq_results <- list()
ORCESTRA_ID = tail(rnaseq_select, n=1)

	  
tools <- grep(pattern = 'Kallisto|Salmon', x = rnaseq_select)
tools <- rnaseq_select[tools]
tools <- gsub("-", "_", tools)
transcriptome <- grep(pattern = 'Gencode|Ensembl', x = rnaseq_select)
transcriptome <- rnaseq_select[transcriptome]
tool_path = expand.grid(a = tools,b = transcriptome)
tool_path = paste0(tool_path$a, "_",tool_path$b)
	  
print(tool_path)
	  

    
    
    
    #match to curations
    
    matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating ids, we have multiple matches")
        }
        if(length(myx) == 0){return(NA_character_)}
        return(tbl[myx, returnColumn])
      })
    }
    
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]"
    
    #get curations
    
    cell_all <- read.csv(file = "/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
    curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
    rownames(curationCell) <- curationCell[ , "unique.cellid"]
    
    curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
    rownames(curationTissue) <- curationCell[ , "unique.cellid"]
    
    drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
    curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
    curationDrug <- curationDrug[,c("unique.drugid","GRAY.drugid")]
    rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
    
    
    load("/pfs/GRAY2017RawSensitivity/drug_norm_post.RData")
    
    # cell information from 2013 GRAY (cell slot)
    
    cellineinfo <- read.xlsx("/pfs/getGRAY2017/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
    cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
    rn <- cellineinfo[-1, 1]
    cn <- t(cellineinfo[1, -1])
    cn <- gsub(badchars, ".", cn)
    cellineinfo <- cellineinfo[-1, -1]
    dimnames(cellineinfo) <- list(rn, cn)
    cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])
    cellineinfo <- cellineinfo[which(!is.na(cellineinfo$Transcriptional.subtype)), ]
    
    #duplicate MB157, need to merge with MDAMB157
    cellineinfo[31,"RNASeq.availability"] <- "1"
    cellineinfo[31,"Transcriptional.subtype"] <- "Claudin-low/Basal"
    cellineinfo[31,"ERBB2.status"] <- "Claudin-low/Basal"
    cellineinfo <- cellineinfo[which(!cellineinfo$cellid == "MB157"),]
    
    cellineinfo$cellid <- as.character(matchToIDTable(ids=cellineinfo$cellid, tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    rownames(cellineinfo) <-  cellineinfo$cellid
    
    #curationCell <- curationCell[rownames(cellineinfo),]
    
    
    #published sensitivity (removed published GR data for now)
    
    
    load("/pfs/gray2017ProfilesAssemble/profiles.RData")
    
    #compile sensitivity profiles
    
    res <- res[rownames(raw.sensitivity),]
    
    sensitivity.profiles <-  data.frame("aac_recomputed" = as.numeric(res[,"AAC"]), "ic50_recomputed"=as.numeric(res[,"IC50"]), "HS"=as.numeric(res[,"HS"]), "E_inf"=as.numeric(res[,"E_inf"]), "EC50"=as.numeric(res[,"EC50"]))
    
    sensitivity.profiles$aac_recomputed <- sensitivity.profiles$aac_recomputed/100
    rownames(sensitivity.profiles) <- rownames(res)
    sensitivity.info <- sensitivity.info[rownames(raw.sensitivity),]
    
    #compute slope and add to sensitivity profiles
                                   
    slope <- NULL
    for(exp in rownames(sensitivity.info)){
      slope <- c(slope, computeSlope(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"])) #computeSlope (returns normalized slope of drug response curve)
    }
    
    names(slope) <- rownames(raw.sensitivity)
    sensitivity.profiles <- cbind(sensitivity.profiles, "slope_recomputed"=slope)
   
    #remove cross referencing cells and drugs
    
    cells.cross.reference <- read.csv("/pfs/getGRAY2017/DS0_crossreferencingCELLS.txt", header=T, stringsAsFactors=FALSE, sep="\t")
    drugs.cross.reference <- read.csv("/pfs/getGRAY2017/DS0_crossreferencingPERTURBAGENS.txt", header=T, stringsAsFactors=FALSE, sep="\t")
    
    remove.items <- which(sensitivity.info[,"cellid"] %in% cells.cross.reference[which(cells.cross.reference$COMMENT == "REMOVE"), "HEISER.NAME"])
    sensitivity.info <- sensitivity.info[-remove.items,]
    sensitivity.profiles <- sensitivity.profiles[-remove.items,]
    raw.sensitivity <- raw.sensitivity[-remove.items,,]
    old_sensitivity_info <- sensitivity.info
    sensitivity.info[, "cellid"] <- cells.cross.reference$LINCS.NAME[match( sensitivity.info[, "cellid"], cells.cross.reference[ , "HEISER.NAME"])]
    
    
    remove.items.dd <- which(sensitivity.info[,"drugid"] %in% drugs.cross.reference[which(drugs.cross.reference$COMMENT == "REMOVE"), "HEISER.NAME"])
    sensitivity.info <- sensitivity.info[-remove.items.dd,]
    sensitivity.profiles <- sensitivity.profiles[-remove.items.dd,]
    raw.sensitivity <- raw.sensitivity[-remove.items.dd,,]
    old_sensitivity_info <-  old_sensitivity_info[-remove.items.dd,]
    sensitivity.info <- as.data.frame(sensitivity.info)
    
    sensitivity.info[, "cellid"] <- gsub("-", "", sensitivity.info[, "cellid"])
    sensitivity.info[, "cellid"] <- gsub(" ", "", sensitivity.info[, "cellid"])
    sensitivity.info[, "cellid"] <- toupper(sensitivity.info[, "cellid"])
    
    sensitivity.info[, "drugid"] <- gsub("\\s*\\([^\\)]+\\)","",sensitivity.info[, "drugid"])
    
    sensitivity.info[, "cellid"] <- as.character(matchToIDTable(ids=sensitivity.info[, "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    sensitivity.info[, "drugid"] <- as.character(matchToIDTable(ids=sensitivity.info[, "drugid"], tbl=curationDrug, column = "GRAY.drugid", returnColumn = "unique.drugid"))
    
                                  
    # drug info (drug slot)

    curationDrug <- curationDrug[as.character(unique(sensitivity.info[,"drugid"])),]
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    
    #summarize rnaseq quantifications into expression sets (Kallisto)
                                   
  summarizeRnaSeq <- function (dir, 
                             features_annotation,
                             samples_annotation,
			      method) {
  library(Biobase)
  library(readr)
  library(tximport)
  
  load(features_annotation)
    
  tx2gene <- as.data.frame(cbind("transcript"=tx2gene$transcripts, "gene"=tx2gene$genes))
  
  files <- list.files(dir, recursive = TRUE, full.names = T)
  if(method=="kallisto"){
  resFiles <- grep("abundance.h5", files
  }else{
  resFiles <- grep("quant.sf", files
  }
  resFiles <- files[resFiles]
  length(resFiles)
  names(resFiles) <- basename(dirname(resFiles))
  
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData" || features_annotation == "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"){
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  } else{
  txi <- tximport(resFiles, type=method, tx2gene=tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)	  
  }
	  
  head(txi$counts[,1:5])
  dim(txi$counts)
  
  xx <- txi$abundance
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- features_gene[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- features_gene[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type="kallisto", txOut=T)
  
  xx <- txii$abundance
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  featureNames(transcript.exp) <- gsub("\\..*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp)[which(featureNames(transcript.exp) %in% rownames(features_transcript))],] 	  
  }else{
  featureNames(transcript.exp) <- gsub("\\|.*","",featureNames(transcript.exp))
  fData(transcript.exp) <- features_transcript[featureNames(transcript.exp),]
  }
  pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
  annotation(transcript.exp) <- "isoforms"
  
  xx <- txii$counts
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  if(features_annotation == "/pfs/downAnnotations/Ensembl.v99.annotation.RData"){
  featureNames(transcript.count) <- gsub("\\..*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count)[which(featureNames(transcript.count) %in% rownames(features_transcript))],]	  
  }else{
  featureNames(transcript.count) <- gsub("\\|.*","",featureNames(transcript.count))
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]	  
  }	  
  fData(transcript.count) <- features_transcript[featureNames(transcript.count),]
  pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
  annotation(transcript.count) <- "isoforms"
  
  return(list("rnaseq"=gene.exp, 
              "rnaseq.counts"=gene.count, 
              "isoforms"=transcript.exp, 
              "isoforms.counts"=transcript.count))
}

    

    rnaseq.sampleinfo <- read.csv(file="/pfs/downAnnotations/JRGraySRRMapping.csv", stringsAsFactors=FALSE, row.names=1)
    
    rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
   
    for (r in 1:length(tool_path)){
  print(tool_path[r])
  if (grep(pattern = 'Kallisto', x = tool_path[r]) > 0){
    tdir = "download_gray_rnaseqkallisto/Kallisto/"
    tool <- gsub("\\_.*","", tool_path[r])
    rnatool="kallisto"	  
  } else {
    tdir = "download_gray_rnaseqsalmon/Salmon/"
    tool <- gsub("\\_.*","", tool_path[r])
    rnatool="salmon"	  
  }
  
  
  if (length(grep(pattern = 'lift37', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33lift37.annotation.RData"
  } else if (length(grep(pattern = 'v33', x = tool_path[r])) > 0){
    annot = "/pfs/downAnnotations/Gencode.v33.annotation.RData"
  } else {
    annot = "/pfs/downAnnotations/Ensembl.v99.annotation.RData"
  }
    print(annot)
  
  rnaseq <- summarizeRnaSeq(dir=file.path(paste0(myDirPrefix, tdir, tool_path[r])),
                            features_annotation=annot,
                            samples_annotation=rnaseq.sampleinfo,
			    method = rnatool)
  rnaseq_results <- c(rnaseq_results,c(
    rnaseq <- setNames(rnaseq,  paste0(tool,".", names(rnaseq)))
  )
  )
}
    

  cellnall <- unionList(rownames(cellineinfo),rnaseq$rnaseq$cellid, sensitivity.info$cellid)
  newcells <- setdiff(cellnall, rownames(cellineinfo))
  newRows <- matrix(NA_character_, nrow=length(newcells), ncol=ncol(cellineinfo))
  # newRows <- cell.info[newcells,]

  rownames(newRows) <- newcells
  colnames(newRows) <- colnames(cellineinfo)
  newRows[,"cellid"] <- newcells

  cellineinfo <- rbind(cellineinfo, newRows)
    
    
    
  cellsPresent <- sort(unionList(sensitivity.info$cellid,rnaseq$rnaseq$cellid))
  cellineinfo <- cellineinfo[cellsPresent,]

    
    
  drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
  drug_all <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
  drug_all <- drug_all[ , c("unique.drugid", "GRAY.drugid","smiles","inchikey","cid","FDA")]
  rownames(drug_all) <- drug_all[ , "unique.drugid"]

  drug_all <- drug_all[rownames(druginfo),]
  druginfo[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]
    
z <- list()

z <- c(z,c(
  rnaseq_results
  )
)
    GRAY2017 <- PharmacoSet(molecularProfiles=z,
                            name="GRAY", 
                            cell=cellineinfo, 
                            drug=druginfo, 
                            sensitivityInfo=sensitivity.info, 
                            sensitivityRaw=raw.sensitivity, 
                            sensitivityProfiles=sensitivity.profiles, 
                            sensitivityN=NULL,
                            curationCell=curationCell, 
                            curationDrug=curationDrug, 
                            curationTissue=curationTissue, 
                            datasetType="sensitivity")
    
    saveRDS(GRAY2017,file="/pfs/out/GRAY_2017.rds")
    
    return (GRAY2017)
    
  }

getGRAYP(verbose=FALSE, nthread=1)
