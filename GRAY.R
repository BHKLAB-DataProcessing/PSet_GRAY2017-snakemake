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
    cellineinfo <- rbind(cellineinfo, c("SUM190PT", "breast", "NA", "NA", "0","0","0","1","1","1","0", "1" )) #add new SUM190PT cell line to this table
    rownames(cellineinfo)[84] <- "SUM190PT"
    
    cellineinfo$cellid <- as.character(matchToIDTable(ids=cellineinfo$cellid, tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
    rownames(cellineinfo) <-  cellineinfo$cellid
    
    curationCell <- curationCell[rownames(cellineinfo),]
    
    
    #published sensitivity (removed published GR data for now)
    
    
    load("/pfs/gray2017ProfilesAssemble/profiles.RData")
    
    #compile sensitivity profiles
               
    sensitivity.profiles <-  data.frame("aac_recomputed" = as.numeric(res[,"AAC"]), "ic50_recomputed"=as.numeric(res[,"IC50"]), "HS"=as.numeric(res[,"HS"]), "E_inf"=as.numeric(res[,"E_inf"]), "EC50"=as.numeric(res[,"EC50"]))
    
    sensitivity.profiles$aac_recomputed <- sensitivity.profiles$aac_recomputed
    rownames(sensitivity.profiles) <- rownames(res)
    
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
                                 tool=c("kallisto", "stringtie", "cufflinks", "rsem", "salmon"), 
                                 features_annotation,
                                 samples_annotation) {
      library(Biobase)
      library(readr)
      library(tximport)
      
      load(features_annotation)
      tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))
      
      files <- list.files(dir, recursive = TRUE, full.names = T)
      resFiles <- grep("abundance.h5", files)
      resFiles <- files[resFiles]
      length(resFiles)
      names(resFiles) <- basename(dirname(resFiles))
      
      txi <- tximport(resFiles, type="kallisto", tx2gene=tx2gene)
      head(txi$counts[,1:5])
      dim(txi$counts)
      
      xx <- txi$abundance
      gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
      fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
      pData(gene.exp) <- samples_annotation[sampleNames(gene.exp),]
      annotation(gene.exp) <- "rnaseq"
      
      xx <- txi$counts
      gene.count <- Biobase::ExpressionSet(log2(xx + 1))
      fData(gene.count) <- toil.genes[featureNames(gene.count),]
      pData(gene.count) <- samples_annotation[sampleNames(gene.count),]
      annotation(gene.count) <- "rnaseq"
      
      txii <- tximport(resFiles, type="kallisto", txOut=T)
      
      xx <- txii$abundance
      transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
      fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),]
      pData(transcript.exp) <- samples_annotation[sampleNames(transcript.exp),]
      annotation(transcript.exp) <- "isoforms"
      
      xx <- txii$counts
      transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
      fData(transcript.count) <- toil.transcripts[featureNames(transcript.count),]
      pData(transcript.count) <- samples_annotation[sampleNames(transcript.count),]
      annotation(transcript.count) <- "isoforms"
      
      return(list("rnaseq"=gene.exp, 
                  "rnaseq.counts"=gene.count, 
                  "isoforms"=transcript.exp, 
                  "isoforms.counts"=transcript.count))
    }
    

    rnaseq.sampleinfo <- read.csv("/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed/JRGraySRRMapping.csv", stringsAsFactors=FALSE, row.names=1)
    
    rnaseq.sampleinfo[ , "cellid"] <- as.character(matchToIDTable(ids=rnaseq.sampleinfo[ , "cellid"], tbl=curationCell, column = "GRAY.cellid", returnColumn = "unique.cellid"))
   
    rnaseq <- summarizeRnaSeq(dir="/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed", 
                                tool="kallisto", 
                                features_annotation="/pfs/downloadrna/Kallisto_0.43.1_processed/Kallisto_0.43.1_processed/Gencode.v23.annotation.RData",
                                samples_annotation=rnaseq.sampleinfo)
    

    
    GRAY2017 <- PharmacoSet(molecularProfiles=rnaseq,
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
    
    save(GRAY2017,file="/pfs/out/GRAY_2017.RData")
    
    return (GRAY2017)
    
  }

getGRAYP(verbose=FALSE, nthread=1)
