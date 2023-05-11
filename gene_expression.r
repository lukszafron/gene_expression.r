#! /usr/bin/env Rscript

cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

args <- commandArgs(trailingOnly=TRUE)

if(length(args)<19 | length(args)>20) {
  stop("Nineteen arguments must be supplied in the following order (Destination dir, suffix, mapping dir, gtf file, strand specificity, alpha value, FC threshold, filtering functions, threads, CSV file, sample column, grouping column, group, TXT file, list of comparisons to perform, logical argumant indicating if the samples used in the analysis are paired, name of a column with pair identifiers, logical indicator determining if the FDR correction should be aplied in the TOST analysis of equivalence, and independent factor variable(s) (one or two)).") }

args.backup <- args
read.args <- function() {
suffix <<- args[2]
inputdir <<- args[3]
gtffile <<- args[4]
stranded <<- args[5]
alpha <<- as.numeric(args[6])
FCthreshold <<- as.numeric(args[7])
Filters <<- args[8]
threads <<- args[9]
csvfile.path <<- args[10]
sample_col <<- args[11]
grouping_var <<- args[12]
group_values <<- args[13]
txt.file <<- args[14]
comparisons <<- unlist(strsplit(args[15], split = ";"))
paired.samples <<- as.logical(args[16])
pair.ident <<- args[17]
FDR <<- as.logical(args[18])
ind.factor <<- args[19]
ind.factor2 <<- args[20]
}
read.args()

if(grepl(suffix, pattern = "BAMs_wo_dups")) 
{dups <- FALSE; bam.pat <- "_sorted\\.no_dups\\.bam$"} else if(grepl(suffix, pattern = "BAMs_w_dups")) 
{dups <- TRUE; bam.pat <- "_sorted\\.bam$"} else {stop("The duplication status cannot be determined.")}
if(dups) {BAM.type <- "BAM files with duplicates"} else {BAM.type <- "BAM files without duplicates"}

if(dups) {
  workdir <- paste(args[1], "BAMS_with_dups", sep = "/")} else {
  workdir <- paste(args[1], "BAMS_without_dups", sep = "/")}

dir.create(workdir)
setwd(workdir)
save.image <- function(file){save(list=grep(ls(all.names = TRUE, envir = .GlobalEnv), pattern = "^args$", value = T, invert = T), file = file)}

save.image(file=paste0(suffix,".RData"))

library("parallel")
library("BiocParallel")
library("Rsamtools")
library("GenomicFeatures")
library("topGO")
library("openxlsx")
library("GenomicAlignments")
library("DESeq2")
library("genefilter")
library("ggplot2")
library("pals")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("ReportingTools")
library("Rsubread")
library("edgeR")
library("ReactomePA")
library("clusterProfiler")
library("pdftools")
library("dplyr")
library("ggfortify")
library("doMC")
library("foreach")
library("EnhancedVolcano")
library("stringi")
library("TOSTER")
library("tidyr")
library("data.table")

options(ggrepel.max.overlaps = Inf)
registerDoMC(threads)
register(BPPARAM = MulticoreParam(threads))

csvfile <- as.data.frame(fread(csvfile.path, header=TRUE, sep = ";", na.strings = c("", "NA"), stringsAsFactors = T))

if (grouping_var != "ALL_SAMPLES") {
  groups <- unlist(strsplit(grouping_var,split=","))
  group.values <- unlist(strsplit(group_values, split=","))
  num=0
  for (group in groups)
  {
    num=num+1
    group.column <- unlist(strsplit(group.values[grep(x=group.values, pattern=paste0("^", group, "(:|$)"))],split=":"))[1]
    group.value <- unlist(strsplit(group.values[grep(x=group.values, pattern=paste0("^", group, "(:|$)"))],split=":"))[2]
    name <- paste("df.subset",num, sep=".")
    assign(x=name,value=subset(csvfile, csvfile[[group.column]]==group.value))
  }
  if (length(groups) == 2) {
    csvfile <- merge(x=df.subset.1, y=df.subset.2)
  } else if (length(groups) == 1) {
    csvfile <- df.subset.1
  } else { 
    save.image(file=paste0(suffix,".RData"))
    stop("The number of grouping variables cannot be higher than two.")}
}

if (length(args) == 19)
{ csvfile <- subset(csvfile,!is.na(csvfile[[ind.factor]]))
} else if (length(args) == 20)
{ csvfile <- subset(csvfile,!is.na(csvfile[[ind.factor]]) & !is.na(csvfile[[ind.factor2]]))}
if (length(csvfile[[sample_col]])<2) {
  save.image(file=paste0(suffix,".RData"))
  stop("There are not enough samples to perform further analyses.")}

filenames <- grep(list.files(path = inputdir, full.names = TRUE), pattern = bam.pat, value = TRUE)
filenames.checked <- NULL
for(i in csvfile[[sample_col]]) {filenames.checked <- append(filenames.checked, grep(filenames, pattern = paste0("\\/",i, bam.pat), value = TRUE))}
filenames <- filenames.checked
# Exclude bam files the size of which does not exceed 1e6 bytes.
filenames <- filenames[file.size(filenames) > 1e6]

if(length(filenames) == 0) {cat("No valid BAM files have been found. The differential expression analysis for", BAM.type, "was not performed.\n", file = paste("Warning", suffix, "txt", sep = "."))
  } else {

bamfiles <- BamFileList(filenames)

csvfile[[sample_col]] <- as.factor(gsub(csvfile[[sample_col]], pattern = "^ *(.*) *$", replacement = "\\1"))
if(!is.na(ind.factor)) {
  csvfile[[ind.factor]] <- as.factor(gsub(csvfile[[ind.factor]], pattern = "^ *(.*) *$", replacement = "\\1"))} else if(! is.na(ind.factor2)) {
  csvfile[[ind.factor2]] <- as.factor(gsub(csvfile[[ind.factor2]], pattern = "^ *(.*) *$", replacement = "\\1"))
  } else {stop("Independent factors are missing.")}

if(grepl(x = names(bamfiles[1]),pattern = ".*_sorted(\\.no_dups)?\\.bam$")) {
  bamsamples <- as.data.frame(gsub(x=names(bamfiles),pattern="_sorted(\\.no_dups)?\\.bam$",replacement=""))} else {stop("BAM samples cannot be identified.")}
names(bamsamples) <- c("names")
sampleTable <- merge(x=bamsamples, y=csvfile, by.x="names", by.y=sample_col, sort = F)
colnames(sampleTable)[colnames(sampleTable) == "names"] <- sample_col
rownames(sampleTable) <- sampleTable[[sample_col]]

if(! (all(rownames(sampleTable) == bamsamples[["names"]]) & nrow(sampleTable) == length(bamsamples[["names"]]))) {
  save.image(file=paste0(suffix,".RData"))
  stop("Sample names in the sampleTable object and the BAM file list do not match.")}
  
if (nrow(sampleTable)<2) {
  save.image(file=paste0(suffix,".RData"))
  stop("At least two samples are necessary to perform the gene expression analysis.")
}

f.topGO <- function(geneList) {
  colMap <- function(x) {
    .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
    return(.col[match(1:length(x), order(x))])
  }
 
  wb <- createWorkbook()
  for(GOtype in c("BP", "CC", "MF")) { 
  GOdata.GOtype <- new("topGOdata",
                   ontology = GOtype,
                   allGenes = geneList,
                   geneSel = topDiffGenes,
                   annotationFun=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
                   nodeSize = 10)
  assign(paste("GOdata", GOtype, sep = "."), value = GOdata.GOtype)
  
  if(length(GOdata.GOtype@graph@nodes) >= 50) {nodesn.GOtype <- 50} else {nodesn.GOtype <- length(GOdata.GOtype@graph@nodes)}
  
  resultFisher.GOtype <- runTest(GOdata.GOtype, algorithm = "classic", statistic = "fisher")
  resultFisher.GOtype
  assign(paste("resultFisher", GOtype, sep = "."), value = resultFisher.GOtype)
  resultKS.GOtype <- runTest(GOdata.GOtype, algorithm = "classic", statistic = "ks")
  resultKS.GOtype
  assign(paste("resultKS", GOtype, sep = "."), value = resultKS.GOtype)
  resultKS.elim.GOtype <- tryCatch(runTest(GOdata.GOtype, algorithm = "elim", statistic = "ks"), error = function(e){NULL})
  resultKS.elim.GOtype
  assign(paste("resultKS.elim", GOtype, sep = "."), value = resultKS.elim.GOtype)
  if(!is.null(resultKS.elim.GOtype)) {
  allRes.GOtype.classicKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                  classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                                  orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
  allRes.GOtype.elimKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                               classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                               orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
  allRes.GOtype.Fisher <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                               classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.GOtype)} else
{
  allRes.GOtype.classicKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                  classicKS = resultKS.GOtype,
                                  orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
  allRes.GOtype.Fisher <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                               classicKS = resultKS.GOtype,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
}
  geneSig <- geneList[topDiffGenes(geneList)]
  
  AnnotatedGenes.GOtype.classicKS <- sapply(allRes.GOtype.classicKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.classicKS[[x]][AnnotatedGenes.GOtype.classicKS[[x]] %in% names(geneSig)]}
  GeneList.GOtype.classicKS <- sapply(allRes.GOtype.classicKS$GO.ID,genes_with_GO_term)
  if (length(args) == 19) {
    sink(paste("GO_top50-significant_genes.", GOtype, ".classicKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.", GOtype, ".classicKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.GOtype.classicKS)
  sink()
if(exists("allRes.GOtype.elimKS")) {  
  AnnotatedGenes.GOtype.elimKS <- sapply(allRes.GOtype.elimKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.elimKS[[x]][AnnotatedGenes.GOtype.elimKS[[x]] %in% names(geneSig)]}
  GeneList.GOtype.elimKS <- sapply(allRes.GOtype.elimKS$GO.ID,genes_with_GO_term)
  if (length(args) == 19) {
    sink(paste("GO_top50-significant_genes.", GOtype, ".elimKS",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.", GOtype, ".elimKS",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.GOtype.elimKS)
  sink()
}
  AnnotatedGenes.GOtype.Fisher <- sapply(allRes.GOtype.Fisher$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
  genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.Fisher[[x]][AnnotatedGenes.GOtype.Fisher[[x]] %in% names(geneSig)]}
  GeneList.GOtype.Fisher <- sapply(allRes.GOtype.Fisher$GO.ID,genes_with_GO_term)
  if (length(args) == 19) {
    sink(paste("GO_top50-significant_genes.", GOtype, ".Fisher",".(",ind.factor,").",suffix,".txt", sep = ""))
  } else {
    sink(paste("GO_top50-significant_genes.", GOtype, ".Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".txt", sep = ""))
  }
  print(GeneList.GOtype.Fisher)
  sink()

      sheet.name <- paste("GO-top_50", GOtype, "classicKS", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)
      writeData(allRes.GOtype.classicKS, wb = wb, sheet = sheet.name, rowNames = F)

      sheet.name <- paste("GO-top_50", GOtype, "elimKS", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)

if(exists("allRes.GOtype.elimKS")) {
      writeData(allRes.GOtype.elimKS, wb = wb, sheet = sheet.name, rowNames = F)
}    
      sheet.name <- paste("GO-top_50", GOtype, "Fisher", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)
      writeData(allRes.GOtype.Fisher, wb = wb, sheet = sheet.name, rowNames = F)
suppressWarnings(rm(resultKS.elim.GOtype, allRes.GOtype.elimKS))
}

  if (length(args) == 19) {
    saveWorkbook(wb, file = paste("GO-top_50",".(",ind.factor,").",suffix,".xlsx", sep = ""), overwrite = T)
  } else {
    saveWorkbook(wb, file = paste("GO-top_50",".(",ind.factor,"+",ind.factor2,").",suffix,".xlsx", sep = ""), overwrite = T)
  }

  if (length(args) == 19) {
    pdf(paste("GO-p-values_comparison",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-p-values_comparison",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
for(GOtype in c("BP", "CC", "MF")) {
if(!is.null(paste("resultKS.elim", GOtype, sep = "."))) {
  pValue.classic <- score(get(paste("resultKS", GOtype, sep = ".")))
  pValue.elim <- score(get(paste("resultKS.elim", GOtype, sep = ".")))[names(pValue.classic)]
  gstat <- termStat(get(paste("GOdata", GOtype, sep = ".")), names(pValue.classic))
  gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  gCol <- colMap(gstat$Significant)
  plot(x=pValue.classic, y=pValue.elim, xlab = "p-value - KS test (classic)", ylab = "p-value - KS test (elim)", pch = 19, cex = gSize, col = gCol, title(main=paste("Gene ontology analysis", GOtype, sep = " - ")))
} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "resultKS.elim ERROR", sep = " - "))}
}
  dev.off()
   
  if (length(args) == 19) {
    pdf(paste("GO-diagrams_top_10_test_KS_classic",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_KS_classic",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  
  for(GOtype in c("BP", "CC", "MF")) {
   tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultKS",GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
}
  dev.off()
   
  if (length(args) == 19) {
    pdf(paste("GO-diagrams_top_10_test_KS_elim",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_KS_elim",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  for(GOtype in c("BP", "CC", "MF")) {
  tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultKS.elim", GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
  }
  dev.off()
  
  if (length(args) == 19) {
    pdf(paste("GO-diagrams_top_10_test_Fisher",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("GO-diagrams_top_10_test_Fisher",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  for(GOtype in c("BP", "CC", "MF")) {
   tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultFisher", GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
           error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
  }
  dev.off()
}

f_Reactome <- function(res.df.alpha) {
  for (exp.type in c("upregulated genes", "downregulated genes")) {
    if (exp.type == "upregulated genes") {
      res.df.alpha.FCthreshold <- subset(res.df.alpha, FC >= FCthreshold)
    } else {
      res.df.alpha.FCthreshold <- subset(res.df.alpha, FC <= 1/FCthreshold)
    }
    geneList <- res.df.alpha.FCthreshold$FC
    names(geneList) <- res.df.alpha.FCthreshold$Entrez
    geneList <- geneList[!is.na(names(geneList))]
    geneList <- sort(geneList, decreasing = TRUE)
    geneList <- geneList[!duplicated(names(geneList))]
    geneList.df <- as.data.frame(geneList)
    colnames(geneList.df) <- "FC"
    if (nrow(geneList.df) > 0) {
      rownames(geneList.df) <- as.vector(mapIds(org.Hs.eg.db,
                                                keys=rownames(geneList.df),
                                                column="SYMBOL",
                                                keytype="ENTREZID",
                                                multiVals="first"))}
    geneList.df.name <- paste("geneList",exp.type, sep = ".")
    assign(geneList.df.name, value = geneList.df)
    
    de <- names(geneList)
    if (length(de) > 0) {
      de.name <- paste("de",exp.type, sep = ".")
      assign(de.name, value = de)
      
      temp.x <- enrichPathway(gene=get(de.name), organism = "human", pAdjustMethod = "BH", pvalueCutoff=alpha, readable=TRUE)
      x.name <- paste("x",exp.type, sep = ".")
      x.name.df <- paste("x",exp.type, "df", sep = ".")
      assign(x.name, value = temp.x)
      rm(temp.x)
      try(p1 <- dotplot(get(x.name), showCategory=50) + labs(title = paste("Pathway enrichment analysis", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      p1.name <- paste('p1', exp.type, sep = ".")
      if(exists('p1')) {if(nrow(p1$data) > 0) {assign(p1.name, value = p1)
        p1.pathways <- as.character(p1$data[order(p1$data$GeneRatio, decreasing = T),][["Description"]][1:if(nrow(p1$data)<3){nrow(p1$data)} else {3}])
        if(exp.type == "upregulated genes") {
          pdf(paste("Reactome analysis", suffix, "01_1", "pdf", sep = "."), height = 10, width = 13)} else if(exp.type == "downregulated genes") {
            pdf(paste("Reactome analysis", suffix, "07_1", "pdf", sep = "."), height = 10, width = 13)}
        for(i in p1.pathways) {try(expr = {pp1 <- viewPathway(i)
        pp1 <- pp1 + labs(title = paste(exp.type, i, sep = ": ")) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
        print(pp1)})}
        dev.off()
      }}
      try(p2 <- heatplot(get(x.name), showCategory = 50, foldChange = geneList) + labs(title = paste("Pathway enrichment heatmap", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5)))
      p2.name <- paste('p2', exp.type, sep = ".")
      if(exists('p2')) {assign(p2.name, value = p2)}
      try(p3 <- emapplot(get(x.name), showCategory = 50, color = "p.adjust") + labs(title = paste("Pathway enrichment map", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      p3.name <- paste('p3', exp.type, sep = ".")
      if(exists('p3')) {assign(p3.name, value = p3)}
      try(p4 <- cnetplot(get(x.name), categorySize="qvalue", foldChange=geneList) + labs(title = paste("Complex associations map", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      p4.name <- paste('p4', exp.type, sep = ".")
      if(exists('p4')) {assign(p4.name, value = p4)}
      
      tmp.df <- as.data.frame(get(x.name))
      assign(x.name.df, value = tmp.df)
      
      try({y <- gsePathway(geneList, organism = "human", pvalueCutoff = alpha, pAdjustMethod = "BH", by = "fgsea"); res.fgsea <- as.data.frame(y); p5 <- emapplot(y, showCategory = 50, color = "p.adjust") + labs(title = paste("Gene set enrichment analysis", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))})
      if(exists("y")) { if(nrow(y) > 0) {
        for(i in seq(1,length(y@result$core_enrichment))) {y@result$core_enrichment[i] <- paste(as.character(mapIds(org.Hs.eg.db, keys = unlist(strsplit(y@result$core_enrichment[i], split = "/")), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")),collapse = "/")}}}
      p5.name <- paste('p5', exp.type, sep = ".")
      if(exists('p5')) {assign(p5.name, value = p5)}
      
      y.name <- paste("y",exp.type, sep = ".")
      y.name.df <- paste("y",exp.type, "df", sep = ".")
      if(exists("y")) {assign(y.name, value = y)}
      if(exists("y")) {assign(y.name.df, value = as.data.frame(get(y.name)))}
      geneList.Symbols <- geneList
      if (length(geneList.Symbols) > 0) {
        names(geneList.Symbols) <- mapIds(org.Hs.eg.db, keys = names(geneList.Symbols), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")}
      if(exists("y")) {if(nrow(get(y.name.df)) > 0) { 
        p6 <- heatplot(get(y.name), showCategory = 50, foldChange = geneList.Symbols) + labs(title = paste("Gene set enrichment heatmap", exp.type, sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5))}}
      p6.name <- paste('p6', exp.type, sep = ".")
      if(exists('p6')) {assign(p6.name, value = p6)}
      rm(p1, p2, p3, p4, p5, p6, y)
    }}
  if(! exists("de.upregulated genes")) {`de.upregulated genes` <- NULL}
  if (! exists("de.downregulated genes")) {`de.downregulated genes` <- NULL}
  de.list.full <- list(`de.upregulated genes`, `de.downregulated genes`)
  names(de.list.full) <- c("upregulated genes", "downregulated genes")
  
  try(compareClustersRes <- compareCluster(de.list.full, fun="enrichPathway", organism = "human", pAdjustMethod = "BH", pvalueCutoff=alpha, readable=TRUE))
  try(p7 <- dotplot(compareClustersRes, showCategory=50) + labs(title = paste("Pathway enrichment analysis - group comparison", sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  try(compareClustersRes.df <- as.data.frame(compareClustersRes))
  
  if (exists("p1.upregulated genes")){
    pdf(paste("Reactome analysis", suffix, "01", "pdf", sep = "."), height = 10, width = 13)
    print(`p1.upregulated genes`)
    dev.off()}
  if(exists('p2.upregulated genes')) { if(nrow(`p2.upregulated genes`$data) > 260){
    pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = nrow(`p2.upregulated genes`$data)/20)} else {
      pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 13)
    }} else {
      pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 10)}
  if(exists('p2.upregulated genes')){
    print(`p2.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
  dev.off()
  
  pdf(paste("Reactome analysis", suffix, "03", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p3.upregulated genes')) {
    print(`p3.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
  dev.off()
  pdf(paste("Reactome analysis", suffix, "04", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p4.upregulated genes')) {
    print(`p4.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
  dev.off()
  pdf(paste("Reactome analysis", suffix, "05", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p5.upregulated genes')) {
    print(`p5.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", "upregulated genes", sep = "-"))}
  dev.off()
  if(exists('p6.upregulated genes')) { if (nrow(`p6.upregulated genes`$data) > 260){
    pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = nrow(`p6.upregulated genes`$data)/20)} else {
      pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 13)
    }} else {
      pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 10)}
  if(exists('p6.upregulated genes')) {
    print(`p6.upregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map.", "upregulated genes", sep = "-"))}
  dev.off()
  if (exists("p1.downregulated genes")){
    pdf(paste("Reactome analysis", suffix, "07", "pdf", sep = "."), height = 10, width = 13)
    print(`p1.downregulated genes`)
    dev.off()}
  if(exists('p2.downregulated genes')) { if(nrow(`p2.downregulated genes`$data) > 260){
    pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = nrow(`p2.downregulated genes`$data)/20)} else {
      pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 13)
    }} else {
      pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 10)}
  if(exists('p2.downregulated genes')){
    print(`p2.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
  dev.off()
  pdf(paste("Reactome analysis", suffix, "09", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p3.downregulated genes')) {
    print(`p3.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
  dev.off()
  pdf(paste("Reactome analysis", suffix, "10", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p4.downregulated genes')) {
    print(`p4.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
  dev.off()
  pdf(paste("Reactome analysis", suffix, "11", "pdf", sep = "."), height = 13, width = 13)
  if(exists('p5.downregulated genes')) {
    print(`p5.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", "downregulated genes", sep = "-"))}
  dev.off()
  if(exists('p6.downregulated genes')) { if(nrow(`p6.downregulated genes`$data) > 260){
    pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = nrow(`p6.downregulated genes`$data)/20)} else {
      pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width =  13)
    }} else {
      pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = 10)}
  if(exists('p6.downregulated genes')){
    print(`p6.downregulated genes`)} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map.", "downregulated genes", sep = "-"))}
  dev.off()
  if(exists("p7")) {
    if(nrow(p7$data)>60) {p7.height = floor(nrow(p7$data)/6)} else {p7.height = 10}
    pdf(paste("Reactome analysis", suffix, "13", "pdf", sep = "."), height = p7.height, width = 13)
    print(p7)
    dev.off()}
  
  pdffiles <- sort(list.files(pattern = paste("Reactome analysis", gsub(suffix, pattern = "\\+", replacement = "\\\\+"), ".*", "pdf", sep = ".")), method = "radix")
  
  if (length(args) == 19) {
    reactome.pdf.name = paste0("Reactome analysis", ".(",ind.factor,").",suffix,".pdf", sep = "")
  } else {
    reactome.pdf.name = paste0("Reactome analysis", ".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = "")
  }
  pdf_combine(pdffiles, output = reactome.pdf.name)
  unlink(pdffiles)
  
  if (length(args) == 19) {
    reactome.xls.name = paste0("Reactome analysis", ".(",ind.factor,").",suffix,".xlsx", sep = "")
  } else {
    reactome.xls.name = paste0("Reactome analysis", ".(",ind.factor,"+",ind.factor2,").",suffix,".xlsx", sep = "")
  }
  
  unlink(reactome.xls.name)
  wb <- createWorkbook()
  
  if(exists("geneList.upregulated genes")) {
    
    sheet.name <- "Upregulated.genes"
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`geneList.upregulated genes`, wb = wb, sheet = sheet.name, rowNames = T)}
  
  if(exists("x.upregulated genes.df")) {
    
    sheet.name <- "Path.enrich.upregulated.genes"
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`x.upregulated genes.df`, wb = wb, sheet = sheet.name, rowNames = F)}
  
  if(exists("y.upregulated genes.df")) {
    
    sheet.name <- "Gene.set.enrich.upregulated.genes"
    sheet.name <- substr(sheet.name, 1, 31)
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`y.upregulated genes.df`, wb = wb, sheet = sheet.name, rowNames = F)}
  
  if(exists("geneList.downregulated genes")) {
    
    sheet.name <- "Downregulated.genes"
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`geneList.downregulated genes`, wb = wb, sheet = sheet.name, rowNames = T)}
  
  if(exists("x.downregulated genes.df")) {
    
    sheet.name <- "Path.enrich.downregulated.genes"
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`x.downregulated genes.df`, wb = wb, sheet = sheet.name, rowNames = F)}
  
  if(exists("y.downregulated genes.df")) {
    
    sheet.name <- "Gene.set.enrich.downregulated.genes"
    sheet.name <- substr(sheet.name, 1, 31)
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(`y.downregulated genes.df`, wb = wb, sheet = sheet.name, rowNames = F)}
  
  if(exists("compareClustersRes.df")) {
    
    sheet.name <- "Path.enrich.group.comparison"
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(compareClustersRes.df, wb = wb, sheet = sheet.name, rowNames = F)}
  
  saveWorkbook(wb, file = reactome.xls.name, overwrite = T)
}

if(grepl(suffix, pattern="DESEQ2")) {
  if(! file.exists(paste0(suffix,".SE.transformed.filtered",".RData"))) {
    if (! file.exists(paste0(suffix,".SE.transformed",".RData"))) {
      
      if(!paired.samples) {
        if (length(args) == 19) { factors.formula <- as.formula(paste0("~",ind.factor)) } else { factors.formula <- as.formula(paste0("~",ind.factor,"+",ind.factor2)) }
      } else {
        sampleTable[[pair.ident]] <- as.factor(sampleTable[[pair.ident]])
        if (length(args) == 19) { factors.formula <- as.formula(paste0("~",ind.factor, "+", pair.ident)) } else { factors.formula <- as.formula(paste0("~",ind.factor,"+",ind.factor2, "+", pair.ident)) }
      }
      txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
      txdb
      tbg <- transcriptsBy(txdb, by="gene")

      if(stranded == "yes") {
        se <- summarizeOverlaps(features=tbg, reads=bamfiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=FALSE,
                                fragments=TRUE )
      } else {
        se <- summarizeOverlaps(features=tbg, reads=bamfiles,
                                mode="Union",
                                singleEnd=FALSE,
                                ignore.strand=TRUE,
                                fragments=TRUE )	
      }
      se
      # dim(se)
      # assayNames(se)
      rowRanges(se)
      str(metadata(rowRanges(se)))
      colnames(se) <- bamsamples$names
      if(is.null(colnames(se))) {stop("Column names of the se object are missing.")}
      if(! all(colnames(se) == rownames(sampleTable))) {stop("The sample names in the se and sampleTable objects do not match.")}
      colData(se) <- DataFrame(sampleTable)
      colData(se)
      head(assay(se), 3)
      colSums(assay(se))

      ds <- DESeqDataSet(se, design = factors.formula)
      if (length(ds[[sample_col]]) <= 30)
      {
        rld.e <- tryCatch(rld <- rlog(ds, blind = FALSE),error=function (e) { rlog(ds, blind = TRUE) })
        if (exists("rld")) { rm(rld.e) } else { rld <- rld.e}
      } else if (length(ds[[sample_col]]) > 30)
      {
        vsd <- vst(ds, blind = FALSE)
      }
      save.image(file=paste0(suffix,".SE.transformed",".RData"))
    } else {
      cat("Loading a pre-existing RData file:", paste0(suffix,".SE.transformed",".RData...\n"))
      load(file=paste0(suffix,".SE.transformed",".RData"))
      read.args() }
      
    if (exists("rld")) 
    {
      if (Filters != "UNFILTERED")
      {
        ffunctions <- unlist(strsplit(Filters,split=";"))
        for(i in seq(from=1, to=length(ffunctions),by=1)) { name <- paste("ff",i,sep=".")
        assign(name,eval(parse(text=ffunctions[i])))}
        ff.number <- length(grep(objects(),pattern="ff.[0-9]+"))
        if (ff.number == 1) 
        { selGenes <- genefilter(assay(rld),filterfun(ff.1))
        } else if (ff.number == 2) 
        { selGenes <- genefilter(assay(rld),filterfun(ff.1,ff.2))
        } else if (ff.number == 3)
        { selGenes <- genefilter(assay(rld),filterfun(ff.1,ff.2,ff.3))
        } else { 
          save.image(file=paste0(suffix,".SE.transformed",".RData"))
          stop("The number of filtering functions is incorrect.") }
        eset <- assay(rld[selGenes,])
        rld[rownames(assay(rld)) %in% rownames(eset), ] -> rld
        ds[rownames(assay(ds)) %in% rownames(eset), ] -> ds
      }
      save.image(file=paste0(suffix,".SE.transformed.filtered",".RData"))
    } else if(exists("vsd")) 
    {
      if (Filters != "UNFILTERED")
      {
        ffunctions <- unlist(strsplit(Filters,split=";"))
        for(i in seq(from=1, to=length(ffunctions),by=1)) { name <- paste("ff",i,sep=".")
        assign(name,eval(parse(text=ffunctions[i])))}
        ff.number <- length(grep(objects(),pattern="ff.[0-9]+"))
        if (ff.number == 1) 
        { selGenes <- genefilter(assay(vsd),filterfun(ff.1))
        } else if (ff.number == 2) 
        { selGenes <- genefilter(assay(vsd),filterfun(ff.1,ff.2))
        } else if (ff.number == 3)
        { selGenes <- genefilter(assay(vsd),filterfun(ff.1,ff.2,ff.3))
        } else { 
          save.image(file=paste0(suffix,".SE.transformed",".RData"))
          stop("The number of filtering functions is incorrect.") }
        eset <- assay(vsd[selGenes,])
        vsd[rownames(assay(vsd)) %in% rownames(eset), ] -> vsd
        ds[rownames(assay(ds)) %in% rownames(eset), ] -> ds
      }
      save.image(file=paste0(suffix,".SE.transformed.filtered",".RData"))
    }	
  } else { 
    cat("Loading a pre-existing RData file:", paste0(suffix,".SE.transformed.filtered",".RData...\n"))
    load(file=paste0(suffix,".SE.transformed.filtered",".RData"))
  read.args() }
  
  DEAPP <- "DESeq2"
  RUNID <- gsub(x=suffix, pattern="(^[^.]*)(.*)", replacement="\\1")
  for (i in grep(ls(), pattern="(^rld$|^vsd$)", value=TRUE)) {
    
    pdf(paste("Read counts histogram-",suffix,".pdf", sep = ""))
    hist(counts(ds), xlim = c(0,quantile(counts(ds), c(0.95))), breaks = 100000, col = "gray", main = paste0("DESeq2 - histogram of cumulative read counts in the ", RUNID, " run\n(centiles: 0-95), median value: ", median(counts(ds))), xlab = "Read count")
    dev.off()
    
    if (length(args) == 19) {
      pdf(paste("PCA-plot.",i,".(",ind.factor,").",suffix,".pdf", sep = ""))
      plot.pca <- plotPCA(get(i), intgroup = ind.factor, returnData=TRUE)
      percentVar <- round(100 * attr(plot.pca, "percentVar"))
      if(NROW(plot.pca$name)<=25) {
        ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = name, shape = group)) +
          geom_point(size =3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(color=sample_col, shape=ind.factor, title = "Principal component analysis (PCA) plot") +
          coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
          theme(plot.title = element_text(hjust = 0.5))
        print(ggplot.1)
        dev.off()
      } else 	{
        ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size =3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(color=ind.factor, title = "Principal component analysis (PCA) plot") +
          coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
          theme(plot.title = element_text(hjust = 0.5))
        print(ggplot.1)
        dev.off()
      }
    } else {
      ind.interactions <- interaction(sampleTable[[ind.factor]], sampleTable[[ind.factor2]])
      pdf(paste("PCA-plot.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
      plot.pca <- plotPCA(get(i), intgroup = c(ind.factor,ind.factor2), returnData=TRUE)
      percentVar <- round(100 * attr(plot.pca, "percentVar"))
      if(NROW(plot.pca$name)<=25) {
        ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = name, if(length(ind.interactions) <= 6) {shape = group})) +
          geom_point(size =3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(color=sample_col, shape=paste(ind.factor,ind.factor2,sep=":"), title = "Principal component analysis (PCA) plot") +
          coord_fixed() + scale_color_manual(values = as.vector(cols25())) +
          theme(plot.title = element_text(hjust = 0.5))
        print(ggplot.1)
        dev.off()
      } else 	{
        ggplot.1 <- ggplot(plot.pca, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size =3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          labs(color=paste(ind.factor,ind.factor2,sep=":"), title = "Principal component analysis (PCA) plot") +
          coord_fixed() + if(length(ind.interactions) <= 25) {scale_color_manual(values = as.vector(cols25()))} +
          theme(plot.title = element_text(hjust = 0.5))
        print(ggplot.1)
        dev.off()
      }
    }
    
    topVarGenes <- head(order(rowVars(assay(get(i))), decreasing = TRUE), 20)
    mat  <- assay(get(i))[ topVarGenes, ]
    mat  <- mat - rowMeans(mat)
    mat.ens <- mat 
    rownames(mat) <- gsub(rownames(mat),pattern="\\..*$", replacement="")
    rownames(mat) <- as.vector(mapIds(org.Hs.eg.db,
                                      keys=rownames(mat),
                                      column="SYMBOL",
                                      keytype="ENSEMBL",
                                      multiVals="first"))
    for(j in seq(from=1, to=nrow(mat))) {rownames(mat)[j] [is.na(rownames(mat)[j])] <- rownames(mat.ens)[j]}
    
    if(ncol(mat)>35) {width <- ceiling(ncol(mat)/5)} else {width <- 7}
  if(! paired.samples) {  
    if (length(args) == 19) {
      anno <- as.data.frame(colData(get(i))[, ind.factor],row.names = as.vector(sampleTable[[sample_col]]))
      colnames(anno) <- ind.factor
      pdf(paste("TopVarGenes.",i,".(",ind.factor,").",suffix,".pdf", sep = ""), width = width)
    } else {
      anno <- as.data.frame(colData(get(i))[, c(ind.factor,ind.factor2)],row.names = as.vector(sampleTable[[sample_col]]))
      colnames(anno) <- c(ind.factor, ind.factor2)
      pdf(paste("TopVarGenes.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = width)
    }
  } else {
    if (length(args) == 19) {
      anno <- as.data.frame(colData(get(i))[, c(ind.factor,pair.ident)],row.names = as.vector(sampleTable[[sample_col]]))
      colnames(anno) <- c(ind.factor, pair.ident)
      pdf(paste("TopVarGenes.",i,".(",ind.factor,").",suffix,".pdf", sep = ""), width = width)
    } else {
      anno <- as.data.frame(colData(get(i))[, c(ind.factor,ind.factor2,pair.ident)],row.names = as.vector(sampleTable[[sample_col]]))
      colnames(anno) <- c(ind.factor, ind.factor2, pair.ident)
      pdf(paste("TopVarGenes.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = width)
    }
  }
    pheatmap(mat, annotation_col = anno, cellwidth = 10, cellheight = 10)
    dev.off()
    
    sampleDists <- dist(t(assay(get(i))))
    sampleDists
    sampleDistMatrix <- as.matrix( sampleDists )
 
    #colnames(sampleDistMatrix) <- NULL
    if(ncol(sampleDistMatrix)>35) {width <- ceiling(ncol(sampleDistMatrix)/5)} else {width <- 7}
    height <- width
    
    if (length(args) == 19) {
      pdf(paste("Sample_distance.",i,".(",ind.factor,").",suffix,".pdf", sep = ""), width = width, height = height)
    } else {
      pdf(paste("Sample_distance.",i,".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = width, height = height)
    }
    
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             annotation_col = anno,
             cellwidth = 10, cellheight = 10
    )
    dev.off()
  }

  for(comparison in comparisons)
  {
    res.dir <- stri_reverse(sub(stri_reverse(comparison), pattern = ":", replacement = "_sv_"))
    ost.dir <- paste(paste0("Comparison:", res.dir), paste(paste("FC_threshold",FCthreshold, sep = "="), paste("Alpha_value", alpha, sep = "="), sep = ";"), sep = ",")
    dir.create(ost.dir, showWarnings = F)
    setwd(dir = ost.dir)
  
  if (length(ds[[sample_col]])<2) {
    stop('The analysis cannot be performed for less than two samples.')
  } 	else if (exists("rld.e"))
  {
    res <- data.frame(
      assay(rld),
      avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
      rLogFC = assay(rld)[,2] - assay(rld)[,1], FC = 2^(assay(rld)[,2] - assay(rld)[,1]) )
    res$Ensembl.ID <- gsub("\\..*","",rownames(res))

    res$symbol <- mapIds(org.Hs.eg.db,
                         keys=res$Ensembl.ID,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    res$Entrez <- mapIds(org.Hs.eg.db,
                         keys=res$Ensembl.ID,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  } else
  {
    dds <- DESeq(ds)
    res <- results(dds, contrast = unlist(strsplit(comparison, split = ":")))
    res <- res[!is.na(res$padj), , drop = F]
    res$FC <- 2**(res$log2FoldChange)
    res$Ensembl.ID <- gsub("\\..*","",rownames(res))
    res$symbol <- mapIds(org.Hs.eg.db,
                         keys=res$Ensembl.ID,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
    res$Entrez <- mapIds(org.Hs.eg.db,
                         keys=res$Ensembl.ID,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
    res.df <- as.data.frame(res)
    res.df.FC <- subset(res.df, FC >= FCthreshold | FC <= 1/FCthreshold)
    geneList <- as.vector(res.df.FC$padj)
    names(geneList) <- as.vector(res.df.FC$symbol)
    geneList <- geneList[!is.na(geneList)]
    topDiffGenes <- function(pvalue) {return(pvalue < alpha) }
    if(sum(topDiffGenes(geneList)) > 0) 
    {
      try(f.topGO(geneList))
    }
    res.df.alpha <- subset(res.df, padj < alpha)
    f_Reactome(res.df.alpha = res.df.alpha)
  }
  print(head(res))
  #summary(res)
  #mcols(res)
  if (! exists("rld.e")) {
    if (length(args) == 19) {
      pdf(paste("MA-plot",".(",ind.factor,").",suffix,".pdf", sep = ""))
    } else {
      pdf(paste("MA-plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
    }
    tryCatch(expr=DESeq2::plotMA(res, ylim = c(-5, 5)),
             error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="MA plot function returned an ERROR")}})
    dev.off()
  }
 
  if (! exists("rld.e")) {
    if (length(args) == 19) {
      pdf(paste("P-value-histogram",".(",ind.factor,").",suffix,".pdf", sep = ""))
    } else {
      pdf(paste("P-value-histogram",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
    }
    hist(res$padj, col = "grey50", border = "white", main = "Histogram of BH-adjusted p-values", xlab = "BH-adjusted p-value")
  } else if (exists("rld.e"))
  {
    if (length(args) == 19) {
      pdf(paste("FC-values-histogram",".(",ind.factor,").",suffix,".pdf", sep = ""))
    } else {
      pdf(paste("FC-values-histogram",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
    }
    hist(res$FC, col = "grey50", border = "white", main = "Histogram of Fold change values", xlab = "Fold change")
  }
  dev.off()
  
if (! exists("rld.e")) {
  if (length(args) == 19) {
    pdf(paste("Volcano plot",".(",ind.factor,").",suffix,".pdf", sep = ""), width = 10)
  } else {
    pdf(paste("Volcano plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = 10)
  }
  if(nrow(res) >=200){volcano.nrow.limit <- 200} else {volcano.nrow.limit <- nrow(res)}
  if(abs(res$log2FoldChange)[order(abs(res$log2FoldChange), decreasing = T)][volcano.nrow.limit] < log2(FCthreshold)) {l2FC.threshold <- log2(FCthreshold)} else 
  {l2FC.threshold <- abs(res$log2FoldChange)[order(abs(res$log2FoldChange), decreasing = T)][volcano.nrow.limit]}
  if(res$padj[order(res$padj)][volcano.nrow.limit] > alpha) {padj.threshold <- alpha} else {padj.threshold <- res$padj[order(res$padj)][volcano.nrow.limit]}
res.volcano <- res
res.volcano$symbol[is.na(res.volcano$symbol)] <- res.volcano$Ensembl.ID[is.na(res.volcano$symbol)]
volcano1 <- EnhancedVolcano(res.volcano,
                  lab = res.volcano$symbol,
                  x = "log2FoldChange",
                  y = "padj",
                  pCutoff = padj.threshold,
                  FCcutoff = l2FC.threshold,
                  pointSize = 1.5,
                  labSize = 2.5,
                  shape = c(6, 6, 19, 16),
                  title = paste(DEAPP, "differential expression analysis results", sep = " - "),
                  subtitle = NULL,
                  caption = paste0("log2FC cutoff = ", formatC(l2FC.threshold, digits = 3), "\nBH-adjusted p-value cutoff = ", formatC(padj.threshold, digits = 3)),
                  legendPosition = "top",
                  legendLabSize = 14,
                  col = c("grey30", "forestgreen", "royalblue", "red2"),
                  colAlpha = 0.9,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5) + 
    theme(plot.title = element_text(hjust = 0.5))
  tryCatch(expr = {print(volcano1)}, error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("ERROR: Generation of the Volcano plot failed."))}})
  dev.off()
}
  if (! exists("rld.e")) {
    resOrdered <- res[order(res$pvalue),]
    head(resOrdered)
    resOrderedDF <- as.data.frame(resOrdered)
    if (length(args) == 19) {
      write.table(resOrderedDF, file = paste("DESeq2_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
      htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""),
                            reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""))
    } else {
      write.table(resOrderedDF, file = paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
      htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
                            reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
    }
    publish(resOrderedDF, htmlRep)
    url <- finish(htmlRep)
    browseURL(url, browser="firefox")
  } else if (exists("rld.e"))
  {
    resOrdered <- res[order(res$FC),]
    head(resOrdered)
    if (length(args) == 19) {
      write.table(resOrdered, file = paste("DESeq2_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
      htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""),
                            reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,").",suffix, sep = ""))
    } else {
      write.table(resOrdered, file = paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
      htmlRep <- HTMLReport(shortName="DESeq2 analysis report", title=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
                            reportDirectory=paste("DESeq2_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
    }
    publish(resOrdered, htmlRep)
    url <- finish(htmlRep)
    browseURL(url, browser="firefox")
  }
setwd("../")
}
  
  if (! exists("rld.e")) {
    
  if(ncol(counts(dds))>35) {width <- ceiling(ncol(counts(dds))/5)} else {width <- 7}
    
  pdf(paste("Read counts normalization-",suffix,".pdf", sep = ""), height = 10, width = width)
  par(mfrow=c(2,1))
  boxplot(counts(dds)+1, col = "lightblue", las = 2, cex.names = 1, log="y")
  title(main=paste0("DESeq2 - raw read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
  boxplot(counts(dds, normalized=TRUE)+1, col = "lightblue", las = 2, cex.names = 1, log="y")
  title(main=paste0("DESeq2 - normalized read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
  dev.off()  

  normalized.counts <- counts(dds, normalized=TRUE)
  rownames(normalized.counts) <- gsub(rownames(normalized.counts),pattern="\\..*$", replacement="")
  rownames(normalized.counts) <- as.vector(mapIds(org.Hs.eg.db,
                                                  keys=rownames(normalized.counts),
                                                  column="SYMBOL",
                                                  keytype="ENSEMBL",
                                                  multiVals="first"))
  normalized.counts <- normalized.counts[!is.na(rownames(normalized.counts)),]
  normalized.counts <- normalized.counts[order(rownames(normalized.counts)),]
  normalized.counts <- t(normalized.counts)
  normalized.counts <- as.data.frame(normalized.counts)
  }
  
  print("ColSums unfiltered:")
  print(colSums(assay(se)))
  print("ColSums filtered:")
  print(colSums(counts(ds)))
  final.RData.name <- paste0(suffix,".SE.transformed.filtered",".RData")
  save.image(file = final.RData.name)
} else 
  if(grepl(suffix, pattern="EDGER")) {
  if (! file.exists(paste0(suffix,".SE.transformed",".RData"))) {
    if (! file.exists(paste0(suffix,".SE.RData"))) {
      if(stranded == "yes") {
        se <- featureCounts(files=filenames,annot.ext=gtffile,
                            isGTFAnnotationFile=TRUE,GTF.featureType="transcript",GTF.attrType="gene_id",
                            isPairedEnd = TRUE,
                            tmpDir=tempdir(),
                            nthreads=threads,
                            strandSpecific = 1 )
      } else {
        se <- featureCounts(files=filenames,annot.ext=gtffile,
                            isGTFAnnotationFile=TRUE,GTF.featureType="transcript",GTF.attrType="gene_id",
                            isPairedEnd = TRUE,
                            tmpDir=tempdir(),
                            nthreads=threads,
                            strandSpecific = 0 )
      }
      
      if(! all(foreach(i = seq(1,nrow(sampleTable)), .combine = c) %do% {grepl(filenames[i], pattern = rownames(sampleTable)[i])})) {
        save.image(file=paste0(suffix,".RData"))
        stop("The BAM file names and sample names do not match.")
      }
      colnames(se$counts) <- rownames(sampleTable)
      edgeR.dge <- DGEList(counts=se$counts)
      edgeR.dge <- calcNormFactors(edgeR.dge)
      if (Filters != "UNFILTERED")
      {
        keep <- filterByExpr(edgeR.dge)
        edgeR.counts.filtered <- edgeR.dge$counts[keep,]
        edgeR.dge <- DGEList(counts=edgeR.counts.filtered)
        edgeR.dge <- calcNormFactors(edgeR.dge)
      }
      save.image(file=paste0(suffix,".SE.RData")) } else {
        cat("Loading a pre-existing RData file:", paste0(suffix,".SE.RData...\n"))
        load(file=paste0(suffix,".SE.RData"))
    read.args() }
    
    for(i in colnames(edgeR.dge$counts)) {temp <- edgeR.dge$counts[,i]/(edgeR.dge$samples$lib.size[rownames(edgeR.dge$samples) %in% i]*edgeR.dge$samples$norm.factors[rownames(edgeR.dge$samples) %in% i])*1e6; assign(value=temp,x=paste("col",i, sep="."))}
    colList <- grep(objects(),pattern="^col\\..+$",value=TRUE)
    edgeR.counts.norm <- sapply(mget(colList),cbind)
    colnames(edgeR.counts.norm) <- colnames(edgeR.dge$counts)
    rownames(edgeR.counts.norm) <- rownames(edgeR.dge$counts)
    edgeR.counts.norm.df <- as.data.frame(edgeR.counts.norm, stringsAsFactors=FALSE)
    rm(list=grep(objects(),pattern="^col\\..+$",value=TRUE))
    
    if(!paired.samples) {
      if (length(args) == 19) { 
        design <- model.matrix(~ 0 + sampleTable[[ind.factor]]) } else {
        design <- model.matrix(~ 0 + sampleTable[[ind.factor]] + sampleTable[[ind.factor2]]) }
    } else {
      sampleTable[[pair.ident]] <- as.factor(sampleTable[[pair.ident]])
      if (length(args) == 19) { 
        design <- model.matrix(~ 0 + sampleTable[[ind.factor]] + sampleTable[[pair.ident]]) } else {
        design <- model.matrix(~ 0 + sampleTable[[ind.factor]] + sampleTable[[ind.factor2]] + sampleTable[[pair.ident]]) }
    }
    
    colnames(design) <- sub(sub(sub(sub(colnames(design), pattern = "\\[\\[ind.factor\\]\\]", replacement = paste0(".", ind.factor, ".")), pattern = "\\[\\[ind.factor2\\]\\]", replacement = paste0(".", ind.factor2, ".")), pattern = "\\[\\[pair.ident\\]\\]", replacement = paste0(".", pair.ident, ".")), pattern = "sampleTable\\.", replacement = "")
    
    edgeR.dge <- estimateDisp(edgeR.dge, design)
    rld <- cpm(edgeR.dge, log=TRUE) ## Shows filtered, normalized results in a log2 scale
    #cpm(edgeR.dge, log=FALSE) ## Shows filtered, normalized results in a linear scale
    save.image(file=paste0(suffix,".SE.transformed",".RData"))
  } else {
    cat("Loading a pre-existing RData file:", paste0(suffix,".SE.transformed",".RData...\n"))
    load(file=paste0(suffix,".SE.transformed",".RData"))
    read.args() }
    
  DEAPP <- "edgeR"
  RUNID <- gsub(x=suffix, pattern="(^[^.]*)(.*)", replacement="\\1")
  
  pdf(paste("Read counts histogram-",suffix,".pdf", sep = ""))
  hist(edgeR.dge$counts, xlim = c(0,quantile(edgeR.dge$counts, c(0.95))), breaks = 100000, col = "gray", main = paste0("edgeR - histogram of cumulative read counts in the ", RUNID, " run\n(centiles: 0-95), median value: ", median(edgeR.dge$counts)), xlab = "Read count")
  dev.off()
  
  if(ncol(edgeR.dge$counts)>35) {width <- ceiling(ncol(edgeR.dge$counts)/5)} else {width <- 7}

    pdf(paste("Read counts normalization-",suffix,".pdf", sep = ""), height = 10, width = width)
    par(mfrow=c(2,1))
    boxplot(edgeR.dge$counts+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
    title(main=paste0("edgeR - raw read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
    boxplot(edgeR.counts.norm.df+1, col = "lightgreen", las = 2, cex.names = 1, log="y")
    title(main=paste0("edgeR - normalized read counts in the ", RUNID, " run."), xlab=sample_col, ylab="Read counts in a log10 scale")
    dev.off()
  
  plot.mds <- try(plotMDS(edgeR.dge))
  
if(class(plot.mds) != "try-error") { 
  
  dev.off()
  unlink('Rplots.pdf')
  plot.mds.df <- as.data.frame(plot.mds$cmdscale.out)
  if (length(args) == 19) {
    pdf(paste("MDS-plot",".(",ind.factor,").",suffix,".pdf", sep = ""))
    if(nrow(plot.mds.df)<=25) {
      ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = rownames(plot.mds.df), shape = sampleTable[[ind.factor]])) + 
        geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
        xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
        labs(title = "Multidimensional scaling (MDS) plot", shape = ind.factor, color = sample_col) + 
        theme(plot.title = element_text(hjust = 0.5))
      print(ggplot.1)
    } else 	{
      ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = sampleTable[[ind.factor]])) + 
        geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
        xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
        labs(title = "Multidimensional scaling (MDS) plot", color = ind.factor) + 
        theme(plot.title = element_text(hjust = 0.5))
      print(ggplot.1)
    }
    dev.off()
  } else {
    ind.interactions <- interaction(sampleTable[[ind.factor]], sampleTable[[ind.factor2]])
    pdf(paste("MDS-plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
    if(nrow(plot.mds.df)<=25) {
      ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = rownames(plot.mds.df), if(length(ind.interactions) <= 6) {shape = ind.interactions})) + 
        geom_point(size = 3) + coord_fixed() + scale_color_manual(values = as.vector(cols25())) + 
        xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
        labs(title = "Multidimensional scaling (MDS) plot", shape = paste(ind.factor, ind.factor2, sep = "+"), color = sample_col) + 
        theme(plot.title = element_text(hjust = 0.5))
      print(ggplot.1)
    } else 	{
      ggplot.1 <- ggplot(plot.mds.df, aes(x = V1, y = V2, color = ind.interactions)) + 
        geom_point(size = 3) + coord_fixed() + if(length(ind.interactions) <= 25) {scale_color_manual(values = as.vector(cols25()))} + 
        xlab("Leading logFC dim1") + ylab("Leading logFC dim2") + 
        labs(title = "Multidimensional scaling (MDS) plot", color = paste(ind.factor, ind.factor2, sep = "+")) + 
        theme(plot.title = element_text(hjust = 0.5))
      print(ggplot.1)
    }
    dev.off()
  }
} else {
  if (length(args) == 19) {
    pdf(paste("MDS-plot",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("MDS-plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=plot.mds[1], cex.main = 1)
  dev.off()
}
  
  topVarGenes <- head(order(rowVars(rld), decreasing = TRUE), 20)
  mat  <- rld[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  mat.ens <- mat
  rownames(mat) <- gsub(rownames(mat),pattern="\\..*$", replacement="")
  rownames(mat) <- as.vector(mapIds(org.Hs.eg.db,
                                    keys=rownames(mat),
                                    column="SYMBOL",
                                    keytype="ENSEMBL",
                                    multiVals="first"))
  for(i in seq(from=1, to=nrow(mat))) {rownames(mat)[i] [is.na(rownames(mat)[i])] <- rownames(mat.ens)[i]}
  
  if(ncol(mat)>35) {width <- ceiling(ncol(mat)/5)} else {width <- 7}
if(! paired.samples) {
    if (length(args) == 19) {
    anno <- as.data.frame(sampleTable[, ind.factor],row.names = as.vector(sampleTable[[sample_col]]))
    colnames(anno) <- ind.factor
    pdf(paste("TopVarGenes.rld",".(",ind.factor,").",suffix,".pdf", sep = ""), width = width)
  } else {
    anno <- as.data.frame(sampleTable[, c(ind.factor, ind.factor2)],row.names = as.vector(sampleTable[[sample_col]]))
    colnames(anno) <- c(ind.factor, ind.factor2)
    pdf(paste("TopVarGenes.rld",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = width)
  }
} else {
  if (length(args) == 19) {
    anno <- as.data.frame(sampleTable[, c(ind.factor, pair.ident)],row.names = as.vector(sampleTable[[sample_col]]))
    colnames(anno) <- c(ind.factor, pair.ident)
    pdf(paste("TopVarGenes.rld",".(",ind.factor,").",suffix,".pdf", sep = ""), width = width)
  } else {
    anno <- as.data.frame(sampleTable[, c(ind.factor, ind.factor2, pair.ident)],row.names = as.vector(sampleTable[[sample_col]]))
    colnames(anno) <- c(ind.factor, ind.factor2, pair.ident)
    pdf(paste("TopVarGenes.rld",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = width)
  }
}
  pheatmap(mat, annotation_col = anno, cellwidth = 10, cellheight = 10)
  dev.off()
  
  sampleDists <- dist(t(rld))
  sampleDists
  sampleDistMatrix <- as.matrix(sampleDists)
 
  #colnames(sampleDistMatrix) <- NULL
  if(ncol(sampleDistMatrix)>35) {width <- ceiling(ncol(sampleDistMatrix)/5)} else {width <- 7}
  height <- width
  
    if (length(args) == 19) {
    pdf(paste("Sample_distance.rld",".(",ind.factor,").",suffix,".pdf", sep = ""), height = height, width = width)
  } else {
    pdf(paste("Sample_distance.rld",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), height = height, width = width)
  }
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           annotation_col = anno,
           cellwidth = 10, cellheight = 10
  )
  dev.off()
  
  if (length(sampleTable[[sample_col]])<3) {
    save.image(file=paste0(suffix,".SE.transformed",".RData"))
    stop('The analysis cannot be performed for less than three samples.')}
  
  for(comparison in comparisons)
  {
    res.dir <- stri_reverse(sub(stri_reverse(comparison), pattern = ":", replacement = "_sv_"))
    ost.dir <- paste(paste0("Comparison:", res.dir), paste(paste("FC_threshold",FCthreshold, sep = "="), paste("Alpha_value", alpha, sep = "="), sep = ";"), sep = ",")
    dir.create(ost.dir, showWarnings = F)
    setwd(dir = ost.dir)
  
  comparison.elements <- unlist(strsplit(comparison, split = ":"))
  contrast <- makeContrasts(paste(paste(comparison.elements[1], comparison.elements[2], sep = "."), paste(comparison.elements[1], comparison.elements[3], sep = "."), sep = "-"), levels = design)
  fit <- glmQLFit(edgeR.dge, design)
  qlf <- glmQLFTest(fit, contrast = contrast)
  res <- topTags(qlf, n=nrow(qlf$table))
  res$table <- res$table[!is.na(res$table$FDR), , drop = F]
  res$table$FC <- 2**res$table$logFC
  res$table$Ensembl.ID <- gsub("\\..*","",rownames(res))

  res$table$symbol <- mapIds(org.Hs.eg.db,
                             keys=res$table$Ensembl.ID,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
  res$table$Entrez <- mapIds(org.Hs.eg.db,
                             keys=res$table$Ensembl.ID,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
  res.df <- as.data.frame(res$table)
  res.df.FC <- subset(res.df, FC >= FCthreshold | FC <= 1/FCthreshold)
  geneList <- as.vector(res.df.FC$FDR)
  names(geneList) <- as.vector(res.df.FC$symbol)
  geneList <- geneList[!is.na(geneList)]
  topDiffGenes <- function(pvalue) {return(pvalue < alpha) }
  if(sum(topDiffGenes(geneList)) > 0) 
  {
    try(f.topGO(geneList))
  }
  res.df.alpha <- subset(res.df, FDR < alpha)
  f_Reactome(res.df.alpha = res.df.alpha)

  print(head(res))
  if (length(args) == 19) {
    pdf(paste("P-value-histogram",".(",ind.factor,").",suffix,".pdf", sep = ""))
  } else {
    pdf(paste("P-value-histogram",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""))
  }
  hist(res$table$FDR, col = "grey50", border = "white", main = "Histogram of BH-adjusted p-values", xlab = "BH-adjusted p-value")
  dev.off()
  
  if (length(args) == 19) {
    pdf(paste("Volcano plot",".(",ind.factor,").",suffix,".pdf", sep = ""), width = 10)
  } else {
    pdf(paste("Volcano plot",".(",ind.factor,"+",ind.factor2,").",suffix,".pdf", sep = ""), width = 10)
  }
  if(nrow(res) >=200){volcano.nrow.limit <- 200} else {volcano.nrow.limit <- nrow(res)}
  if(abs(res$table$logFC)[order(abs(res$table$logFC), decreasing = T)][volcano.nrow.limit] < log2(FCthreshold)) {l2FC.threshold <- log2(FCthreshold)} else 
  {l2FC.threshold <- abs(res$table$logFC)[order(abs(res$table$logFC), decreasing = T)][volcano.nrow.limit]}
  if(res$table$FDR[order(res$table$FDR)][volcano.nrow.limit] > alpha) {FDR.threshold <- alpha} else {FDR.threshold <- res$table$FDR[order(res$table$FDR)][volcano.nrow.limit]}
  res.volcano <- res
  res.volcano$table$symbol[is.na(res.volcano$table$symbol)] <- res.volcano$table$Ensembl.ID[is.na(res.volcano$table$symbol)]
  volcano1 <- EnhancedVolcano(res.volcano$table,
                              lab = res.volcano$table$symbol,
                              x = "logFC",
                              y = "FDR",
                              pCutoff = FDR.threshold,
                              FCcutoff = l2FC.threshold,
                              pointSize = 1.5,
                              labSize = 2.5,
                              shape = c(6, 6, 19, 16),
                              title = paste(DEAPP, "differential expression analysis results", sep = " - "),
                              subtitle = NULL,
                              caption = paste0("log2FC cutoff = ", formatC(l2FC.threshold, digits = 3), "\nBH-adjusted p-value cutoff = ", formatC(FDR.threshold, digits = 3)),
                              legendPosition = "top",
                              legendLabSize = 14,
                              col = c("grey30", "forestgreen", "royalblue", "red2"),
                              colAlpha = 0.9,
                              drawConnectors = TRUE,
                              widthConnectors = 0.5) + 
    theme(plot.title = element_text(hjust = 0.5))
  tryCatch(expr = {print(volcano1)}, error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("ERROR: Generation of the Volcano plot failed."))}})
  dev.off()

  resOrderedDF <- res$table[order(res$table$FDR),]
  head(resOrderedDF)
  if (length(args) == 19) {
    write.table(resOrderedDF, file = paste("edgeR_analysis_results",".(",ind.factor,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
    htmlRep <- HTMLReport(shortName="edgeR analysis report", title=paste("edgeR_analysis_results",".(",ind.factor,").",suffix, sep = ""),
                          reportDirectory=paste("edgeR_analysis_results",".(",ind.factor,").",suffix, sep = ""))
  } else {
    write.table(resOrderedDF, file = paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = ""), sep = ";", dec = ".", row.names = FALSE)
    htmlRep <- HTMLReport(shortName="edgeR analysis report", title=paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""),
                          reportDirectory=paste("edgeR_analysis_results",".(",ind.factor,"+",ind.factor2,").",suffix, sep = ""))
  }
  publish(resOrderedDF, htmlRep)
  url <- finish(htmlRep)
  browseURL(url, browser="firefox")
  
  setwd("../")
  }
  
  normalized.counts <- cpm(edgeR.dge, log=FALSE)
  rownames(normalized.counts) <- gsub(rownames(normalized.counts),pattern="\\..*$", replacement="")
  rownames(normalized.counts) <- as.vector(mapIds(org.Hs.eg.db,
                                                  keys=rownames(normalized.counts),
                                                  column="SYMBOL",
                                                  keytype="ENSEMBL",
                                                  multiVals="first"))
  normalized.counts <- normalized.counts[!is.na(rownames(normalized.counts)),]
  normalized.counts <- normalized.counts[order(rownames(normalized.counts)),]
  normalized.counts <- t(normalized.counts)
  normalized.counts <- as.data.frame(normalized.counts)
  
  print("ColSums unfiltered:")
  print(colSums(se$counts))
  print("ColSums filtered:")
  print(colSums(edgeR.dge$counts))
  final.RData.name <- paste0(suffix,".SE.transformed",".RData")
  
  save.image(file = final.RData.name)
  }

group.comparison <- function(normalized.counts, txt.file = NULL){

  if(! is.null(txt.file)) {
  gene.signature <- gsub(txt.file, pattern = "^.*\\/", replacement = "")
  con <- file(txt.file)
  gene.list1 <- unique(gsub(readLines(con), pattern = "\\s", replacement = ""))
  close(con)
  
  gene.list1.conv <- NULL
  for(i in gene.list1) {new.symbol <- alias2Symbol(i, species = "Hs"); if(length(new.symbol) == 0) {gene.list1.conv <- append(gene.list1.conv, values = i)} else {gene.list1.conv <- append(gene.list1.conv, values = new.symbol)}}
  gene.list1.conv <- unique(gene.list1.conv)
  
  if(grepl(suffix, pattern = "DESEQ2")){
    if(exists('vsd')) {rld.mx <- assay(vsd)} else 
      if(exists('rld')) {rld.mx <- assay(rld)}} else 
  if(grepl(suffix, pattern = "EDGER")) {rld.mx <- rld}
  
  rownames(rld.mx) <- sub(rownames(rld.mx), pattern = "\\.[0-9]+$", replacement = "")
  
  rownames(rld.mx) <- as.vector(mapIds(org.Hs.eg.db,
                                       keys=rownames(rld.mx),
                                       column="SYMBOL",
                                       keytype="ENSEMBL",
                                       multiVals="first"))
  
  rld.mx <- rld.mx[!is.na(rownames(rld.mx)),]
  rld.mx.sel <- rld.mx[rownames(rld.mx) %in% gene.list1.conv, , drop = FALSE]

if(!is.null(nrow(rld.mx.sel))) {
  if(nrow(rld.mx.sel) > 1) {
  rld.mx.sel <- rld.mx.sel[!duplicated(rownames(rld.mx.sel)), , drop = FALSE]

  if(ncol(rld.mx.sel)>35) {width <- ceiling(ncol(rld.mx.sel)/5)} else {width <- 7}
  if(nrow(rld.mx.sel)>35) {height <- ceiling(nrow(rld.mx.sel)/5)} else {height <- 7}

  pdf(paste("1", gene.signature, "gene signature comparison heatmap.pdf", sep = "-"), height = height, width = width)
  heatmap0 <- pheatmap(rld.mx.sel, annotation_col = anno, main = paste("Heatmap for the", gene.signature, "gene signature"), cellwidth = 10, cellheight = 10)
  dev.off()
  
  rld.mx.sel.txp <- t(rld.mx.sel)
  rld.mx.sel.txp.merged <- merge(x = rld.mx.sel.txp, y = anno, by.x = 0, by.y = 0)
  rownames(rld.mx.sel.txp.merged) <- rld.mx.sel.txp.merged$Row.names
if(! paired.samples) {
  if(length(args) == 19){
  delcols <- c("Row.names", ind.factor)} else {
  delcols <- c("Row.names", ind.factor, ind.factor2)}
  } else {
  if(length(args) == 19){
  delcols <- c("Row.names", ind.factor, pair.ident)} else {
  delcols <- c("Row.names", ind.factor, ind.factor2, pair.ident)}
  }
  
  rld.mx.sel.txp.merged.pca.data <- prcomp(rld.mx.sel.txp.merged %>% dplyr::select(!delcols))

  rld.mx.sel.txp.merged.df <- as.data.frame(rld.mx.sel.txp.merged)
  rld.mx.sel.txp.merged.df <- rld.mx.sel.txp.merged.df[match(rownames(rld.mx.sel.txp.merged.pca.data$x), rownames(rld.mx.sel.txp.merged.df)), , drop = F]
  
  if(length(args) == 19){
    rld.mx.sel.txp.merged.df[[sample_col]] <- rownames(rld.mx.sel.txp.merged.df)
    if(nrow(rld.mx.sel.txp.merged.df) > 25) {
      ind.length <- length(levels(rld.mx.sel.txp.merged.df[[ind.factor]]))} else {
      ind.length <- length(rld.mx.sel.txp.merged.df[[sample_col]])
      }
    pdf(paste("0", gene.signature, "gene signature comparison PCA plot.pdf", sep = "-"), width = 10, height = 10)
      if(ncol(rld.mx.sel.txp.merged.pca.data$x) > 1) {
      ap1 <- autoplot(rld.mx.sel.txp.merged.pca.data, data = rld.mx.sel.txp.merged.df, size = 3, colour = if(nrow(rld.mx.sel.txp.merged.df) > 25) {ind.factor} else {sample_col}, shape = if(nrow(rld.mx.sel.txp.merged.df) > 25) {pty=20} else {ind.factor}) + 
      scale_color_manual(values = if (ind.length > 25) {rainbow(ind.length)} else {as.vector(cols25())}) + coord_fixed() + 
      labs(title = paste("PCA plot for the", gene.signature, "gene signature")) + 
      theme(plot.title = element_text(hjust = 0.5))
      try(print(ap1))} else {
        plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="ERROR: Too few dimensions to draw a PCA plot.")
      }
    dev.off()
} else {
  ind.length <- c(length(levels(rld.mx.sel.txp.merged.df[[ind.factor]])), length(levels(rld.mx.sel.txp.merged.df[[ind.factor2]])))
  pdf(paste("0", gene.signature, "gene signature comparison PCA plot.pdf", sep = "-"),  width = 10, height = 10)
  if(ncol(rld.mx.sel.txp.merged.pca.data$x) > 1) {
  if(ind.length[1] <= 6) {
    colors.n <- ind.length[2]
    ap1 <- autoplot(rld.mx.sel.txp.merged.pca.data, data = rld.mx.sel.txp.merged.df, size = 3, shape = ind.factor, colour = ind.factor2)} else
  if(ind.length[2] <= 6) {
    colors.n <- ind.length[1]
    ap1 <- autoplot(rld.mx.sel.txp.merged.pca.data, data = rld.mx.sel.txp.merged.df, size = 3, shape = ind.factor2, colour = ind.factor)} else {
    colors.n <- ind.length[1] * ind.length[2]
    ap1 <- autoplot(rld.mx.sel.txp.merged.pca.data, data = rld.mx.sel.txp.merged.df, size = 3, colour = interaction(ind.factor, ind.factor2))}
    ap1 <- ap1 + scale_color_manual(values = if (colors.n > 25) {rainbow(colors.n)} else {as.vector(cols25())}) + coord_fixed() + 
    labs(title = paste("PCA plot for the", gene.signature, "gene signature")) + 
    theme(plot.title = element_text(hjust = 0.5))
    try(print(ap1))} else {
      plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main="ERROR: Too few dimensions to draw a PCA plot.")
    }
  dev.off()}
  gene.signature.pdf.name <- paste(gene.signature, "gene_signature_comparison_summary", suffix, "pdf", sep = ".")
  unlink(gene.signature.pdf.name)
  pdffiles2 <- list.files(pattern = paste0(paste(gene.signature, "gene signature comparison", sep = "-"), ".*\\.pdf$"))
  if(length(pdffiles2) >0) {
    pdf_combine(pdffiles2, output = gene.signature.pdf.name)
    unlink(pdffiles2)}
}}} else {
  gene.signature <- "ALL_GENES" }
  
  if(! is.null(txt.file)) {
  normalized.counts <- normalized.counts[, colnames(normalized.counts) %in% gene.list1.conv, drop = F]
  }
if(ncol(normalized.counts) >0) {
  
  if (length(args) == 19) {
    file.name <- paste("Normalized_counts.gene_signature:", gene.signature, ".(",ind.factor,").",suffix,".csv", sep = "")
  } else {
    file.name <- paste("Normalized_counts.gene_signature:", gene.signature, ".(",ind.factor,"+",ind.factor2,").",suffix,".csv", sep = "")
  }
  write.table(x = normalized.counts, row.names = T, col.names = NA, quote = F, sep = ";", dec = ".", file = file.name)
  
  if(! is.na(ind.factor2)) {
  anno <- anno %>% unite(col = "Combined", ind.factor:ind.factor2, sep = ".")
  anno[["Combined"]] <- as.factor(anno[["Combined"]])
  combined.var.name <- paste(ind.factor, ind.factor2, sep = ".")
  colnames(anno)[colnames(anno) == "Combined"] <- combined.var.name
  indfactor <- combined.var.name} else {
    indfactor <- ind.factor
  }

  if(length(levels(anno[[indfactor]])) >1) {
  
  gene.signature.xlsx.name <- paste(gene.signature, "gene_signature_comparison_summary", paste("alpha", alpha, sep = ":"), paste("FDR", FDR, sep = ":"), suffix, "xlsx", sep = ".")
  
  wb <- createWorkbook()
  combinations <- combn(levels(anno[[indfactor]]), m = 2)
  for(iter in seq(1,ncol(combinations))) {
    comparison.group <- substr(paste(combinations[,iter], collapse = "_vs_"),1,31)
    addWorksheet(wb = wb, sheetName = comparison.group)
    anno.sel <- anno[anno[[indfactor]] %in% combinations[, iter], , drop = F]
    anno.sel[[indfactor]] <- factor(anno.sel[[indfactor]])
    normalized.counts.sel <- normalized.counts[match(rownames(anno.sel), rownames(normalized.counts)), , drop = F]
    if(! all(rownames(normalized.counts.sel) == rownames(anno.sel))) {stop("Sample names in the normalized counts and anno objects do not match.")}

  missing.group <- combinations[,iter][! combinations[,iter] %in% anno.sel[[indfactor]]]
  if(length(missing.group) == 0) {
  min.N <- min(table(anno.sel[[indfactor]]))
  min.N.limit <- 2
  if(min.N >= min.N.limit) {
    all.tested.genes <- unique(colnames(normalized.counts.sel))
     
    f.data.stats <- function(data) {rbind(colMeans(data[,-1], na.rm = T), 
                                          Sds <- colSds(as.matrix(data[,-1]), na.rm = T), 
                                          pooledSd <- sqrt(sum(sapply(Sds, function(x) {x**2}))/2),
                                          colSums(!is.na(data[,-1])),
                                          if(paired.samples) {
                                            sddif <- sd(data[,2] - data[,3])
                                            cor.res <- cor.test(x = data[,2], data[,3])
                                            if(is.na(cor.res$estimate)) {cor.res$estimate <- 0} else
                                            if(cor.res$estimate == 1) {cor.res$estimate <- 0.9999999}
                                            rbind(rep(sddif,2),
                                            rep(cor.res$estimate,2)) } )
                                          }
    ff.10 <- pOverA(0.101,0, na.rm = T)
    
    value.conv <- function(value) {
      if(!is.na(value)) {
        value <- as.numeric(value)
        if(abs(value) < 0.001) {
          value <- formatC(value, format = "e", digits = 3)} else 
          {
            value <- round(value,4)}
      return(as.numeric(value))} else {return(value)}}
    
    TOST.res <- foreach(gene.name = all.tested.genes, .combine = rbind) %do% {
      if(! paired.samples) {
      df.tmp <- data.frame(sample = rownames(anno.sel), gene = normalized.counts.sel[[gene.name]], indfactor = anno.sel[[indfactor]])
      df.tmp <- reshape(data = df.tmp, timevar = "indfactor", direction = "wide", idvar = "sample")
      } else {
      df.tmp <- data.frame(sample = rownames(anno.sel), gene = normalized.counts.sel[[gene.name]], indfactor = anno.sel[[indfactor]], Pairs = anno.sel[[pair.ident]])
      df.tmp <- reshape(df.tmp[,-1], idvar = "Pairs", timevar = "indfactor", direction = "wide")
      }
      
      if(any(sapply(df.tmp[,-1],ff.10))) {
        df.tmp.stats.two <- f.data.stats(data = df.tmp)
        df.tmp.stats.two[2,][df.tmp.stats.two[2,] == 0] <- 1e-50
        
        sink("/dev/null")
        sink(type = "m")
        
        if(paired.samples) {
          eqbound <- suppressMessages(round(max(powerTOSTpaired.raw(alpha=alpha, statistical_power=0.8, N = min.N, sdif = df.tmp.stats.two[5,1])),2)) } else {
          eqbound <- suppressMessages(round(max(powerTOSTtwo.raw(alpha=alpha, statistical_power=0.8, N = min.N, sdpooled = df.tmp.stats.two[3,1])),2)) }
        
        sink()
        
        if(eqbound == 0) {eqbound <- 0.1}
        
        if(paired.samples) {
        t1 <- tsum_TOST(m1=df.tmp.stats.two[1,1],
                        m2=df.tmp.stats.two[1,2],
                        sd1=df.tmp.stats.two[2,1],
                        sd2=df.tmp.stats.two[2,2],
                        n1=df.tmp.stats.two[4,1],
                        n2=df.tmp.stats.two[4,2],
                        low_eqbound=-eqbound,
                        high_eqbound=eqbound,
                        alpha = alpha,
                        var.equal=FALSE,
                        eqbound_type = "raw",
                        paired = TRUE,
                        r12 = df.tmp.stats.two[6,1])
        } else {
          t1 <- tsum_TOST(m1=df.tmp.stats.two[1,1],
                          m2=df.tmp.stats.two[1,2],
                          sd1=df.tmp.stats.two[2,1],
                          sd2=df.tmp.stats.two[2,2],
                          n1=df.tmp.stats.two[4,1],
                          n2=df.tmp.stats.two[4,2],
                          low_eqbound=-eqbound,
                          high_eqbound=eqbound,
                          alpha = alpha,
                          var.equal=FALSE,
                          eqbound_type = "raw",
                          paired = FALSE)
                        }

        group.names <- gsub(colnames(df.tmp.stats.two), pattern = "gene\\.", replacement = "")        
        v1 <- c(t1$method, group.names, df.tmp.stats.two[4,], df.tmp.stats.two[1,], df.tmp.stats.two[2,], -eqbound, eqbound, if(paired.samples) {df.tmp.stats.two[6,1]} else {NA}, t1$TOST["TOST Lower","p.value"], t1$TOST["TOST Upper","p.value"], t1$TOST["t-test","p.value"])
        v1[4:length(v1)] <- sapply(v1[4:length(v1)], value.conv)
        
        if(max(c(t1$TOST["TOST Lower","p.value"], t1$TOST["TOST Upper","p.value"])) < alpha & t1$TOST["t-test","p.value"] >= alpha) {
          v1 <- c(gene.name, v1, "Yes")
        } else {
          v1 <- c(gene.name, v1, "No")
        }
        v1
    }
    }
    if(!is.null(TOST.res)) {
      if(is.vector(TOST.res)) {TOST.res <- as.data.frame(t(TOST.res), stringsAsFactors = F)} else 
      {TOST.res <- as.data.frame(TOST.res, stringsAsFactors = F)}
      
      rownames(TOST.res) <- NULL
      colnames(TOST.res) <- c("Gene name", "Test type", "Group1 name", "Group2 name", "Group1 N", "Group2 N", "Group1 mean", "Group2 mean", "Group1 SD", "Group2 SD", "Lower equivalence bounds", "Upper equivalence bounds", "Pearson's correlation R2", "TOST Lower p-value", "TOST Upper p-value", "NHST t-test p-value","Equivalence")
      TOST.res <- TOST.res %>% mutate(across(!c(1:4,17), as.numeric))
      
      if(FDR) {
      TOST.res$`NHST t-test BH-adjusted p-value` <- p.adjust(TOST.res$`NHST t-test p-value`, method = "BH")
      TOST.res.no.NA <- TOST.res[rowSums(is.na(TOST.res[,-13])) == 0,, drop = F]
      TOST.res.no.NA[["Difference"]][(TOST.res.no.NA$`TOST Lower p-value` >= alpha | TOST.res.no.NA$`TOST Upper p-value` >= alpha) & TOST.res.no.NA$`NHST t-test BH-adjusted p-value` < alpha] <- "Yes"
      TOST.res <- merge(x = TOST.res, y = TOST.res.no.NA %>% dplyr::select(c("Gene name", "Difference")), by.x = "Gene name", by.y = "Gene name", all = T)
      TOST.res <- distinct(TOST.res)
      TOST.res$Difference[is.na(TOST.res$Difference)] <- "No"
      TOST.res <- TOST.res[,c(1:16,18,17,19)]
      } else {
        TOST.res.no.NA <- TOST.res[rowSums(is.na(TOST.res[,-13])) == 0,, drop = F]
        TOST.res.no.NA[["Difference"]][(TOST.res.no.NA$`TOST Lower p-value` >= 0.05 | TOST.res.no.NA$`TOST Upper p-value` >= 0.05) & TOST.res.no.NA$`NHST t-test p-value` < 0.05] <- "Yes"
        TOST.res <- merge(x = TOST.res, y = TOST.res.no.NA %>% dplyr::select(c("Gene name", "Difference")), by.x = "Gene name", by.y = "Gene name", all = T)
        TOST.res <- distinct(TOST.res)
        TOST.res$Difference[is.na(TOST.res$Difference)] <- "No"
      }
      
      TOST.res <- TOST.res[order(as.character(TOST.res[["Equivalence"]]), as.character(TOST.res[["Difference"]]), as.character(TOST.res[["Gene name"]]), decreasing = c(T,T,F), method = "radix"),]
      
      evalGenes.N <- if(gene.signature == "ALL_GENES") {ncol(normalized.counts.sel)} else {length(gene.list1.conv)}
      equivalent.genes.N <- sum(TOST.res$Equivalence == "Yes")
      non.equivalent.genes.N <- sum(TOST.res$Equivalence == "No")
      different.genes.N <- sum(TOST.res$Difference == "Yes")
      non.different.genes.N <- sum(TOST.res$Difference == "No")
      
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes evaluated to determine the gene signature:", evalGenes.N), startRow = 1)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes passing the initial DESeq2/edgeR filters:", length(all.tested.genes)), startRow = 2)      
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes passing the expression frequency filter (more than 10% of samples with more than 0 normalized counts in at least one group):", nrow(TOST.res)), startRow = 3)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes with non-equivalent expression profiles:", non.equivalent.genes.N), startRow = 5)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes with equivalent expression profiles:", equivalent.genes.N), startRow = 6)
      writeData(wb = wb, sheet = comparison.group, x = paste("The percentage of genes with equivalent expression profiles:", round(equivalent.genes.N/nrow(TOST.res)*100,2), "%"), startRow = 7)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes with non-different expression profiles:", non.different.genes.N), startRow = 9)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes with different expression profiles:", different.genes.N), startRow = 10)
      writeData(wb = wb, sheet = comparison.group, x = paste("The percentage of genes with different expression profiles:", round(different.genes.N/nrow(TOST.res)*100,2), "%"), startRow = 11)
      writeData(wb = wb, sheet = comparison.group, x = paste("The table below contains TOST and NHST statistical results for all the genes that passed the expression frequency filtering."), startRow = 13)

      writeData(wb = wb, sheet = comparison.group, x = TOST.res, startRow = 15, keepNA = TRUE, na.string = "NA")
      } else {
      evalGenes.N <- if(gene.signature == "ALL_GENES") {ncol(normalized.counts.sel)} else {length(gene.list1.conv)}
        
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes evaluated to determine the gene signature:", evalGenes.N), startRow = 1)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes passing the initial DESeq2/edgeR filters:", length(all.tested.genes)), startRow = 2)
      writeData(wb = wb, sheet = comparison.group, x = paste("The number of genes passing the expression frequency filter (more than 10% of samples with more than 0 normalized counts in at least one group):", 0), startRow = 3)
      }
  } else {
    writeData(wb = wb, sheet = comparison.group, x = paste0("The number of samples in one of the analyzed groups: ", min.N, " is lower than the predefined limit: ", min.N.limit, ". Therefore, the gene signature comparison has not been performed."))
  }
  } else {
    writeData(wb = wb, sheet = comparison.group, x = paste0("There are no samples in the following analyzed group(s): ", paste(missing.group, collapse = ", "), "."))
  } 
  }
  saveWorkbook(wb = wb, file = gene.signature.xlsx.name, overwrite = T)
}
}
}

if(exists("normalized.counts")) {
group.comparison(normalized.counts = normalized.counts)

if(txt.file != "NA") {
  group.comparison(normalized.counts = normalized.counts, txt.file = txt.file)
}
}

save.image(file = final.RData.name)
}
sessionInfo()
proc.time()
date()

cat("All done.\n")
