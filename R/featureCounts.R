samExplore <- function(files,annot.inbuilt="mm9",annot.ext=NULL,
isGTFAnnotationFile=FALSE,GTF.featureType="exon",
GTF.attrType="gene_id",useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,
minFragLength=50,maxFragLength=600,nthreads=1,strandSpecific=0,
minMQS=0,readExtension5=0,readExtension3=0,
read2pos=NULL,minReadOverlap=1,countSplitAlignmentsOnly=FALSE,
countMultiMappingReads=FALSE,countPrimaryAlignmentsOnly=FALSE,
countChimericFragments=TRUE,ignoreDup=FALSE,chrAliases=NULL,
reportReads=FALSE, subsample_d=1, N_boot=1){
 .printParam(subsample_d, N_boot) 
 resultlist <- list()
 for(i in 1:N_boot){
 cat("//========================== Sub sampling ========================\\\\ \n")
 cat(paste("||		Sub sample N: ", i,"			||","\n",sep=""))
 cat("\\\\===================================================================//")
 resultlist[[i]] <- .featureCounts(files,annot.inbuilt,annot.ext,
  isGTFAnnotationFile,GTF.featureType,GTF.attrType,useMetaFeatures,
  allowMultiOverlap,isPairedEnd,requireBothEndsMapped,checkFragLength,
  minFragLength,maxFragLength,nthreads,strandSpecific,minMQS,
  readExtension5,readExtension3,read2pos,minReadOverlap,
  countSplitAlignmentsOnly,countMultiMappingReads,
  countPrimaryAlignmentsOnly,countChimericFragments,ignoreDup,
  chrAliases,reportReads, subsample_d)	
 }
 resultlist
}


.featureCounts <- function(files,annot.inbuilt="mm9",annot.ext=NULL,
isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",
useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
isPairedEnd=FALSE,requireBothEndsMapped=FALSE,
checkFragLength=FALSE,minFragLength=50,
maxFragLength=600,nthreads=1,strandSpecific=0,
minMQS=0,readExtension5=0,readExtension3=0,
read2pos=NULL,minReadOverlap=1,countSplitAlignmentsOnly=FALSE,
countMultiMappingReads=FALSE,countPrimaryAlignmentsOnly=FALSE,
countChimericFragments=TRUE,ignoreDup=FALSE,chrAliases=NULL,
reportReads=FALSE, subsample_d=1)
{
flag <- FALSE
if(is.null(annot.ext)){
 switch(tolower(as.character(annot.inbuilt)),
 mm9={
 ann <- system.file("annot","mm9_RefSeq_exon.txt",package="samExploreR")
 cat("NCBI RefSeq annotation for mm9 (build 37.2) is used.\n")
 },
 mm10={
 ann <- system.file("annot","mm10_RefSeq_exon.txt",package="samExploreR")
 cat("NCBI RefSeq annotation for mm10 (build 38.1) is used.\n")
 },
 hg19={
 ann <- system.file("annot","hg19_RefSeq_exon.txt",package="samExploreR")
 cat("NCBI RefSeq annotation for hg19 (build 37.2) is used.\n")
 },
 {
 stop("In-built annotation for ", annot.inbuilt, " is not available.\n")
 }
 ) # end switch
}
else{
 if(is.character(annot.ext)){
 ann <- annot.ext
}
else{
 annot_df <- as.data.frame(annot.ext,stringsAsFactors=FALSE)
 if(sum(c("geneid","chr","start","end", "strand") %in% 
    tolower(colnames(annot_df))) != 5)
      stop("One or more required columns are missing in the provided 
      annotation data. Please refer to help page 
      for annotation format.\n")
 colnames(annot_df) <- tolower(colnames(annot_df))
 annot_df <- data.frame(geneid=annot_df$geneid,chr=annot_df$chr,
 start=annot_df$start,end=annot_df$end,strand=annot_df$strand,
 stringsAsFactors=FALSE)
 annot_df$chr <- as.character(annot_df$chr)
 fout_annot <- file.path(".",paste(".Rsubread_UserProvidedAnnotation_pid",
 Sys.getpid(),sep=""))
 oldScipen <- options(scipen=999)
 write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
 options(oldScipen)
 ann <- fout_annot
 flag <- TRUE
}
}

fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))
files_C <- paste(files,collapse=";")
#if(nchar(files_C) == 0) stop("No read files provided!")
 chrAliases_C <- chrAliases
 if(is.null(chrAliases))
 chrAliases_C <- " "
 countMultiMappingReads_C <- countMultiMappingReads
 if(countPrimaryAlignmentsOnly) countMultiMappingReads_C <- 2
 read2pos_C <- read2pos
 if(is.null(read2pos)) read2pos_C <- 0
 cmd <- paste("readSummary",ann,files_C,fout,as.numeric(isPairedEnd),
 minFragLength,maxFragLength,
 0,as.numeric(allowMultiOverlap),as.numeric(useMetaFeatures),nthreads,
 as.numeric(isGTFAnnotationFile),strandSpecific,as.numeric(reportReads),
 as.numeric(requireBothEndsMapped),as.numeric(!countChimericFragments),
 as.numeric(checkFragLength),GTF.featureType,GTF.attrType,minMQS,
 as.numeric(countMultiMappingReads_C),chrAliases_C," ",as.numeric(FALSE),
 14,readExtension5,readExtension3,minReadOverlap,
 as.numeric(countSplitAlignmentsOnly),
 read2pos_C," ",as.numeric(ignoreDup)," ", as.numeric(subsample_d),sep=",")
 n <- length(unlist(strsplit(cmd,",")))
 C_args <- .C("R_readSummary_wrapper",as.integer(n),as.character(cmd),
 PACKAGE="samExploreR")

 x <- read.delim(fout,stringsAsFactors=FALSE)
 colnames(x)[1:6] <- c("GeneID","Chr","Start","End","Strand","Length")
 x_summary <- read.delim(paste(fout,".summary",sep=""), stringsAsFactors=FALSE)
 file.remove(fout)
 file.remove(paste(fout,".summary",sep=""))
 if(flag) 
  file.remove(fout_annot)
  if(ncol(x) == 6){
   stop("No count data were generated.")
  }
 y <- as.matrix(x[,-c(1:6)])
 colnames(y) <- colnames(x)[-c(1:6)]
 rownames(y) <- x$GeneID
 z <- list(counts=y,annotation=x[,1:6],targets=colnames(y),stat=x_summary)
 z
}
.printParam <- function(depth=1, boot=1){
 cat("\n\n\n		samExplorer \n")
 cat("		designed to simulate reduces sequencing depth \n\n")
 cat("		with using of Rsubread \n\n\n")

 cat("//===================== samExplore setting =================\\\\ \n")
 cat("||							     ||\n")
 cat("||							     ||\n")
 cat("\\\\=========================================================//\n\n")


 cat("//========================== samExplore setting ============\\\\ \n")
 cat("||							     ||\n")
 cat("||							     ||\n")
 cat(paste("||	       Simulated Depth: ", depth*100,"%  ||", "\n",sep=""))
 cat(paste("||	       Number of sub samples: ", boot,"	  ||","\n",sep=""))
 cat("||							    || \n")
 cat("||					      	            || \n")
 cat("\\\\=========================================================//\n\n")
}
