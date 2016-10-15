samExplore <- function (files, annot.inbuilt = "mm10", annot.ext = NULL, isGTFAnnotationFile = FALSE, 
    GTF.featureType = "exon", GTF.attrType = "gene_id", chrAliases = NULL, 
    useMetaFeatures = TRUE, allowMultiOverlap = FALSE, minOverlap = 1, 
    largestOverlap = FALSE, readExtension5 = 0, readExtension3 = 0, 
    read2pos = NULL, countMultiMappingReads = FALSE, fraction = FALSE, 
    minMQS = 0, splitOnly = FALSE, nonSplitOnly = FALSE, primaryOnly = FALSE, 
    ignoreDup = FALSE, strandSpecific = 0, juncCounts = FALSE, 
    genome = NULL, isPairedEnd = FALSE, requireBothEndsMapped = FALSE, 
    checkFragLength = FALSE, minFragLength = 50, maxFragLength = 600, 
    countChimericFragments = TRUE, autosort = TRUE, nthreads = 1, 
    maxMOp = 10, reportReads = FALSE, subsample_d=1, N_boot=1,
    countboot=c("all","Assigned", "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",   "Unassigned_NoFeatures",
    "Unassigned_Unmapped",       "Unassigned_MappingQuality",
    "Unassigned_FragmentLength", "Unassigned_Chimera",
    "Unassigned_Secondary",      "Unassigned_Nonjunction",
    "Unassigned_Duplicate" )){

 .printParam(subsample_d, N_boot) 
 resultlist <- list()
 
fc <- featureCounts(
	files=files, annot.inbuilt =annot.inbuilt , annot.ext = annot.ext,
 isGTFAnnotationFile = isGTFAnnotationFile,
    GTF.featureType = GTF.featureType, GTF.attrType =GTF.attrType, chrAliases =chrAliases,
    useMetaFeatures = useMetaFeatures, allowMultiOverlap = allowMultiOverlap, minOverlap = minOverlap,
    largestOverlap = largestOverlap, readExtension5 = readExtension5, readExtension3 = readExtension3,
    read2pos = read2pos, countMultiMappingReads = countMultiMappingReads, fraction = fraction,
    minMQS = minMQS, splitOnly = splitOnly, nonSplitOnly = nonSplitOnly, primaryOnly = primaryOnly,
    ignoreDup = ignoreDup, strandSpecific = strandSpecific, juncCounts = juncCounts,
    genome = genome, isPairedEnd = isPairedEnd, requireBothEndsMapped = requireBothEndsMapped,
    checkFragLength = checkFragLength, minFragLength = minFragLength, maxFragLength = maxFragLength,
    countChimericFragments = TRUE, autosort = TRUE, nthreads = 1,
    maxMOp = maxMOp, reportReads=reportReads)
    
 ts <- c()
for(i in 1:ncol(fc$counts)){
  fctmp <- list()
  fctmp[["counts"]] <- as.matrix(fc$counts[,i])
  fctmp[["stat"]] <- fc$stat[,c(1,(i+1))]	
 cat("//========================== Sub sampling ========================\\\\ \n")
 cat(paste("||		Sub sample N: ", i,"			||","\n",sep=""))
 cat("\\\\===================================================================//")
 restmp <- .featureCountsBoot(fctmp,
  subsample_d=subsample_d, N_boot=N_boot, countboot=countboot)	
 
 resultlist[[i]] <- restmp[[1]]
 ts <- c(ts, restmp[[2]])
 }
names(resultlist) <- fc$targets
names(ts) <- fc$targets
list(bootres=resultlist, target.size=ts, featuremain=fc)
}


.featureCountsBoot <- function(flcount,
subsample_d=1, N_boot=10, countboot="all")
{
    ts <- seq <- NULL
	totalcount <- as.numeric(flcount$stat[,2])
	names(totalcount) <- as.vector(flcount$stat[,1])
	if(countboot[1]=="all"){
            ts <- round(sum(totalcount)*subsample_d)
	    seqtmp <- as.vector(flcount$counts)
	    names(seqtmp) <- rownames(flcount$counts)		
	    seq <- c(seqtmp, totalcount[-which(names(totalcount)=="Assigned")])	
        } 
	else{
	    countboot <- unique(c("Assigned",countboot))
            ts <- round(sum(totalcount[countboot])*subsample_d)
	    seqtmp <- as.vector(flcount$counts)
            names(seqtmp) <- rownames(flcount$counts)
	    tmpx <- totalcount[countboot]
	    tmpx <- tmpx[-which(names(tmpx)=="Assigned")]
	    seq <- c(seqtmp, tmpx)	
	}
	bootmat <- c()
   	
    for(i in 1:N_boot){	
        tmp <- thinCounts(seq, target.size=ts)
    	bootmat <- cbind(bootmat, tmp[rownames(flcount$counts),]) 
    }
    list(bootmat, target.size=unique(ts,sum(tmp)))

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
