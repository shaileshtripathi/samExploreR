\name{samExplore}
\alias{samExplore}
\title{samExplore: }
\description{ samExplore: This function assigns mapped sequencing reads to 
 genomic features and simulates a sample with reduced sequencing depth}
\usage{
samExplore(files,annot.inbuilt="mm9",annot.ext=NULL,isGTFAnnotationFile=FALSE,
GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,
allowMultiOverlap=FALSE,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,
checkFragLength=FALSE,minFragLength=50,maxFragLength=600,
nthreads=1,strandSpecific=0,minMQS=0,
readExtension5=0,readExtension3=0,read2pos=NULL,
minReadOverlap=1,countSplitAlignmentsOnly=FALSE,
countMultiMappingReads=FALSE,countPrimaryAlignmentsOnly=FALSE,
countChimericFragments=TRUE,ignoreDup=FALSE,chrAliases=NULL,
reportReads=FALSE, subsample_d=1, N_boot=1,
countboot=c("all","Assigned", "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",   "Unassigned_NoFeatures",
    "Unassigned_Unmapped",       "Unassigned_MappingQuality",
    "Unassigned_FragmentLength", "Unassigned_Chimera",
    "Unassigned_Secondary",      "Unassigned_Nonjunction",
    "Unassigned_Duplicate"
)
}
\arguments{
 \item{files}{ a character vector giving names of input files containing 
 read mapping results. The files can be in either SAM format or BAM format. 
 The file format is automatically detected by the function.}
 \item{annot.inbuilt}{ a character string specifying an in-built 
 annotation used for read summarization.
 It has three possible values including \code{mm9}, \code{mm10} and 
 \code{hg19}, corresponding to the NCBI RefSeq annotations for genomes 
 'mm9', 'mm10' and 'hg19', respectively. 
 \code{mm9} by default. The in-built annotation has a SAF format (see below).}
  \item{annot.ext}{ A character string giving name of a user-provided 
    annotation file or a data frame including user-provided annotation data. 
    If the annotation is in GTF format, it can only be provided as a file. 
    If it is in SAF format, it can be provided as a file or a data frame. 
    See below for more details about SAF format annotation. \code{annot.ext} 
    will override \code{annot.inbuilt} if they are both provided.}
  \item{isGTFAnnotationFile}{ logical indicating whether the annotation 
    provided via the \code{annot.ext} argument is in GTF format or not. 
    \code{FALSE} by default. This option is only applicable when 
    \code{annot.ext} is not \code{NULL}.}
  \item{GTF.featureType}{ a character string giving the feature type used 
    to select rows in the GTF annotation which will be used for read 
    summarization. \code{exon} by default. This argument is only applicable 
    when \code{isGTFAnnotationFile} is \code{TRUE}.}
  \item{GTF.attrType}{ a character string giving the attribute type in the 
    GTF annotation which will be used to group features (eg. exons) into 
    meta-features (eg. genes). \code{gene_id} by default. This argument is 
    only applicable when \code{isGTFAnnotationFile} is \code{TRUE}.}
  \item{useMetaFeatures}{ logical indicating whether the read summarization 
    should be performed at the feature level (eg. exons) or meta-feature level
    (eg genes). If \code{TRUE}, features in the annotation 
    (each row is a feature) will be grouped into meta-features 
    using their values in the "GeneID" column in the SAF-format annotation 
    file or using the "gene_id" attribute in the GTF-format annotation file, 
    and reads will assiged to the meta-features instead of the features. 
    See below for more details.}
  \item{allowMultiOverlap}{ logical indicating if a read is allowed to be 
    assigned to more than one feature (or meta-feature) if it is found to 
    overlap with more than one feature (or meta-feature). \code{FALSE} by 
    default.}
  \item{isPairedEnd}{ logical indicating if paired-end reads are used. 
    If \code{TRUE}, fragments (templates or read pairs) will be counted 
    instead of individual reads. \code{FALSE} by default.}
  \item{requireBothEndsMapped}{ logical indicating if both ends from 
    the same fragment are required to be successfully aligned before 
    the fragment can be assigned to a feature or meta-feature. 
    This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}.}
  \item{checkFragLength}{ logical indicating if the two ends from the same 
    fragment are required to satisify the fragment length criteria before 
    the fragment can be assigned to a feature or meta-feature. 
    This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}. 
    The fragment length criteria are specified via \code{minFragLength} 
    and \code{maxFragLength}.}
  \item{minFragLength}{ \code{integer} giving the minimum fragment length 
    for paired-end reads. \code{50} by default.}
  \item{maxFragLength}{ \code{integer} giving the maximum fragment length 
    for paired-end reads. \code{600} by default. \code{minFragLength} and 
    \code{maxFragLength} are only applicable when \code{isPairedEnd} is 
    \code{TRUE}. Note that when a fragment spans two or more exons, 
    the observed fragment length might be much bigger than the nominal 
    fragment length.}
  \item{nthreads}{ \code{integer} giving the number of threads used for 
    running this function. \code{1} by default.}
  \item{strandSpecific}{ \code{integer} indicating if strand-specific read 
    counting should be performed. It has three possible values: 
    \code{0} (unstranded), \code{1} (stranded) and 
    \code{2} (reversely stranded). \code{0} by default.}
  \item{minMQS}{ \code{integer} giving the minimum mapping quality score a 
    read must satisfy in order to be counted. For paired-end reads, at least 
    one end should satisfy this criteria. \code{0} by default.}
  \item{readExtension5}{ \code{integer} giving the number of bases extended 
    upstream from 5' end of each read. \code{0} by default.}
  \item{readExtension3}{ \code{integer} giving the number of bases extended 
    downstream from 3' end of each read. \code{0} by default.}
  \item{read2pos}{ Specifying whether each read should be reduced to its 
    5' most base or 3' most base. It has three possible values: \code{NULL}, 
    \code{5} (denoting 5' most base) and \code{3} (denoting 3' most base). 
    The default value is \code{NULL}. With the default value, the whole read 
    is used for summarization. When \code{read2pos} is set to \code{5} 
    (or \code{3}), read summarization will be performed based on the 5' (or 3') 
    most base position. \code{read2pos} can be used together with 
    \code{readExtension5} and \code{readExtension3} parameters to set any 
    desired length for reads.}
  \item{minReadOverlap}{ \code{integer} giving the minimum number of 
    overlapped bases between a read and a feature required for the read 
    to be assigned to the feature. Negative values are also accepted, 
    indicating a gap being allowed between a read and a feature. 
    \code{1} by default.}
  \item{countSplitAlignmentsOnly}{ logical indicating whether only split 
    alignments (their CIGAR strings containing letter 'N') should be used 
    for summarization. \code{FALSE} by default. Example split alignments 
    are exon-spanning reads from RNA-seq data. \code{useMetaFeatures} should 
    be set to \code{FALSE} and \code{allowMultiOverlap} should be set to 
    \code{TRUE}, if the purpose of summarization is to assign exon-spanning 
    reads to all their overlapping exons.}
  \item{countMultiMappingReads}{ logical indicating if multi-mapping 
    reads/fragments should be counted, \code{FALSE} by default. 
    If \code{TRUE}, a multi-mapping read will be counted up to 
    N times if it has N reported mapping locations. 
    This function uses the 'NH' tag to find multi-mapping reads.}
  \item{countPrimaryAlignmentsOnly}{ logical indicating if only primary 
  alignments should be counted. Primary and secondary alignments are 
  identified using bit 0x100 in the Flag field of SAM/BAM files. 
  If \code{TRUE}, all primary alignments in a dataset will be counted 
  no matter they are from multi-mapping reads or not 
  (ie. \code{countMultiMappingReads} is ignored).}
  \item{countChimericFragments}{ logical indicating whether a chimeric fragment, 
   which has its two reads mapped to different chromosomes, should be counted 
   or not. If this fragment overlaps with only one feature 
   (or one meta-feature), typically by one of its two read, this fragment will 
   be assigned to that feature (or meta-feature). 
   If it is found to overlap more than one feature (or meta-feature), 
   for example one of its two reads overlaps meta-feature A and the other 
   overlaps meta-feature B, and \code{allowMultiOverlap} is \code{FALSE}, 
   then this fragment will not be counted. This parameter is only appliable 
   when \code{isPairedEnd} is \code{TRUE}.}
  \item{ignoreDup}{ logical indicating whether reads marked as duplicates 
   should be ignored. \code{FALSE} by default. 
   Read duplicates are identified using bit Ox400 in the FLAG field 
   in SAM/BAM files. The whole fragment (read pair) will be ignored 
   if paired end.}
  \item{chrAliases}{ a character string giving the name of a file that 
   contains aliases of chromosome names. The file should be a comma 
   delimited text file that includes two columns. 
   The first column gives the chromosome names used in the annotation 
   and the second column gives the chromosome names used by reads. 
   This file should not contain header lines. Names included in 
   this file are case sensitive.}
  \item{reportReads}{ logical indicating if read counting result for 
   each read/fragment is saved to a file. If \code{TRUE}, 
   read counting results for reads/fragments will be saved to a 
   tab-delimited file that contains four columns including name of 
   read/fragment, status(assigned or the reason if not assigned), 
   name of target feature/meta-feature and number of hits if the read/fragment 
   is counted multiple times. Name of the file is the same as name of the input 
   read file except a suffix 'featureCounts' is added. Multiple files will be 
   generated if there is more than one input read file.}
 \item{subsample_d}{numeric value which describes fraction of reads to be 
  remained in subsampling. }
 \item{N_boot}{integer value for number of resample procedures to be run.}
\item{countboot}{is a character vector which contains following options:
    "all","Assigned", "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",   "Unassigned_NoFeatures",
    "Unassigned_Unmapped",       "Unassigned_MappingQuality",
    "Unassigned_FragmentLength", "Unassigned_Chimera",
    "Unassigned_Secondary",      "Unassigned_Nonjunction",
    "Unassigned_Duplicate  
A user can select any of theses options for resampling
if user selects "all" then the resampling procedure will consider all
"assigned" and "unassigned" reads. If a user selects "Assigned" option
then resampling procedure will consider "Assigned" reads only for resampling.
If a user selects any other option it will consider those unmapped reads along with
"Assigned" reads. A user can selects more than one choices and input as a vector 

}
  }
\value{returns a list object.}
\details{
\code{samExplore} See featureCounts for details.
Output is a list objects which has three components.
1) "bootres": is a list object of size of input files, each list object contains
a resampling matrix of features.

2) "target.size": it is a numeric vector contains total feature counts of a 
certain sequence depth for each input file.

3) "feature main": returns a list object which is the ouptput of 'featureCounts' function
of Rsubread package.

}
\examples{
## Not run:

# Simulate a sample with sequencing depth 80% of initial for SAM format 
# single-end reads using built-in RefSeq annotation for hg19:
#### Consider all mapped and unmapped reads for resampling##
inpf <- system.file("testdata", package="samExploreR")
inpf1 <- paste(inpf,"/","Test1.sam",sep="") 
inpf2 <- paste(inpf,"/","Test2.sam",sep="")
inpf <- c(inpf1, inpf2) 
res1 <- samExplore(files=inpf,annot.inbuilt="hg19", subsample_d = 0.8)

#### Consider Assigned and Unassigned Unmapped reads for resampling##
res2 <- samExplore(inpf, N_boot=10, subsample_d=.8,
countboot=c("Assigned","Unassigned_Unmapped"))


#### Consider only Assigned reads for resampling##
res3 <- samExplore(inpf, N_boot=10, subsample_d=.8,
countboot="Assigned")

}
