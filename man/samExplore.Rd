\name{samExplore}
\alias{samExplore}
\title{samExplore: }
\description{ samExplore: This function assigns mapped sequencing reads to 
 genomic features and simulates a sample with reduced sequencing depth}
\usage{
samExplore(..., subsample_d=1, N_boot=1,
    countboot=c("all","Assigned", "Unassigned_Ambiguity",
    "Unassigned_MultiMapping",   "Unassigned_NoFeatures",
    "Unassigned_Unmapped",       "Unassigned_MappingQuality",
    "Unassigned_FragmentLength", "Unassigned_Chimera",
    "Unassigned_Secondary",      "Unassigned_Nonjunction",
    "Unassigned_Duplicate" ))
}
\arguments{
 \item{...}{These are the same arguments of \code{featureCounts} function of 
 \code{Rsubread} package, for more details check \code{featureCounts} function.}
 \item{subsample_d}{numeric value which describes fraction of reads to be 
  remained in subsampling. }
 \item{N_boot}{integer value for number of resample procedures to be run.}
 \item{countboot}{is a character vector which contains following options:
 \code{all},\code{Assigned}, \code{Unassigned_Ambiguity},
    \code{Unassigned_MultiMapping},   \code{Unassigned_NoFeatures},
    \code{Unassigned_Unmapped},       \code{Unassigned_MappingQuality},
    \code{Unassigned_FragmentLength}, \code{Unassigned_Chimera},
    \code{Unassigned_Secondary},      \code{Unassigned_Nonjunction},
    \code{Unassigned_Duplicate}  
A user can select any of theses options for resampling
if user selects \code{all} then the resampling procedure will consider all
\code{assigned} and \code{unassigned} reads. If a user selects \code{Assigned} option
then resampling procedure will consider \code{Assigned} reads only for resampling.
If a user selects any other option it will consider those unmapped reads along with
\code{Assigned} reads. A user can selects more than one choices and input as a vector 
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
# Simulate a sample with sequencing depth 80% of initial for SAM format 
# single-end reads using built-in RefSeq annotation for hg19:
#### Consider all mapped and unmapped reads for resampling##
inpf <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
res1 <- samExplore(files=inpf,annot.inbuilt="hg19", subsample_d = 0.8)

#### Consider Assigned and Unassigned Unmapped reads for resampling##
res2 <- samExplore(inpf, N_boot=10, subsample_d=.8,
    countboot=c("Assigned","Unassigned_Unmapped"))

#### Consider only Assigned reads for resampling##
res3 <- samExplore(inpf, N_boot=10, subsample_d=.8,
    countboot="Assigned")
}
