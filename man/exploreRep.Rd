\name{exploreRep}
\alias{exploreRep}
\title{exploreRep: function to explore the reproducibility}
\description{This function explores the reproducibility of analysis with annotation altering}
\usage{
exploreRep(df_d, lbl_vect, f)
}
\arguments{
  \item{df_d}{ a dataframe containing the dataset to explore with 3 columns: 
   label, f ratio, value to compare (e.g. number of differentially 
   expressed genes)}
  \item{lbl_vect}{a vector of character strings specifing the labels for 
   which the analysis should be run}
  \item{f}{ A numeric value of f for which the analysis should be run}
  }
\details{
\code{exploreRep} function to explore the reproducibility of the analysis 
 with altering of annotation. It runs ANOVA test for values to compare 
 (e.g. number of differentially expressed genes) corresponding to different 
 Annotation labels (i.e. analysis' run for different annotation types)

This function takes as input a dataframe containing the dataset to explore.

Here is the example of the dataframe

\preformatted{
...
AnnotA	0.1	13
AnnotB	0.1	101	
AnnotC	0.1	36
AnnotA	0.1	13
AnnotB	0.1	101	
AnnotC	0.1	36
AnnotA	0.4	40
AnnotB	0.4	153	
AnnotC	0.4	62
AnnotA	0.8	71
AnnotB	0.8	203
AnnotC	0.8	160
...
}

\code{exploreRob}
Thired column gives the values to compare (here number of differentially
 expressed genes).

\code{exploreRep} function subsets the dataset to consider only valyes for one 
f and runs ANOVA test for groups corresponding to annotations of interest.
}
\value{ 
An output of aov function

}
\author{Alexey Stupnikov and Shailesh Tripathi}
%\note{}
\examples{
%\dontrun{
#library(samExploreR)
data("df_sole")
#run ANOVA for annotation types labeled 'New, Gene' and 'New, Exon' and 
#f value 0.9
exploreRep(df_sole, lbl_vect = c('New, Gene', 'Old, Gene'), f = 0.9)

#run ANOVA for annotation type labeled 'Old' and 'New' and f value 0.5
exploreRep(df_sole, lbl_vect = c('New, Gene', 'Old, Gene'), f = 0.5)
%}
}
