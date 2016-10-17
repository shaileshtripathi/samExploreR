\name{exploreRob}
\alias{exploreRob}
\title{exploreRob: function to explore the robustness}
\description{This function explores the robustness of analysis with 
 sequencing depth altering}
\usage{
exploreRob(df_d, lbl, f_vect)
}
\arguments{
  \item{df_d}{ a dataframe containing the dataset to explore with 3 
  columns : label, f ratio, value to compare 
  (e.g. number of differentially expressed genes)}
  \item{lbl}{ a character string specifing the label 
   for which the analysis should be run}
  \item{f_vect}{ A numeric vector containing the values 
   of f for which the analysis should be run}
}
\details{
\code{exploreRob} function to explore the robustness of the analysis 
 with altering of sequencing depth. It runs ANOVA test for values to 
 compare (e.g. number of differentially expressed genes) corresponding 
 to different f ratio values (i.e. values of sequencing depth)

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

%exploreRob
%Third column gives the values to compare 
%(here number of differentially expressed genes)

 \code{exploreRob} function subsets the dataset to consider only valyes for 
 one type of annotation and runs ANOVA test for groups corresponding to f 
 values of interest.


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
#run ANOVA for annotation type labeled 'New, Gene' and f values 0.9, 0.95
exploreRob(df_sole, lbl = 'New, Gene', f_vect = c(0.9, 0.95))


#run ANOVA for annotation type labeled 'Old' and f values 0.5, 0.95
exploreRob(df_sole, lbl = 'Old, Gene', f_vect = c(0.5, 0.95))


%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

