\name{plotsamExplorer}
\alias{plotsamExplorer}
\title{
    Plots the results of output dataframe object.
}
\description{
  Boxplot results between sequence-depth and number of differentially 
  expressed genes.
}
\usage{
plotsamExplorer(dat, save = FALSE, filename = NULL, p.depth = 0.9, 
font.size = 3.5, anova = TRUE, x.lab=NULL, y.lab=NULL, leg.lab=NULL)
}
\arguments{
  \item{dat}{
    is a dataframe object, which consists three columns strictly labelled as:
     "Label", "Variable" and "Value". 
}
  \item{save}{
    is a logical value to save plot as a pdf. 
}
  \item{filename}{
     is a character to assign filename, if a user want to save the plot. 
}
  \item{p.depth}{
    is a numeric value for anova test to be performed for number 
differentially expressed genes of different sequence-depths.
}
  \item{font.size}{
    is a numeric value to set font size of the plot.
}
  \item{anova}{
    is a logical value for anova test to be performed for number 
differentially expressed genes of different sequence-depths.
}
  \item{x.lab}{
    is a string value to assign a label for x-axis.
}
  \item{y.lab}{
    is a string value to assign a label for y-axis.
}
  \item{leg.lab}{
    is a string vector assigns lables for legends in the plot.
}
}

\value{
    Generates a plot in a pdf format.
}
\author{
    Frank-Emmert Streib, Shailesh Tripathi, Aleksei sputnikov
}

\examples{

    data("df_sole")
    data("df_intersect")

    plotsamExplorer(df_sole,save=TRUE,filename="ss",p.depth=.9,
    font.size=4, anova=TRUE)
    plotsamExplorer(df_intersect,save=TRUE,filename="ss",p.depth=.9,
    font.size=4, anova=FALSE)


}
