\name{plotsamExplorer}
\alias{plotsamExplorer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Plots the results of output dataframe object.
}
\description{
  Boxplot results between sequence-depth and number of differentially expressed genes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotsamExplorer(dat, save = FALSE, filename = NULL, p.depth = 0.9, font.size = 3.5, anova = TRUE, x.lab=NULL, y.lab=NULL, leg.lab=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    is a dataframe object, which consists three columns strictly labelled as:
     "Label", "Variable" and "Value". 
%%     ~~Describe \code{dat} here~~
}
  \item{save}{
    is a logical value to save plot as a pdf 
%%     ~~Describe \code{save} here~~
}
  \item{filename}{
     is a character to assign filename, if a user want to save the plot   
%%     ~~Describe \code{filename} here~~
}
  \item{p.depth}{
        is a numeric value for anova test to be performed for number differentially
         expressed genes of different sequence-depths.
%%     ~~Describe \code{p.depth} here~~
}
  \item{font.size}{
        is a numeric value to set font size of the plot.
%%     ~~Describe \code{font.size} here~~
}
  \item{anova}{
        is a logical value for anova test to be performed for number differentially
         expressed genes of different sequence-depths.
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
        Frank-Emmert Streib, Shailesh Tripathi, Aleksei sputnikov
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

    data("df_sole")
    data("df_intersect")

    plotsamExplorer(df_sole,save=TRUE,filename="ss",p.depth=.9,font.size=4, anova=TRUE)
    plotsamExplorer(df_intersect,save=TRUE,filename="ss",p.depth=.9,font.size=4, anova=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
