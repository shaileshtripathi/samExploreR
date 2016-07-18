### R code from vignette source 'Manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: preliminaries
###################################################
library(samExploreR)



###################################################
### code chunk number 3: Manual.Rnw:288-293
###################################################
#Loading library
library(samExploreR)
data("df_sole")
 #Performing robustness analysis
exploreRob(df_sole, lbl = 'New, Gene', f_vect = c(0.85, 0.9, 0.95))


###################################################
### code chunk number 4: Manual.Rnw:329-335
###################################################
#Loading library
library(samExploreR)
data("df_sole")
 #Performing robustness analysis
 t = exploreRep(df_sole, lbl_vect = c('New, Gene', 'Old, Gene', 'New, Exon'), f = 0.9)



###################################################
### code chunk number 5: plotunif
###################################################

require(samExploreR)
########## Loading the example data
data("df_sole")
data("df_intersect")
head(df_sole)
#head(df_intersect)



###################################################
### code chunk number 6: fig1
###################################################
### Generation of the plot:
require(samExploreR)
data("df_sole")
plotsamExplorer(df_sole,p.depth=.9,font.size=4, anova=TRUE)


###################################################
### code chunk number 7: fig2
###################################################
###Generation of the plot:
require(samExploreR)
data("df_intersect")
plotsamExplorer(df_intersect,font.size=3, anova=FALSE)


###################################################
### code chunk number 8: Manual.Rnw:401-402
###################################################
sessionInfo()


