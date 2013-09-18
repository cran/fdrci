\name{FDRplot}
\alias{FDRplot}

\title{
 Plot FDR Table Results
}
\description{
 This function plots FDR point and CI estimates over a sequence of possible significance thresholds. Results from fdrTbl() can be plotted directly as input to FDRplot.
}
\usage{
FDRplot(plotdat, lowerbound, upperbound, mn, lpos = "bottomleft", outfile = FALSE)
}

\arguments{
  \item{plotdat}{
  a table that is returned from fdrTbl(), or results formated in the same way.
}
  \item{lowerbound}{
  -log10(p-value) lower bound for the x-axis of the plot.
}
  \item{upperbound}{
  -log10(p-value) upper bound for the x-axis of the plot.
}
  \item{mn}{
  text for main title
}
  \item{lpos}{
  legend position, one of, "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
}
  \item{outfile}{
  Either 'FALSE' or the path/name of a file. If FALSE then the plot is displayed rather than saved to a file.
}
}

\references{
 Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
}
\author{
Joshua Millstein
}

\examples{
nrow_=100
ncol_=100
X = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
Y = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
nperm = 10

myanalysis = function(X,Y){
	ntests = ncol(X)
	rslts = as.data.frame(matrix(NA,nrow=ntests,ncol=2))
	names(rslts) = c("ID","pvalue")
	rslts[,"ID"] = 1:ntests
	for(i in 1:ntests){
		fit = cor.test(X[,i],Y[,i],na.action="na.exclude",
			alternative="two.sided",method="pearson")
		rslts[i,"pvalue"] = fit$p.value
	}
	return(rslts)
} # End myanalysis

## Generate observed results
obs = myanalysis(X,Y)

## Generate permuted results
perml = vector('list',nperm)
for(p_ in 1:nperm){
	X1 = X[order(runif(ncol_)),]
	perml[[p_]] = myanalysis(X1,Y)
}

## FDR results table
myfdrtbl = fdrTbl(obs$pvalue,perml,"pvalue",ncol_,0,3)
## Plot results
FDRplot(myfdrtbl,0,3,"My FDR Plot",lpos = "bottomleft")
}

\keyword{ htest }
\keyword{ nonparametric }
