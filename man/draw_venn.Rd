% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{draw_venn}
\alias{draw_venn}
\title{draw a venn plot}
\usage{
draw_venn(x, y, z, name)
}
\arguments{
\item{x, y, z}{a character vector ,diffriencial analysis genes by Deseq2 ,edgeR ,limma respectively}

\item{name}{main of the plot}
}
\value{
a venn plot according to \code{x}, \code{y} and.\code{z} named "name" paramter
}
\description{
print a venn plot for deg result created by three packages
}
\examples{
x = sample(1:100,30);y = sample(1:100,30);z = sample(1:100,30)
draw_venn(x,y,z,"test")
}
\seealso{
\code{\link{draw_pca}};\code{\link{draw_volcano}};\code{\link{draw_heatmap}}
}
\author{
Xiaojie Sun
}
