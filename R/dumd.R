##' count unique values in every colunms for data.frame
##'
##' simplify pdata, delete columns with unique values,which can't be used as group vector
##' @usage dumd(x = iris)
##' @param x A data.frame.
##' @return The simple data.frame of columns unique values count in \code{x}
##' @examples
##' dumd(iris)
##' data(ToothGrowth)
##' x = ToothGrowth
##' dumd(ToothGrowth)
##' @section just :
##' See what are you doing

dumd <- function(x){
  suppressMessages(library(dplyr))
  colname <- vector("character")
  count <- vector("integer")
  for(i in 1:ncol(x)){
    colname[i] = colnames(x)[[i]]
    count[i]=nrow(x[!duplicated(x[,i]),])
  }
  suppressMessages(library(tibble))
  df <- tibble(colname,count) %>%
    arrange(desc(count))
  print(df)
}
