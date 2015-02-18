#' Read CSV
#' 
#' Wrapper for read.csv
#' 
#' @param ... arguments passed on to read.csv
#' @return a data frame
#' @examples mytemp <- tempfile();
#' write.csv(iris, mytemp, row.names=FALSE);
#' newiris <- read.csv.new(mytemp);
#' stopifnot(identical(iris, newiris));
#' @export
read.csv.new <- function(...){
	mydata <- utils::read.csv(na.strings=c("","NA", " ", "#NULL!"),...);
  ssp.id.var <- as.numeric(rownames(mydata))
  mydata <- cbind(ssp.id.var, mydata)
  
	factors <- list();
	factors <- unname(sapply(lapply(lapply(mydata, class), "==", "factor"),any));  
	mydata[,c(factors)] <- as.data.frame(sapply(mydata[,c(factors)],gsub,pattern="\"",replacement="\'"))
  mydata[,c(factors)] <- as.data.frame(sapply(mydata[,c(factors)], gsub,pattern="“", replacement="\'"))
	mydata[,c(factors)] <- as.data.frame(sapply(mydata[,c(factors)], gsub,pattern="”", replacement="\'"))
  
  
  return(mydata)
}

#' Read.table
#' 
#' Wrapper for read.table
#' 
#' @param ... arguments passed to read.table 
#' @return data frame
#' @examples mytemp <- tempfile();
#' write.table(iris, mytemp, row.names=FALSE);
#' newiris <- read.table.new(mytemp, header=TRUE);
#' stopifnot(identical(iris, newiris));
#' @export
read.table.new <- function(...){
	mydata <-utils::read.table(na.strings=c("","NA", " ", "#NULL!"),...);
	ssp.id.var <- as.numeric(rownames(mydata))
	mydata <- cbind(ssp.id.var, mydata)
	mydata <- as.data.frame(sapply(mydata,gsub,pattern="\"",replacement="\'"))
	mydata <- as.data.frame(sapply(mydata,gsub,pattern="“",replacement="\'"))
	mydata <- as.data.frame(sapply(mydata,gsub,pattern="”",replacement="\'"))
	return(mydata)  
}