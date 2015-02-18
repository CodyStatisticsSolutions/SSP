#' Read SPSS (new)
#' 
#' Special wrapper for read.spss in order to preserve both orriginal encodings and factors.
#' 
#' @param file a valid spss data file (.sav)
#' @param ... arguments passed to read.spss
#' @param to.data.frame (depricated). Will be ignored. 
#' @aliases testdata
#' @return data frame
#' @import foreign
#' @examples spssfile <- system.file(package="StatSolutions", "files/1991GS.sav");
#' mydata <- read.spss.new(spssfile);
#' @export
read.spss.new <- function(file, ..., to.data.frame=FALSE){
	#read the dataset both with and without converting to factors
	aslabels <- foreign::read.spss(file=file, ..., to.data.frame=FALSE, use.missings=TRUE, use.value.labels=TRUE);
	nolabels <- foreign::read.spss(file=file, ..., to.data.frame=FALSE, use.missings=TRUE, use.value.labels=FALSE);
	
	#manually convert to dataframe because of stringsAsFactors
	aslabels <- as.data.frame(aslabels, stringsAsFactors=FALSE);
	nolabels <- as.data.frame(nolabels, stringsAsFactors=FALSE);
	
	#set the default
	mydata <- aslabels;
	
	#attach both 'orriginal' datasets.
	attr(mydata, "nolabels") <- nolabels;
	attr(mydata, "aslabels") <- aslabels;
	
	#return
	class(mydata) <- c("spss.savelabels", class(aslabels));
	ssp.id.var <- as.numeric(rownames(mydata));
	mydata <- cbind(ssp.id.var, mydata);  
	return(mydata);
}




