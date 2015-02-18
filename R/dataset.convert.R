#' Convert variables in data
#' 
#' Special function to convert variables in a dataset from nominal to numeric.
#' When used on datasets imported from SPSS, it uses orriginal codings. 
#' 
#' @param data dataset
#' @param asnumeric vector of variable names or indices to be converted to numbers
#' @param asfactor vector of variable names or indices to be converted to factors
#' @param asordered vector of variable names or indices to be converted to ordered factors
#' @return a new dataset
#' @export
#' 
dataset.convert <- function(data, asnumeric = vector(), asfactor = vector(), asordered = vector()){
	
	mydata <- data;
	
	#special function for datasets created by read.spss.savelabels
	if("spss.savelabels" %in% class(mydata)){
		
		nolabels <- attr(mydata, "nolabels");
		aslabels <- attr(mydata, "aslabels");
		
		#to convert to numeric, we use the orriginal spss codes.
		for(thisvar in asnumeric){
			mydata[[thisvar]] <- as.numeric(nolabels[[thisvar]]);
		}	
		
		#to convert variables to factors we use the factors as created by read.spss
		for(thisvar in asfactor){
			mydata[[thisvar]] <- as.factor(aslabels[[thisvar]]);
		}	
		
		#to convert variables to ordered factor we use the factors as created by read.spss
		for(thisvar in asordered){
			mydata[[thisvar]] <- as.ordered(aslabels[[thisvar]]);
		}			
		
	} else {
		#data created elsewise (read.csv, etc)
		
		for(thisvar in asnumeric){
			mydata[[thisvar]] <- as.numeric(mydata[[thisvar]]);
		}	
		
		for(thisvar in asfactor){
			mydata[[thisvar]] <- as.factor(mydata[[thisvar]]);
		}	
		
		for(thisvar in asordered){
			mydata[[thisvar]] <- as.ordered(mydata[[thisvar]]);
		}			
		
	}
	
	return(mydata);
}
