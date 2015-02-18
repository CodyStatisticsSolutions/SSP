#' Dataset info
#' 
#' Get some summary and descriptive statistics on a dataframe to show in a GUI.
#'  
#' @param data the dataframe
#' @return list with all kind of summary stuff
#' @examples spssfile <- system.file(package="StatSolutions", "files/1991GS.sav");
#' mydata <- read.spss.new(spssfile);
#' dataset.info(mydata);
#' @export
dataset.info <- function(data){
	
	mydf <- data;

	metadata <- list();
	metadata$nrow <- nrow(mydf);
	metadata$ncol <- ncol(mydf);
	
	allvars <- list();
	allvars$colnames <- colnames(mydf);
	allvars$class <- unname(sapply(sapply(mydf, class), head, 1));
	allvars$factor <- unname(sapply(lapply(lapply(mydf, class), "==", "factor"),any));
	allvars$ordered <- unname(sapply(lapply(lapply(mydf, class), "==", "ordered"),any));
	allvars$string <- unname(sapply(lapply(lapply(mydf, class), "==", "character"),any));
	allvars$numeric <- unname(sapply(lapply(lapply(mydf, class), "%in%", c("numeric", "integer")),any));
	allvars$missing <- unname(sapply(lapply(mydf, is.na), sum));
	
	factordf <- mydf[allvars$factor];
	allfactors <- list();
	if(length(factordf) > 0){
		for(i in 1:length(factordf)){
			allfactors[[i]] <- list();
			allfactors[[i]]$name <- names(factordf[i]);
			allfactors[[i]]$levels <- levels(factordf[[i]]);
			allfactors[[i]]$ordered <- is.ordered(factordf[[i]]);
			allfactors[[i]]$missing <- sum(is.na(factordf[[i]])); 
			allfactors[[i]]$frequencies <- as.list(table(factordf[[i]]));
			allfactors[[i]]$percentages <- as.list( table(factordf[[i]]) / length(factordf[[i]][!is.na(factordf[[i]])]) );
		}	
	}
	
	numericdf <- mydf[allvars$numeric];
	allnumerics <-list();	
	if(length(numericdf) > 0){
		for(i in 1:length(numericdf)){
			allnumerics[[i]] <- list();
			allnumerics[[i]]$name <- names(numericdf[i]);
			allnumerics[[i]]$max <- max(numericdf[[i]], na.rm=T);
			allnumerics[[i]]$min <- min(numericdf[[i]], na.rm=T);
			allnumerics[[i]]$mean <- mean(numericdf[[i]], na.rm=T);
			allnumerics[[i]]$median <- median(numericdf[[i]], na.rm=T);
			allnumerics[[i]]$stdev <- sd(numericdf[[i]], na.rm=T);
		}
	}

	return(list(metadata=metadata, variables=allvars, factorvariables=allfactors, numericvariables=allnumerics));
}
