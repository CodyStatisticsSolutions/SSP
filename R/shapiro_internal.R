shapiro_internal <- function(x, y, data){
	splitdata <- split(data[[y]], data[[x]]);
	shtests <- lapply(splitdata, shapiro.test);
	datanames <- names(shtests);
	for(i in 1:length(datanames)){
		shtests[[i]]$group <- datanames[i]
		shtests[[i]]$data.name <- NULL;
	}
	names(shtests) <- NULL
	return(shtests)
}


#x <- "Species"
#y <- "Sepal.Width"
#data <- iris
#shapiro_internal(x,y, data)