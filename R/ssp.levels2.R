ssp.levels2 <- function (data){
        
        names <- names(data)
        levels <- list()
        
        for(i in names){
                levels[[i]] <- levels(as.factor(data[[i]]))
        }
        
        return(levels)

}

