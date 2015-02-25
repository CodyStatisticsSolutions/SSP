# ggplot2 of various statistical graphs
# required libraries - ggplot2, GGally, 

#' @param df - data fram that contains the data to be plotted
#' @param plottye - type of plots available.  Currently,
#'        Barplot, Histogram, Boxplot, Scatter Plot, line plot (CI for time series forecast), violin plot, density plot #
#'        Plot Matrix (MXN variables), Parallel Coordinates Plot (Multiple Continuous Variables) 
#' @param vpanels - variables defining vertical panels, Default to NULL, factor or character, can be multiple variables, panel label automatically generated from the variable values
#' @param hpanels - variables defining horizontal panels, Default to NULL, factor or character, can be multiple variables, panel label automatically generated from the variable values
#' @param freey - y axis range in panel plots, default to fixed (every panel has the same y-range), can be free_y (variable range)
#' @param plottitle - Plot Title, Default "Quick Plot"
#' @param titlefsize - Plot Title font size, Default to 9
#' @param xval - variable used for x-axis
#' @param yval - variable used for y-axis
#' @param varlist - only applicable to matrix plot and parallel plot
#' @param histdensity - Plot percentage/density instead of raw count in barplot and histogram, Default to TRUE
#' @param linegroup - Create multiple line (x-axis is time dimension and y is any continuous variable) by group in line plots, Default to NULL
#' @param vref - Add vertical reference line(s) at x = a0, a1, ..... Default to x = 0 (y-axis)
#' @param href - Add horizontal reference line(s) at y = b0, b1, ....., Default to y = 0 (x-axis)
#' @param xlabel - X-axis Label
#' @param ylabel - Y-axis Label
#' @param xlimit - X-axis Range, list of two numbers
#' @param ylimit - y-axis Range, list of two numbers
#' @param ylogsqrt - Take log or square root of the Y-axis, Default to "none"
#' @param grpcolour - Variable that defines the color of subgroup in plots, factor/character variable except density plot (can be expanded to numberic in the future)
#' @param grpsize - Variable that defines the Size of subgroups in plots, factor/character/numeric variable, mainly for use in scatterplot 
#' @param pointsize - fixed size of observation dots in scatterplots or thickness in lineplots, Default to 2
#' @param plotpos - specifically for bar plots and density plots, stacking/filling/side-by-side, Default is stack
#' @param transparency - Transparency of plots, range from 0 (completely transparent) to 1 (solid, Default), 
#'                       useful for scatterplots and parallel coordinate plots where overplotting can occur
#' @param colpal - Colour palette used for grpcolour, Default to "Set1" for distinct groups
#' @param bwbackground - black and white background - Default to TRUE
#' @param flipaxis - transpose x and y axis - Default to FALSE
#' @param add.density - Overlay density estimation on top of histogram - Default to FALSE
#' @param add.smooth - Overlay smoothing average on top of scatterplots - Default to FALSE     
#' @param add.points - Overlay jittered points on top of boxplots - Default to FALSE  
#' @param add.normal - Add standard normal density to histogram or density plot - Default to FALSE

#library(ggplot2)
#library(GGally)
#library(plyr)

quickplot.ssp <- function(df, 
                          plottype = c("barplot", "histogram", "boxplot", "scatterplot", "lineplot", "violinplot", "densityplot", "matrixplot", "parallelplot"), 
                          vpanels = list(), hpanels = list(), freey = c("fixed", "free_y"), 
                          plottitle = "Quick Plot", titlefsize = 16,
                          xval, yval, varlist = list(), histdensity = TRUE, linegroup = NULL, 
                          vref = list(), href = list(), xlabel = "", ylabel = "", xlimit = list(), ylimit = list(), ylogsqrt = c("none","log","sqrt"),
                          grpcolour = NULL, grpsize = NULL, pointsize = 1, 
                          plotpos = c("dodge", "stack", "fill"), transparency = 1, colpal = "Set1", 
                          bwbackground = TRUE, flipaxis = FALSE, 
                          add.density = FALSE, add.smooth = FALSE, add.points = FALSE, add.normal = FALSE) {
        
        # Bar Plot #  
        if (plottype == "barplot") {
                xvalstr <- deparse(substitute(xval))
                fillstr <- deparse(substitute(grpcolour)) 
                if (fillstr == "NULL") {fillstr <- NULL}
                
                # Create plot layer #
                if (histdensity == TRUE) {
                        p <- ggplot(data = df, aes_string(x = xvalstr, fill = fillstr)) + 
                                geom_bar(aes(y = (..count..)/sum(..count..)), position=plotpos[1]) + 
                                ylab("Relative Frequency")
                } else {
                        p <- ggplot(data = df, aes_string(x = xvalstr, fill = fillstr)) + geom_bar(position = plotpos[1])
                } 
                if (plotpos[1] == "fill") {
                        p <- p + ylab("Cumulative Percent by X-value (x 100)")
                }
        } 
        
        # Histogram - overlaying with density plot is best used with "dodge" position #
        else 
                if (plottype == "histogram") {
                        xvalstr <- deparse(substitute(xval))
                        fillstr <- deparse(substitute(grpcolour))
                        if (fillstr == "NULL") {fillstr <- NULL}
                        
                        # Create plot layer #
                        if (histdensity == TRUE) {
                                p <- ggplot(data = df) + aes_string(x = xvalstr, fill = fillstr) + 
                                        geom_histogram(aes(y = ..density..), position = plotpos[1])
                                if (add.density == TRUE) {
                                        p <- p + geom_density(aes_string(colour = fillstr), alpha=0, size = pointsize) 
                                } 
                                if (add.normal == TRUE) {
                                        m <- median(df[,xvalstr], na.rm = TRUE) 
                                        p <- p + stat_function(fun = dnorm, colour = "black", size = pointsize, 
                                                               args = list(mean = m)) 
                                }
                        } else {
                                p <- ggplot(data = df) + aes_string(x=xvalstr, fill=fillstr) + geom_histogram(position = plotpos[1])
                        } 
                } 
        
        # Box Plot #
        else 
                if (plottype == "boxplot") {
                        xvalstr <- deparse(substitute(xval))
                        yvalstr <- deparse(substitute(yval))       
                        fillstr <- deparse(substitute(grpcolour))
                        if (fillstr == "NULL") {fillstr <- NULL}
                        
                        # Create plot layer #
                        p <- ggplot(data = df) + aes_string(x = xvalstr, y = yvalstr, fill = fillstr) + geom_boxplot(outlier.colour = "black", outlier.size = 2)      
                        if (add.points == TRUE) {
                                p <- p + geom_jitter() 
                        }    
                }
        
        # SCatter Plot #
        else 
                if (plottype == "scatterplot") {
                        xvalstr <- deparse(substitute(xval))
                        yvalstr <- deparse(substitute(yval))       
                        colourstr <- deparse(substitute(grpcolour))
                        if (colourstr == "NULL") {colourstr <- NULL}
                        sizestr <- deparse(substitute(grpsize))
                        if (sizestr == "NULL") {sizestr <- NULL}
                        
                        # Create plot layer #
                        if (is.null(sizestr)) {
                                p <- ggplot(data = df) + aes_string(x = xvalstr, y = yvalstr, colour = colourstr) + geom_point(alpha = transparency, size = pointsize)
                        } else {
                                p <- ggplot(data = df) + aes_string(x=xvalstr, y = yvalstr, colour = colourstr, size = sizestr) + geom_point(alpha = transparency)         
                        }
                        # Add Smoother - Default is loess (< 1000 obs) or gam (>= 1000 obs) #
                        if (add.smooth == TRUE) {
                                p <- p + stat_smooth(se = FALSE, size = 1)
                        }
                }
        
        # Line Plot #
        else
                if (plottype == "lineplot") {
                        xvalstr <- deparse(substitute(xval))
                        yvalstr <- deparse(substitute(yval))       
                        sizestr <- deparse(substitute(grpsize))
                        if (sizestr == "NULL") {sizestr <- NULL}       
                        groupstr <- deparse(substitute(linegroup))
                        if (groupstr == "NULL") {groupstr <- NULL}
                        
                        # Create plot layer #
                        if (is.null(sizestr)) {
                                p <- ggplot(data = df) + aes_string(x = xvalstr, y = yvalstr, colour = groupstr, group = groupstr) + 
                                        geom_line(alpha=transparency, size = pointsize) 
                        } else {
                                p <- ggplot(data = df) + aes_string(x = xvalstr, y = yvalstr, colour = groupstr, group = groupstr, size = sizestr) + 
                                        geom_line(alpha=transparency)          
                        }
                        if (!is.null(groupstr)) {
                                if (is.numeric(df[,groupstr])) {
                                        p <- p + scale_colour_gradient(low="red") 
                                }          
                        }
                }
        
        # Violin Plot #
        else 
                if (plottype == "violinplot") {
                        xvalstr <- deparse(substitute(xval))
                        yvalstr <- deparse(substitute(yval))       
                        fillstr <- deparse(substitute(grpcolour))
                        if (fillstr == "NULL") {fillstr <- NULL}
                        
                        # Create plot layer #      
                        p <- ggplot(data = df) + aes_string(x = xvalstr, y = yvalstr, fill = fillstr) +  
                                geom_violin(scale = "count") + geom_boxplot(width = 0.1) 
                        if (add.points == TRUE) {
                                p <- p + geom_jitter() 
                        }  
                }
        
        # Density Plot #
        else
                if (plottype == "densityplot") {
                        xvalstr <- deparse(substitute(xval))  
                        groupstr <- deparse(substitute(linegroup))
                        if (groupstr == "NULL") {groupstr <- NULL}
                        
                        # Create plot layer #
                        p <- ggplot(data = df) + aes_string(x = xvalstr, colour = groupstr, group = groupstr) + 
                                geom_density(size = pointsize)
                }
        
        # Matrix Plot #
        else   
                if (plottype == "matrixplot") {
                        varlistp <- paste(substitute(varlist))[-1]
                        grpcolourColumn <- deparse(substitute(grpcolour))
                        if (grpcolourColumn == "NULL") {grpcolourColumn <- NULL} 
                        
                        # Create plot layer #
                        p <- ggpairs(df[,varlistp], 
                                     color = grpcolourColumn, 
                                     legends = TRUE, 
                                     upper=list(params=list(corSize=1)),
                                     diag=list(continuous="density",discrete="bar"), 
                                     axisLabels = "show",
                                     title =  plottitle) + 
                                theme(strip.text.x = element_text(face="bold",size=8),
                                      strip.text.y = element_text(face="bold",size=8))
                }
        
        # Parallel Coordinate Plot #
        else 
                if (plottype == "parallelplot") {
                        varlistp <- paste(substitute(varlist))[-1]
                        grpcolourColumn <- deparse(substitute(grpcolour))
                        if (grpcolourColumn == "NULL") 
                        {grpcolourColumn <- NULL} else 
                        {grpcolourColumn <- which(names(df) == grpcolourColumn)}
                        
                        # Create plot layer #   
                        # Default Variable Scaling - Standard Normal Scores (can be expanded in the future) # 
                        p <- ggparcoord(data = df, 
                                        columns = which(names(df) %in% varlistp),
                                        groupColumn = grpcolourColumn,
                                        alphaLines = transparency,
                                        title = plottitle)
                        if (!xlabel == "") {p <- p + xlab(xlabel)}
                        if (!ylabel == "") {p <- p + ylab(ylabel)}       
                }
        
        # Add black and white background #  
        if (bwbackground == TRUE) {
                p <- p + theme_bw()
        } 
        
        # Add themes (except matrixplots and parallel coordinate plots) #
        if (!plottype %in% c("matrixplot", "parallelplot")) {  
                
                # Add vertical and/or horizontal reference line(s) if applicable #    
                vreflines <- as.numeric(paste(substitute(vref))[-1])
                hreflines <- as.numeric(paste(substitute(href))[-1])
                if (length(vreflines) > 0) {
                        vreflines.df <- data.frame(vreflines=vreflines)
                        p <- p + geom_vline(aes(xintercept=vreflines), data = vreflines.df) 
                }
                if (length(hreflines) > 0) {
                        hreflines.df <- data.frame(hreflines=hreflines)
                        p <- p + geom_hline(aes(yintercept=hreflines), data = hreflines.df) 
                }
                
                # Split plots into panels by vpanels and hpanels #
                vpanelstr <- paste(substitute(vpanels))[-1]
                hpanelstr <- paste(substitute(hpanels))[-1] 
                if (length(vpanelstr)>0 | length(hpanelstr)>0) {
                        if (length(vpanelstr)==0) 
                        {vpanelstr <- "."} else {
                                vpanelstr <- paste(vpanelstr,collapse="+")
                        }
                        if (length(hpanelstr)==0) 
                        {hpanelstr <- "."} else {
                                hpanelstr <- paste(hpanelstr,collapse="+")       
                        }
                        p <- p + facet_grid(as.formula(paste(hpanelstr,"~", vpanelstr)), scales = freey[1])
                        plottitle <- paste(plottitle, "by", vpanelstr, "(Vertical) and", hpanelstr, "(Horizontal)")
                }
                
                # Modify x- and y-axis labels #
                if (!xlabel == "") {
                        p <- p + xlab(xlabel) + theme(axis.title.x = element_text(vjust = 0))
                }
                
                if (!ylabel == "") {
                        p <- p + ylab(ylabel) + theme(axis.title.y = element_text(vjust = 1))
                }
                
                # y in log scale #
                if (ylogsqrt[1] == "log") {
                        p <- p + scale_y_log10()
                } else if (ylogsqrt[1] == "sqrt") {
                        p <- p + scale_y_sqrt()    
                }
                
                p <- p + 
                        ggtitle(plottitle) + 
                        theme(axis.title.x = element_text(face="bold",size=14),
                              axis.title.y = element_text(face="bold",size=14),
                              axis.text.x = element_text(size=10),
                              axis.text.y = element_text(size=10),
                              plot.title = element_text(face="bold",size=titlefsize),
                              legend.title = element_text(size=14),
                              legend.text = element_text(size=14),
                              strip.text.x = element_text(face="bold",size=11),
                              strip.text.y = element_text(face="bold",size=11)
                        )
                
                # Flip axes #
                if (flipaxis == TRUE) {
                        p <- p + coord_flip()
                }
        }
        
        print(p)
        
} 


########################################################################################
#################################### End of Scripts ####################################
########################################################################################


data <- read.csv(file.choose())

#quickplot.ssp(df=data, plottype="barplot", xval=Group.3)
#quickplot.ssp(df=data, plottype="barplot", xval=Group.3, plottitle="Barplot Example", xlabel="Group 3", ylabel="Frequency")
#quickplot.ssp(df=data, plottype="barplot", xval=Group.3, plottitle="Barplot Example", xlabel="Group 3", ylabel="Frequency", grpcolour = Different.Group.3, colpal="Set2")
#quickplot.ssp(df=data, plottype="barplot", xval=Group.3, plottitle="Barplot Example", xlabel="Group 3", ylabel="Percentage", grpcolour = Different.Group.3, colpal="Set2", histdensity=TRUE, add.density=TRUE, vpanels=list(Different.Group.2), hpanels=list(Chi.Group.2.Large), href=list(0.10), flipaxis=TRUE)


#quickplot.ssp(df=data, plottype="histogram", xval=Normal)
#quickplot.ssp(df=data, plottype="histogram", xval=Normal, plottitle="Histogram Example", xlabel="Normal", ylabel="Percentage", histdensity=TRUE, href=list(0.10), flipaxis=TRUE, grpcolour=Different.Group.3)

#quickplot.ssp(df=data, plottype="boxplot", xval=0, yval=Not.Normal.Dependent)
#quickplot.ssp(df=data, plottype="boxplot", xval=Group.3, yval=Not.Normal.Dependent, plottitle="Boxplot Example", xlabel="Group 3", ylabel="Not Normal Dependent", flipaxis=TRUE, grpcolour=Different.Group.3, add.points=TRUE)

#quickplot.ssp(df=data, plottype="scatterplot", xval=Normal, yval=Not.Normal.Dependent)
#quickplot.ssp(df=data, plottype="scatterplot", xval=Normal, yval=Not.Normal.Dependent, plottitle="Scatterplot Example", xlabel="Normal", ylabel="Not Normal Dependent", flipaxis=TRUE, grpcolour=Different.Group.3, add.smooth=TRUE)


#quickplot.ssp(df=data, plottype="lineplot", xval=Normal, yval=Not.Normal.Dependent)
#quickplot.ssp(df=data, plottype="lineplot", xval=Normal, yval=Not.Normal.Dependent, plottitle="Lineplot Example", xlabel="Normal", ylabel="Not Normal Dependent", flipaxis=TRUE, linegroup=Different.Group.3)


#quickplot.ssp(df=data, plottype="violinplot", xval=Chi.Group.2.Large, yval=Not.Normal.Dependent)
#quickplot.ssp(df=data, plottype="violinplot", xval=Chi.Group.2.Large, yval=Not.Normal.Dependent, plottitle="Violinplot Example", xlabel="Chi Group 2 Large", ylabel="Not Normal Dependent", flipaxis=TRUE, grpcolour=Group.3)


#quickplot.ssp(df=data, plottype="densityplot", xval=Normal)
#quickplot.ssp(df=data, plottype="densityplot", xval=Normal, plottitle="Density Plot Example", xlabel="Normal", ylabel="Percentage", flipaxis=TRUE, linegroup=Group.3, pointsize=7)

#quickplot.ssp(df=data, plottype="matrixplot", varlist=list(Normal, Group.3, Not.Normal.Dependent, Chi.Group.2.Large))
#quickplot.ssp(df=data, plottype="matrixplot", varlist=list(Normal, Group.3, Not.Normal.Dependent, Chi.Group.2.Large), plottitle="Matrix Plot Example", grpcolour=Group.3)

#quickplot.ssp(df=data, plottype="parallelplot", varlist=list(Normal, Not.Normal.Dependent, Normal.Dependent))
#quickplot.ssp(df=data, plottype="parallelplot", varlist=list(Normal, Not.Normal.Dependent, Normal.Dependent), plottitle="Parallel Plot Example", grpcolour=Group.3, xlabel="Variables", ylabel="Standardized Score")




# Testing using R example data #
#mtcars[,"gear2"] <- factor(mtcars[,"gear"])
#mtcars[,"cyl2"] <- factor(mtcars[,"cyl"])
#diamonds.samp <- diamonds[sample(1:dim(diamonds)[1],200),]
#movies[,"Drama2"] <- factor(movies[,"Drama"])
#movies[,"decade"] <- round_any(movies$year, 10)
#mry <- do.call(rbind, by(movies, round(movies$rating), function(df) {
#  nums <- tapply(df$length, df$year, length)
#  data.frame(rating=round(df$rating[1]), year = as.numeric(names(nums)), number=as.vector(nums))
#}))
#mry[,"rating2"] <- factor(mry[,"rating"])#


#quickplot.ssp(df = mtcars, plottype="matrixplot", varlist = list(mpg,cyl,disp))
#quickplot.ssp(df = mtcars, plottype="matrixplot", varlist = list(mpg,cyl,disp,gear2), grpcolour=gear2)
#quickplot.ssp(df = diamonds.samp, plottype="matrixplot", varlist = list(carat, cut, color, clarity, depth), grpcolour=cut)

#quickplot.ssp(df = mtcars, plottype="parallelplot", varlist = list(mpg,cyl,disp), grpcolour=NULL, transparency = 0.3)
#quickplot.ssp(df = mtcars, plottype="parallelplot", varlist = list(mpg,cyl,disp), grpcolour=gear2, transparency = 0.3)

#quickplot.ssp(df = diamonds.samp, plottype="barplot", xval = cut, href = list(0.1,0.3,0.5))
#quickplot.ssp(df = diamonds.samp, plottype="barplot", xval = carat, grpcolour = cut)
#quickplot.ssp(df = diamonds.samp, plottype="barplot", xval = carat, grpcolour = cut, plotpos = "fill", histdensity=FALSE)
#quickplot.ssp(df = diamonds.samp, plottype="barplot", xval = carat, grpcolour = cut, plotpos = "dodge", histdensity=FALSE, flipaxis = TRUE)

#quickplot.ssp(df = movies, plottype="histogram", xval = rating)
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=TRUE)
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=TRUE, plotpos="fill", plottitle = "Test PLotting")
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=TRUE, grpcolour=Drama2, plotpos="fill")
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=FALSE, grpcolour=Drama2, plotpos="fill")
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=FALSE, grpcolour=Drama2, plotpos="dodge", pointsize=2)
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=TRUE, grpcolour=Drama2, plotpos="dodge", pointsize=2)
#quickplot.ssp(df = movies, plottype="histogram", xval = rating, add.density=TRUE, grpcolour=Drama2, plotpos="dodge", pointsize=2, vref = list(5,10))

#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price)
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = cut)
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = cut, add.points = TRUE)
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = cut, add.points = TRUE, flipaxis = TRUE)
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = clarity, hpanels = list(cut, color))
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = clarity, hpanels = list(cut), vpanels=(color))
#quickplot.ssp(df = diamonds.samp, plottype="boxplot", xval = cut, yval = price, grpcolour = clarity, hpanels = list(cut), vpanels=(color), flipaxis = TRUE)

#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, pointsize=3, add.smooth = TRUE)
#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, pointsize = 3)
#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, grpsize = cut, add.smooth = TRUE)
#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, grpsize = carat)
#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, grpsize = carat, 
#              add.smooth = FALSE, href = list(55,60,65), vref = list(5000), vpanels = NULL, hpanels = list(clarity), transparency=0.4)
#quickplot.ssp(df = diamonds.samp, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, grpsize = carat, 
#              add.smooth = TRUE, href = list(55,60,65), vref = list(5000), vpanels = NULL, hpanels = list(clarity), transparency=0.4)
#quickplot.ssp(df = diamonds, plottype="scatterplot", xval = price, yval = depth, grpsize = carat, 
#              add.smooth = TRUE, href = list(55,60,65), vref = list(5000), vpanels = NULL, hpanels = list(clarity), transparency=0.5)
#quickplot.ssp(df = diamonds, plottype="scatterplot", xval = price, yval = depth, grpcolour = color, grpsize = carat, 
#              add.smooth = TRUE, href = list(55,60,65), vref = list(5000), vpanels = NULL, hpanels = list(clarity), transparency=0.85, freey = "free_y")

#quickplot.ssp(df = mry, plottype="lineplot", xval = year, yval = number, ylogsqrt = "log")
#quickplot.ssp(df = mry, plottype="lineplot", xval = year, yval = number, linegroup = rating, ylogsqrt = "log")
#quickplot.ssp(df = mry, plottype="lineplot", xval = year, yval = number, linegroup = rating, grpsize = rating)
#quickplot.ssp(df = mry, plottype="lineplot", xval = year, yval = number, linegroup = rating2)

#quickplot.ssp(df = mtcars, plottype="violinplot", xval = cyl2, yval = mpg)
#quickplot.ssp(df = mtcars, plottype="violinplot", xval = cyl2, yval = mpg, grpcolour = cyl2)
#quickplot.ssp(df = mtcars, plottype="violinplot", xval = cyl2, yval = mpg, grpcolour = cyl2, add.points = TRUE)

#quickplot.ssp(df = diamonds.samp, plottype="violinplot", xval = cut, yval = price, hpanels = list(clarity))
#quickplot.ssp(df = diamonds.samp, plottype="violinplot", xval = cut, yval = price, hpanels = list(clarity), flipaxis = TRUE)
#quickplot.ssp(df = diamonds.samp, plottype="violinplot", xval = cut, yval = price, grpcolour = clarity, hpanels = list(clarity), freey = "free_y")
#quickplot.ssp(df = subset(diamonds.samp, clarity == "SI1" & cut == "Ideal"), plottype="violinplot", 
#              xval = cut, yval = price, grpcolour = clarity, hpanels = list(clarity))
#quickplot.ssp(df = subset(diamonds.samp, clarity == "VS1" & cut == "Ideal"), plottype="violinplot", xval = cut, yval = price, grpcolour = clarity)
#quickplot.ssp(df = diamonds, plottype="violinplot", xval = cut, yval = price, grpcolour = clarity, hpanels = list(clarity))

#quickplot.ssp(df = movies, plottype="densityplot", xval = rating)
#quickplot.ssp(df = movies, plottype="densityplot", xval = rating, grpcolour = decade, linegroup = decade, pointsize = 2)

################
# Use testdata #
################

#testdata <- read.csv("testdata.csv")
#testdata$normalx <- rnorm(150)
#testdata2 <- data.frame(x=rnorm(1e5), y=factor(rep(letters[1:5], each=1e4)))

#quickplot.ssp(df = testdata, plottype="matrixplot", varlist = list(NotNormalDependent, NotDifferentGroup3, BigGroup, Paired5), bwbackground = FALSE)
#quickplot.ssp(df = testdata, plottype="matrixplot", bwbackground = TRUE, 
#              varlist = list(NotNormalDependent, NotDifferentGroup3, BigGroup, Paired5), grpcolour=BigGroup, href=list(20))
#quickplot.ssp(df = testdata, plottype="parallelplot", bwbackground = TRUE, 
 #             varlist = list(NotNormalDependent,Paired5,Outliers2, IVModerate, Paired6), grpcolour=NULL, transparency = 0.1)

#quickplot.ssp(df = testdata, plottype="barplot", xval = DifferentGroup2, grpcolour = DifferentGroup3, 
#              plotpos = "fill", histdensity=TRUE, href = list(0.1,0.3,0.5))
#quickplot.ssp(df = testdata, plottype="barplot", xval = DifferentGroup2, grpcolour = DifferentGroup3, 
#              plotpos = "dodge", histdensity=TRUE, href = list(0.1,0.3,0.5), flipaxis = TRUE)


#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=FALSE, plotpos="dodge", pointsize=2, flipaxis = TRUE)
#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=FALSE, plotpos="dodge", pointsize=2, histdensity = FALSE)
#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=TRUE, plotpos="dodge", pointsize=2, histdensity = FALSE)
#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=FALSE, plotpos="dodge", pointsize=2, histdensity = TRUE)
#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=TRUE, grpcolour=DifferentGroup3, plotpos="dodge", pointsize=2)

#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=FALSE, grpcolour=DifferentGroup3, plotpos="dodge", pointsize=2, href = list(0.05,0.10,0.25))
#quickplot.ssp(df = testdata, plottype="histogram", xval = Paired2, add.density=TRUE, grpcolour=DifferentGroup3, plotpos="dodge", pointsize=2, href = list(0.05,0.10,0.25))

#quickplot.ssp(df = testdata, plottype="barplot", xval = Paired2, grpcolour=DifferentGroup3, plotpos="dodge", pointsize=2, href = list(0.05,0.10,0.25), vref = list(20))
#quickplot.ssp(df = testdata, plottype="barplot", xval = Paired2, grpcolour=DifferentGroup3, plotpos="fill", pointsize=2, href = list(0.05))

#quickplot.ssp(df = testdata, plottype="histogram", xval = normalx, grpcolour=DifferentGroup3, add.density = TRUE, 
#              plotpos="dodge", pointsize=2, histdensity = TRUE, add.normal = TRUE) 
#quickplot.ssp(df = testdata, plottype="histogram", xval = normalx, grpcolour=DifferentGroup3, add.density = TRUE, 
#              plotpos="dodge", pointsize=2, histdensity = FALSE, add.normal = TRUE, vpanels = list(DifferentGroup2)) 
#quickplot.ssp(df = testdata2, plottype="histogram", xval = x, grpcolour=y, add.density = TRUE, 
#              plotpos="dodge", pointsize=2, histdensity = TRUE, add.normal = TRUE, vpanels = list(y)) 

#quickplot.ssp(df = testdata, plottype="boxplot", xval = Group2, yval = NotNormalDependent, grpcolour = BigGroup, hpanels = list(ChiSquareNotSig), vpanels=(ChiGroup2Small))
#quickplot.ssp(df = testdata, plottype="boxplot", xval = Group2, yval = NotNormalDependent, grpcolour = BigGroup, href=list(20), flipaxis = TRUE)
#quickplot.ssp(df = testdata, plottype="boxplot", xval = Group2, yval = Paired2NotDifferent, grpcolour = BigGroup, flipaxis = TRUE)

#quickplot.ssp(df = testdata, plottype="scatterplot", xval = Paired2Different, yval = NoCorrelate, grpcolour = DifferentGroup3, grpsize = NotNormalDependent, 
#              add.smooth = FALSE, href = list(55,60,65), hpanels = list(ChiSquareNotSig), transparency=0.4)
#quickplot.ssp(df = testdata, plottype="scatterplot", xval = Paired2Different, yval = NoCorrelate, grpcolour = DifferentGroup3, grpsize = NotNormalDependent, 
#              add.smooth = TRUE, href = list(55,60,65), vpanels = list(ChiSquareNotSig), transparency=0.7)
#quickplot.ssp(df = testdata, plottype="scatterplot", xval = Paired2Different, yval = NoCorrelate, grpcolour = DifferentGroup3, grpsize = NotNormalDependent, 
#              add.smooth = TRUE, href = list(60), hpanels = list(ChiSquareNotSig), transparency=0.7, freey = "free_y")
#quickplot.ssp(df = testdata, plottype="scatterplot", xval = Paired2Different, yval = NoCorrelate, grpsize = NotNormalDependent, 
#              add.smooth = TRUE, href = list(60), hpanels = list(ChiSquareNotSig), transparency=0.7, freey = "free_y", ylogsqrt = "log")
#quickplot.ssp(df = testdata, plottype="scatterplot", xval = Paired2Different, yval = normalx, grpsize = NotNormalDependent, 
#              add.smooth = TRUE, href = list(60), hpanels = list(ChiSquareNotSig), transparency=0.7, freey = "free_y", ylogsqrt = "sqrt", flipaxis = TRUE)

#quickplot.ssp(df = testdata, plottype="lineplot", xval = Paired2Different, yval = NoCorrelate)
#quickplot.ssp(df = testdata, plottype="lineplot", xval = Paired2Different, yval = NoCorrelate, ylogsqrt = "log", flipaxis = TRUE)

#quickplot.ssp(df = testdata, plottype="violinplot", xval = ChiSquareNotSig, yval = Paired2Different, grpcolour = ChiSquareNotSig, flipaxis = TRUE)

#quickplot.ssp(df = testdata, plottype="densityplot", xval = Paired2Different, linegroup = ChiSquareNotSig, hpanels = list(ChiSquareNotSig),
#              pointsize = 2, xlabel = "Paired2Different without missing values", ylabel = "Density of Paired2Different", flipaxis = TRUE)












