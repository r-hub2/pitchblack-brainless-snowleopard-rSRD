utilsRangescaleColumn <- function(x){
  return((x-min(x)) / (max(x)- min(x))
  )}

utilsNormalizeColumn <- function(x){return(x / sqrt(sum(x^2)))}

utilsScaleToMax <- function(x){x / max(x)}

#' @title utilsCreateReference
#' @name utilsCreateReference
#' @aliases utilsCreateReference
#' @export utilsCreateReference
#' @author Ali Tugay Sen, Linus Olsson \email{linusmeol@@gmail.com}
#' @description Adds a new reference column based on the input DataFrame df and the given method. 
#' This function iterates over the rows and applies the given method to define the value of the reference. 
#' Available options are: max, min, median, mean and mixed. This column is appended to the DataFrame.
#' When "mixed" is specified the function will consider the refVector for creating the reference column. 
#' @param df A DataFrame.
#' @param method A string value specifying the reference creating method. Available options: max, min, median, mean and mixed.
#' @param refVector A vector of strings that specifies a method for each row.
#' Vector size should be equal to the number of rows in the DataFrame df.
#' @return Returns a new DataFrame appended with the reference column created by the method.
#' @examples
#' SRDInput <- data.frame(
#' A=c(32, 52, 44, 44, 47),
#' B=c(73, 75, 65, 76, 70),
#' C=c(60, 59, 57, 55, 60),
#' D=c(35, 24, 44, 83, 47),
#' E=c(41, 52, 46, 50, 65))
#' proc_data <- rSRD::utilsPreprocessDF(SRDInput)
#' ref <- c("min","max","min","max","mean")
#' rSRD::utilsCreateReference(proc_data, method = "mixed", ref)
utilsCreateReference <- function (df, method = "max", refVector = c()) {
  
  availableMethods = c("max", "min", "median", "mean")
  
  if(method == "mixed" && (nrow(df) == length(refVector))) {
    refCol = c()
    for (i in 1:nrow(df)){
      
      if(!is.element(refVector[i], availableMethods)) {
        stop("Method not defined. Available options are: max, min, median, mean. If you want to use 'mixed' methods,
           the refCol vector has to contain a valid method for each row of the given DataFrame.")
      }
      
      if(refVector[i] == "max") {
        refCol <- c(refCol, max(as.numeric(df[i,])))
      }
      else if(refVector[i] == "min") {
        refCol <- c(refCol, min(as.numeric(df[i,])))
      }
      else if(refVector[i] == "median") {
        refCol <- c(refCol, median(as.numeric(df[i,])))
      }
      else {
        refCol <- c(refCol, mean(as.numeric(df[i,])))
      }
     
    }
    df <- cbind(df, refCol)
    
  }
  else {
    
    if(!is.element(method, availableMethods)) {
      stop("Method not defined. Available options are: max, min, median, mean. If you want to use 'mixed' methods,
           the refCol vector has to contain a valid method for each row of the given DataFrame.")
    }
  
    if(method== "max") {
      df <- df %>% tibble::add_column("refCol":=apply(df,1,max))
    }
    else if(method == "min") {
      df <- df %>% tibble::add_column("refCol":=apply(df,1,min))
    }
    else if(method == "median") {
      df <- df %>% tibble::add_column("refCol":=apply(df,1,median))
    }
    else {
      df <- df %>% tibble::add_column("refCol":=apply(df,1,mean))
    }
  }
  return(df)
}

#' @title utilsPreprocessDF
#' @name utilsPreprocessDF
#' @aliases utilsPreprocessDF
#' @export utilsPreprocessDF
#' @author Ali Tugay Sen, Dennis Horn \email{dennishorn@@hotmail.de}, Linus Olsson \email{linusmeol@@gmail.com}
#' @description This function preprocesses the DataFrame depending on the \code{method}. 
#' @param df A DataFrame.
#' @param method A string that should contain "scale_to_unit", "standardize", "range_scale" or "scale_to_max".
#' @return Returns a new \code{df} that has a Distance Column based on the \code{nameCol}.
#' @examples
#' SRDInput <- data.frame(
#' A=c(32, 52, 44, 44, 47),
#' B=c(73, 75, 65, 76, 70),
#' C=c(60, 59, 57, 55, 60),
#' D=c(35, 24, 44, 83, 47),
#' E=c(41, 52, 46, 50, 65))
#' method <- "standardize"
#' utilsPreprocessDF(SRDInput,method)
utilsPreprocessDF <- function(df, method="range_scale"){
  if(method == "scale_to_unit"){
    df <- mapply(utilsNormalizeColumn, df)  
  }
  else if(method =="standardize"){
    df <- data.frame(scale(df))
  }
  else if(method == "range_scale"){
    df <- mapply(utilsRangescaleColumn, df)
  }
  else if (method == "scale_to_max") {
    df <- mapply(utilsScaleToMax, df)
  }
  else {
    stop("This method is not available, so please use one of the following: scale_to_unit, standardize, range_scale or scale_to_max.")
  }
  return (data.frame(df))
}
