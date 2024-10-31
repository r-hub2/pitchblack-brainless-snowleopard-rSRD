#' @name utilsMaxSRD
#' @title utilsMaxSRD
#' @export utilsMaxSRD
#' @description Calculates the maximum distance between two rankings of size n.
#' This function is used to normalize SRD values.
#' @aliases utilsMaxSRD
#' @author Dennis Horn
#' @param rowsCount The number of rows in the SRD calculation.
#' @return The maximum achievable SRD value.
#' @examples 
#' maxSRD <- rSRD::utilsMaxSRD(5)
utilsMaxSRD <- function(rowsCount) {
  k <- floor(rowsCount / 2)
  if ((rowsCount %% 2) == 0) {
    # even number
    return(2 * k * k)
  } else {
    # odd number
    return(2 * k * (k + 1))
  }
}

#' @title plotPermTest
#' @name plotPermTest
#' @aliases plotPermTest
#' @export plotPermTest
#' @author Linus Olsson \email{linusmeol@@gmail.com}
#' @description Plots the permutation test for the given data frame by
#' using the simulation data created by the calculateSRDDistribution() function. 
#' @param df A DataFrame.
#' @param simulationData The output of the calculateSRDDistribution() function.
#' @param densityToDistr Flag to display the cumulative distribution function instead of the probability density.
#' @return None.
#' @examples 
#' \donttest{
#' df <- data.frame(
#' A=c(32, 52, 44, 44, 47),
#' B=c(73, 75, 65, 76, 70),
#' C=c(60, 59, 57, 55, 60),
#' D=c(35, 24, 44, 83, 47),
#' E=c(41, 52, 46, 50, 65))
#' 
#' simulationData <- rSRD::calculateSRDDistribution(df)
#' plotPermTest(df, simulationData)
#' }
plotPermTest  <- function(df, simulationData, densityToDistr = FALSE) {
  
  
  srd_dist <- simulationData$SRD_Distribution
  srd_values <- calculateSRDValues(df)
  
  if(densityToDistr) {
    freq <- cumsum(srd_dist$relative_frequency)
  }
  else {
    freq <- srd_dist$relative_frequency
  }
  
  gDf <- data.frame(
    x = srd_dist$SRD_value,
    y = freq
  )
  
  vertLines <- c(simulationData$xx1, simulationData$median, simulationData$xx19)
  
  prob_levels <- data.frame(Ref  = c("XX1","Med", "XX19"),
                            vals = vertLines,
                            stringsAsFactors = FALSE)
  
  srd_solution <- data.frame(Ref = names(df)[-length(df)],
                             vals = srd_values,
                             stringsAsFactors = FALSE)
  
  graph <- ggplot2::ggplot(gDf, ggplot2::aes_string(x="SRD_Value", y="Relative_frequency")) +
    #ggplot2::geom_vline(xintercept = vertLines, linetype = 'solid', colour="red") +
    geom_vline(mapping = aes(xintercept = .data$vals),
               colour = "gray20",
               linetype = 'dashed',
               data = prob_levels,
               show.legend = FALSE) +
    ggrepel::geom_text_repel(mapping = aes(x = .data$vals,
                            y = 0,
                            label = .data$Ref,
                            hjust = -0.5,
                            vjust = -0.5),
              colour ="gray20",
              show.legend = FALSE,
              data = prob_levels) +
    geom_vline(mapping = aes(xintercept = .data$vals,
                             colour = .data$Ref),
               data = srd_solution,
               show.legend = TRUE) +
    ggrepel::geom_text_repel(mapping = aes(x = .data$vals,
                            y = sample(freq, length(df)-1, replace = TRUE) ,
                            label = .data$Ref,
                            colour = .data$Ref,
                            hjust = -0.5,
                            vjust = 1),
              show.legend = TRUE,
              data = srd_solution) +
    #ggplot2::geom_vline(xintercept = srd_values, linetype = 'dashed', colour="blue") +
    ggplot2::geom_point(mapping = ggplot2::aes_string("x", "y"))
  
  plot(graph)
  
}