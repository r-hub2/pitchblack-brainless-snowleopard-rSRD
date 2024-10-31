#' @name plotHeatmapSRD
#' @title plotHeatmapSRD
#' @aliases plotHeatmapSRD
#' @author Attila Gere \email{gereattilaphd@@gmail.com}, Linus Olsson \email{linusmeol@@gmail.com}, Jochen Staudacher \email{jochen.staudacher@@hs-kempten.de}
#' @description Heatmap is generated based on the pairwise distance - measured in SRD - of the columns. 
#' Each column is set as reference once, then SRD values are calculated for the other columns. 
#' @param df A DataFrame.
#' @param output_to_file Logical. If true, the distance matrix will be saved to the hard drive. 
#' @param color Vector of colors used for the image. Defaults to colors \code{utilsColorPalette}.
#' @return Returns a heatmap and the corresponding distance matrix.
#' @export plotHeatmapSRD
#' @examples 
#' srdInput <- data.frame(
#' A=c(32, 52, 44, 44, 47),
#' B=c(73, 75, 65, 76, 70),
#' C=c(60, 59, 57, 55, 60),
#' D=c(35, 24, 44, 83, 47),
#' E=c(41, 52, 46, 50, 65))
#' 
#' plotHeatmapSRD(srdInput)
#' 
#' mycolors<- c("#e3f2fd", "#bbdefb", "#90caf9","#64b5f6","#42a5f5",
#'              "#2196f3","#1e88e5","#1976d2","#1565c0","#0d47a1")
#' plotHeatmapSRD(srdInput, color=mycolors)
#'
plotHeatmapSRD <- function (df, output_to_file=FALSE, color=utilsColorPalette){
cnames<-colnames(df)
SRDs<-matrix(ncol=length(cnames)-1, nrow=length(cnames))
for(i in 1:length(cnames)) { 
  SRDdf2<-utilsDetailedSRDNoChars(df,referenceCol=cnames[i])
  SRDs[i,]<-as.matrix(SRDdf2[nrow(SRDdf2), !apply(is.na(SRDdf2), 2, any)])
  }
rownames(SRDs)<-colnames(df)
SRDs<-t(SRDs)
m<-matrix(ncol=ncol(SRDs),nrow=ncol(SRDs))
diag(m)<-0
m[upper.tri(m)]<-SRDs[upper.tri(SRDs)]
m[lower.tri(m)]<-SRDs[lower.tri(SRDs, TRUE)]
m<-m/utilsMaxSRD(nrow(df))
colnames(m)<-colnames(df)
rownames(m)<-colnames(df)
hmap<-heatmap.2(m,
                distfun = function (x) {as.dist(x)},
                symm = TRUE,
                dendrogram = "both",
                srtCol = 45,
                key = FALSE,
                trace = "none",
                margins = c(10, 10),
                symkey = TRUE,
                col = color)
invisible(hmap)
if (output_to_file==TRUE){
  dist_reordered<-m[rev(hmap$rowInd),hmap$colInd]
  write.table(m,"SRDdistanceMatrix.csv",sep=";", dec='.',col.names=NA)
  } 
}


