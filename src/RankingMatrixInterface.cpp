#include "SRDInterface.h"
#include "SRD_Utility.h"
#include "Ranking_Matrix.h"

#include <vector>

//' @title utilsRankingMatrix
//' @name utilsRankingMatrix
//' @aliases utilsRankingMatrix
//' @author Balázs R. Sziklai \email{sziklai.balazs@@krtk.hu}, Linus Olsson \email{linusmeol@@gmail.com}
//' @description R interface to perform the rank transformation on the columns of the input data frame. 
//' Ties are resolved by fractional ranking.
//' @param data_matrix A DataFrame.
//' @return A DataFrame containing the ranking matrix.
//' @examples
//' df <- data.frame(
//' A=c(32, 52, 44, 44, 47),
//' B=c(73, 75, 65, 76, 70),
//' C=c(60, 59, 57, 55, 60),
//' D=c(35, 24, 44, 83, 47),
//' E=c(41, 52, 46, 50, 65))
//' 
//' utilsRankingMatrix(df)
// [[Rcpp::export]]
Rcpp::DataFrame utilsRankingMatrix(Rcpp::DataFrame data_matrix) {

  std::vector<std::vector<double>> data_matrix_c;
  data_matrix_c.clear();
  
  Rcpp::CharacterVector data_matrix_names;
  
  data_matrix_names = data_matrix.names();
  
  convertToCMatrix(data_matrix, data_matrix_c);
  
  std::vector<std::vector<double>>ranking_matrix_c = Ranking_Matrix::transform_input_matrix(data_matrix_c, -1);
  
  Rcpp::DataFrame ranking_matrix;
  
  convertToDataFrame(ranking_matrix_c, ranking_matrix, data_matrix_names);
  
  return ranking_matrix;
}
