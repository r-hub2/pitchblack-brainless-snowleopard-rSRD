#include "SRD_Utility.h"

const char* crs(Rcpp::String element) {
  // convert Rcpp string to c
  return (element.get_cstring());
}

void convertToCMatrix(Rcpp::DataFrame matrix, std::vector<std::vector<double>>& matrix_c) {
  int nCols = matrix.size();
  int nRows = matrix.nrows();
  
  for(int i=0; i<nRows; i++) {
    Rcpp::NumericVector row(nCols);
    
    for (int j=0; j<nCols;j++) {
      Rcpp::NumericVector column = matrix[j];
      row[j] = column[i] ;
    }
    std::vector<double> tmp = Rcpp::as<std::vector<double>>(row);
    matrix_c.push_back(tmp);
  }
}

void convertToDataFrame(std::vector<std::vector<double>> matrix_c, Rcpp::DataFrame& matrix, Rcpp::CharacterVector matrix_names) {
  
  if (matrix_names.size() != 0 && (matrix_names.size() == matrix_c[0].size())) {
    for (int i = 0; i < matrix_c[0].size(); i++)
    {
      Rcpp::String name = matrix_names[i];
      Rcpp::NumericVector tmp_r;
      for (int j = 0; j < matrix_c.size(); j++)
      {
       tmp_r.push_back(matrix_c[j][i]);
      }
      matrix.push_back(tmp_r, name);
    }
  }
  else {
    for (int i = 0; i < matrix_c[0].size(); i++)
    {
      Rcpp::String name = i;
      Rcpp::NumericVector tmp_r;
      for (int j = 0; j < matrix_c.size(); j++)
      {
        tmp_r.push_back(matrix_c[j][i]);
      }
      matrix.push_back(tmp_r, name);
    }
  }
}
