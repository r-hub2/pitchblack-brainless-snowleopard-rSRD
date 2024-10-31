#ifndef SRD_UTILITY_H
#define SRD_UTILITY_H

#include "SRDInterface.h"


const char* crs(Rcpp::String element);

void convertToCMatrix(Rcpp::DataFrame matrix, std::vector<std::vector<double>>& matrix_c);

void convertToDataFrame(std::vector<std::vector<double>> matrix_c, Rcpp::DataFrame& matrix, Rcpp::CharacterVector matrix_names);

#endif /* SRD_UTILITY_H */
