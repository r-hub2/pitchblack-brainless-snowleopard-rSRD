#include "SRDInterface.h"
#include "SRD_Utility.h"
#include "Cross_Validation.h"

#include <vector>
#include <fstream>


bool validateCVMethod(std::string cvmethod);

// [[Rcpp::export]]
Rcpp::List calculateCrossValidationAdapter(Rcpp::DataFrame data_matrix, Rcpp::String method, int number_of_folds, int precision, bool output_to_file) {

    std::vector<std::vector<double>> data_matrix_c;
    data_matrix_c.clear();
      
    convertToCMatrix(data_matrix, data_matrix_c);

    std::string s_method = crs(method);    
    
    if (!validateCVMethod(s_method)){
        Rcpp::Rcout << "The following Cross Validation methods are available:" << std::endl;
        Rcpp::Rcout << "Wilcoxon" << std::endl;
        Rcpp::Rcout << "Alpaydin" << std::endl;
        Rcpp::Rcout << "Dietterich" << std::endl;
        Rcpp::stop("Invalid Cross Validation method.");
    }
    
    std::vector<int> result;
    result.clear();
    
    std::vector<double> test_statistics;
    test_statistics.clear();
    
    std::vector<std::vector<double>> cv_srd;
    cv_srd.clear();
    
    std::vector<std::vector<double>> boxplot_values;
    boxplot_values.clear();
    
    if (s_method == "Wilcoxon") {
        result = Cross_Validation::Wilcoxon(test_statistics, cv_srd, boxplot_values, data_matrix_c, number_of_folds, output_to_file, precision);
    }
    else if (s_method == "Alpaydin") {
        result = Cross_Validation::Alpaydin(test_statistics, cv_srd, boxplot_values, data_matrix_c, number_of_folds, output_to_file, precision);
    }
    else { //Dietterich
        result = Cross_Validation::Dietterich(test_statistics, cv_srd, boxplot_values, data_matrix_c, number_of_folds, output_to_file, precision);
    }
    
    Rcpp::NumericVector r_result = Rcpp::wrap(result);
    Rcpp::NumericVector r_test_statistics = Rcpp::wrap(test_statistics);
    
    Rcpp::DataFrame r_cv_srd;
    Rcpp::DataFrame r_boxplot_values;
    Rcpp::CharacterVector data_matrix_names;
    
    
    
    convertToDataFrame(cv_srd, r_cv_srd, data_matrix_names);
    convertToDataFrame(boxplot_values, r_boxplot_values, data_matrix_names);
    
    Rcpp::List results;
    
    results["new_column_order_based_on_folds"] = r_result;
    results["test_statistics"] = r_test_statistics;
    results["SRD_values_of_different_folds"] = r_cv_srd;
    results["boxplot_values"] = r_boxplot_values;
    
    
    return results;
}


bool validateCVMethod(std::string s_cvmethod) {
    
    if (s_cvmethod == "Wilcoxon") {
        return true;
    }
    else if (s_cvmethod == "Alpaydin") {
        return true;
    }
    else if (s_cvmethod == "Dietterich") {
        return true;
    }
    else {
        return false;
    }  
}
