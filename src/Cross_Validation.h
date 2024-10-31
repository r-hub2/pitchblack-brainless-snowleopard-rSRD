#ifndef CROSS_VALIDATION_H
#define CROSS_VALIDATION_H

#include <iostream>
#include <vector>
#include <random>




class Cross_Validation
{

public:

	static std::vector<int> Wilcoxon(std::vector<double>& test_statistics, std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold_spec, bool print_it, int precision);
	static std::vector<int> Alpaydin(std::vector<double>& test_statistics, std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold_spec, bool print_it, int precision);
	static std::vector<int> Dietterich(std::vector<double>& test_statistics,std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold_spec, bool print_it, int precision);

};

#endif /* CROSS_VALIDATION_H */
