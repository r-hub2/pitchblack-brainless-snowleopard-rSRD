#ifndef RANKING_MATRIX_H
#define RANKING_MATRIX_H


#include <iostream>
#include <vector>
#include <string>


class Ranking_Matrix {

public:


	static std::vector<std::vector<double>> transform_input_matrix(std::vector<std::vector<double>> data_matrix, int precision);
};

#endif /* RANKING_MATRIX_H */
