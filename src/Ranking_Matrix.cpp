#include "Ranking_Matrix.h"
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

std::vector<std::vector<double>> Ranking_Matrix::transform_input_matrix(std::vector<std::vector<double>> data_matrix, int precision)
{
	if (precision >= 0)
	{
		for (int i = 0; i < data_matrix.size(); ++i)
		{
			for (int j = 0; j < data_matrix[0].size(); ++j)
			{
				data_matrix[i][j] = floor(pow(10, precision) * data_matrix[i][j]);
			}
		}
	}

	

	std::vector<std::vector<double>> ranking_matrix(data_matrix.size(), std::vector<double>(data_matrix[0].size(),0));

	for (int i = 0; i < data_matrix[0].size(); ++i)
	{
		double max_value = data_matrix[0][i]; 
		double min_value = data_matrix[0][i];

		for (int j = 0; j < data_matrix.size(); ++j)
		{
			if (data_matrix[j][i] < min_value)min_value = data_matrix[j][i];
			if (data_matrix[j][i] > max_value)max_value = data_matrix[j][i];
		}


		int ranks_filled = data_matrix.size();
		while (ranks_filled > 0)
		{

			std::vector<int> indices;
			double tmp_max = min_value;
			double number_of_tied_values = 0;

			for (int j = 0; j < data_matrix.size(); ++j)
			{
				if (data_matrix[j][i] >= tmp_max)
				{
					tmp_max = data_matrix[j][i];
				}
			}

			for (int j = 0; j < data_matrix.size(); ++j)
			{
				if (data_matrix[j][i] == tmp_max)
				{
					indices.push_back(j);
					++number_of_tied_values;
					data_matrix[j][i] = min_value - 1;
				}
			}

			for (int j = 0; j < indices.size(); ++j)
			{
				ranking_matrix[indices[j]][i] = (ranks_filled + ranks_filled - (number_of_tied_values - 1)) / 2;
			}

			ranks_filled -= number_of_tied_values;
		}

	}


	return ranking_matrix;
}
