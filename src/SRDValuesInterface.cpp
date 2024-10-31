#include "SRDInterface.h"
#include "SRD_Utility.h"
#include "Ranking_Matrix.h"
#include "Manhattan_Dist.h"

#include <vector>
#include <fstream>

//' @title calculateSRDValues
//' @name calculateSRDValues
//' @aliases calculateSRDValues
//' @author Balázs R. Sziklai \email{sziklai.balazs@@krtk.hu}, Linus Olsson \email{linusmeol@@gmail.com}
//' @description R interface to calculate SRD values.
//' To test the results' significance run calculateSRDDistribution(). 
//' For more information about SRD scores and their validation 
//' see Héberger and Kollár-Hunek (2011).
//' @param data_matrix A DataFrame.
//' @param output_to_file Boolean flag to enable file output.
//' @return A vector containing the SRD values.
//' @references Héberger K., Kollár-Hunek K. (2011) 
//' "Sum of ranking differences for method discrimination and its validation:
//' comparison of ranks with random numbers", Journal of Chemometrics, 25(4), pp. 151–158.
//' @examples
//' df <- data.frame(
//' A=c(32, 52, 44, 44, 47),
//' B=c(73, 75, 65, 76, 70),
//' C=c(60, 59, 57, 55, 60),
//' D=c(35, 24, 44, 83, 47),
//' E=c(41, 52, 46, 50, 65))
//' 
//' calculateSRDValues(df)
// [[Rcpp::export]]
std::vector<double> calculateSRDValues(Rcpp::DataFrame data_matrix, bool output_to_file = false) {
  
  std::vector<std::vector<double>> data_matrix_c;
  data_matrix_c.clear();
  
  Rcpp::CharacterVector data_matrix_names = data_matrix.names();
    
  convertToCMatrix(data_matrix, data_matrix_c);

  std::vector<std::vector<double>>ranking_matrix_c = Ranking_Matrix::transform_input_matrix(data_matrix_c, -1);
  
  Rcpp::DataFrame ranking_matrix;
  
  convertToDataFrame(ranking_matrix_c, ranking_matrix, data_matrix_names);
  
  std::vector<double> reference_column;
  std::vector<double> SRD_values;
      for (int j = 0; j < ranking_matrix_c.size(); ++j)
      {
          reference_column.push_back(ranking_matrix_c[j][ranking_matrix_c[0].size() - 1]);
      }
  
      for (int j = 0; j < ranking_matrix_c[0].size() - 1; ++j)
      {
          std::vector<double> actual_column;
          for (int k = 0; k < ranking_matrix_c.size(); ++k)
          {
              actual_column.push_back(ranking_matrix_c[k][j]);
          }
          
          SRD_values.push_back(Manhattan_Dist::manhattan(actual_column, reference_column));
  
      }
      
      if (output_to_file) {
      
        std::ofstream file_out;
        std::string file_name = "results.txt";
        file_out.open(file_name);
        if (!file_out.is_open())
        {
          Rcpp::stop("Could not open output file.");
        }
        
        file_out << "SRD values\n";
        for (int j = 0; j < ranking_matrix_c[0].size() - 1; ++j)
        {
          file_out << "Col " << j + 1 << "; ";
        }
        file_out << std::endl;
        
        for (int j = 0; j < ranking_matrix_c[0].size() - 1; ++j)
        {
          
          file_out << SRD_values[j] << "; ";
        }
    
        file_out << "\nRanking Matrix\n";
        for (int i = 0; i < ranking_matrix_c.size(); ++i)
        {
            for (int j = 0; j < ranking_matrix_c[0].size() - 1; ++j)
            {
                file_out << ranking_matrix_c[i][j] <<";";
            }
            file_out <<  reference_column[i] << "\n";
        }
    
        file_out.close();	  
      }
  
  return SRD_values;
}
