#include "SRDInterface.h"
#include "SRD_Utility.h"
#include "Ranking_Matrix.h"
#include "Ranking_generator.h"
#include "Manhattan_Dist.h"
#include "Dynamic_Stats.h"
#include "Distribution_Stats.h"

#include <vector>
#include <fstream>

//' @title calculateSRDDistribution
//' @name calculateSRDDistribution
//' @aliases calculateSRDDistribution
//' @author Balázs R. Sziklai \email{sziklai.balazs@@krtk.hu}, Linus Olsson \email{linusmeol@@gmail.com}
//' @description R interface to calculate the SRD distribution that corresponds to the data. 
//' @param data_matrix A DataFrame.
//' @param option A char to specify how ties are generated in the simulation.
//'        The following options are available:
//'        \itemize{
//'        \item 'n': There are no ties for the solution vectors, the reference vector is fixed.
//'        \item 'r': There are no ties. Both the column vector and the reference are generated randomly.
//'        \item 't': Ties occur with a fixed probability specified by the user for both the solution vectors and the reference vector.
//'        \item 'p': Ties occur with a fixed probability specified by the user for the solution vectors, the reference vector is fixed.
//'        \item 'd': Tie distribution reflects the tie frequencies displayed by the solution vectors, the reference vector is fixed.
//'        \item As default (recommended): Tie distribution reflects the tie frequencies displayed in the reference, the reference vector is fixed.
//'        }
//' @param tie_probability The probability with which ties can occur.
//' @param output_to_file Boolean flag to enable file output.
//' @return A List containing the SRD distribution and related descriptive statistics. xx1 value 
//' indicates the 5 percent significance threshold. SRD values falling between xx1 and xx19 are not 
//' distinguishable from SRD scores of random rankings, while an SRD score higher than xx19 indicates 
//' that the solution ranks the objects in a reverse order (with 5 percent significance). 
//' @examples
//' \donttest{
//' df <- data.frame(
//' A=c(32, 52, 44, 44, 47),
//' B=c(73, 75, 65, 76, 70),
//' C=c(60, 59, 57, 55, 60),
//' D=c(35, 24, 44, 83, 47),
//' E=c(41, 52, 46, 50, 65))
//' 
//' calculateSRDDistribution(df, option = 'p', tie_probability = 0.5)
//' }
// [[Rcpp::export]]
Rcpp::List calculateSRDDistribution(Rcpp::DataFrame data_matrix, char option = 'f', double tie_probability = 0, bool output_to_file = false) {

    std::vector<std::vector<double>> data_matrix_c;
    data_matrix_c.clear();

    convertToCMatrix(data_matrix, data_matrix_c);

    if ((option == 't' || option == 'p') && (tie_probability > 1 || tie_probability < 0)) {
        Rcpp::stop("Options 't' and 'p' require a tie probability between 0 and 1.");
    }

    std::vector<std::vector<double>>ranking_matrix_c = Ranking_Matrix::transform_input_matrix(data_matrix_c, -1);

    int size_of_ranking = data_matrix_c.size();
    std::vector<double>fixed_reference;

    for (int i = 0; i < ranking_matrix_c.size(); ++i)
    {
        fixed_reference.push_back(ranking_matrix_c[i][ranking_matrix_c[0].size() - 1]);
    }

    double tie_frequency = 0;
    std::vector<double> tie_distribution;
    int round = 1;
    Distribution_Stats Simulation_results;
    Dynamic_Stats Simulation_results_dyn; //The mean and std. deviation are calculated in a dynamic fashion.
    double sim_size = 1000000;

    switch (option)
    {
    case 'n': //There are no ties for the solution vectors. Reference vector is fixed.
        while (round <= sim_size)
        {
            std::vector<double>solution_ranking = Ranking_generator::randomRanking_noties(size_of_ranking);
            double distance = Manhattan_Dist::manhattan(fixed_reference, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
        break;
    case 'r': //There are no ties. Both the column vector and the reference is generated randomly.
        while (round <= sim_size)
        {
            std::vector<double>reference_ranking = Ranking_generator::randomRanking_noties(size_of_ranking);
            std::vector<double>solution_ranking = Ranking_generator::randomRanking_noties(size_of_ranking);
            double distance = Manhattan_Dist::manhattan(reference_ranking, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
        break;
    case 't'://Ties occur with a fixed probability specified by the user for both the solution vectors and the reference vector.
        while (round <= sim_size)
        {
            std::vector<double>reference_ranking = Ranking_generator::randomRanking(size_of_ranking, tie_probability);
            std::vector<double>solution_ranking = Ranking_generator::randomRanking(size_of_ranking, tie_probability);
            double distance = Manhattan_Dist::manhattan(reference_ranking, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
        break;
    case 'p'://Ties occur with a fixed probability specified by the user for the solution vectors. Reference vector is fixed.
        while (round <= sim_size)
        {
            std::vector<double>solution_ranking = Ranking_generator::randomRanking(size_of_ranking, tie_probability);
            double distance = Manhattan_Dist::manhattan(fixed_reference, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
        break;
    case 'd': //Tie distribution reflects the tie frequencies displayed by the solution vectors. Reference is fixed.
        for (int i = 0; i < data_matrix_c[0].size() - 1; ++i) //The last column is discarded (we are not interested in the number of ties in the reference)
        {
            std::vector<std::vector<double>>sorted_matrix(data_matrix_c); //We need to sort the input to detect tie frequency.
            tie_frequency = 0;
            bool needs_sort = true;
            
            while (needs_sort)
            {
              needs_sort = false;
          
              for (int j = 0; j < sorted_matrix.size()-1; ++j)
              {
                  if (sorted_matrix[j][i] > sorted_matrix[j + 1][i]) {
                      std::swap(sorted_matrix[j][i], sorted_matrix[j + 1][i]);
                      needs_sort = true;
                  }
              }
            }
        
            for (int j = 0; j < sorted_matrix.size() - 1; ++j)
            {
                if (sorted_matrix[j][i] == sorted_matrix[j + 1][i])++tie_frequency;
            }
            
            tie_distribution.push_back(tie_frequency / (sorted_matrix.size() - 1));
        }
      
        while (round <= sim_size)
        {
            int i = int(round) % tie_distribution.size();
            tie_probability = tie_distribution[i];
            std::vector<double>solution_ranking = Ranking_generator::randomRanking(size_of_ranking, tie_probability);
        
            double distance = Manhattan_Dist::manhattan(fixed_reference, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
        break;
    default: //Tie distribution reflects the tie frequencies displayed in the reference. Reference is fixed.
         std::vector<double> sorted_reference;
        for (int j = 0; j < data_matrix_c.size(); ++j)
        {
            sorted_reference.push_back(data_matrix_c[j][data_matrix_c[0].size() - 1]);
        }
        bool needs_sort = true;
        while (needs_sort)
        {
            needs_sort = false;
            for (int j = 0; j < sorted_reference.size() - 1; ++j)
            {
                if (sorted_reference[j] > sorted_reference[j + 1])
                {
                    std::swap(sorted_reference[j], sorted_reference[j + 1]);
                    needs_sort = true;
                }
            }
        }

        for (int j = 0; j < data_matrix_c.size() - 1; ++j)
        {
            if (sorted_reference[j] == sorted_reference[j + 1])++tie_frequency;
        }

        tie_probability = tie_frequency / (data_matrix_c.size() - 1);
        while (round <= sim_size)
        {
            std::vector<double>solution_ranking = Ranking_generator::randomRanking(size_of_ranking, tie_probability);

            double distance = Manhattan_Dist::manhattan(fixed_reference, solution_ranking);
            Simulation_results.push_data(distance);
            Simulation_results_dyn.push_data(distance);
            round++;
        }
    }

    std::vector<std::vector<double>>distribution_c = Simulation_results.get_distribution();

    if(output_to_file) {

        std::ofstream file_out;
        std::string file_name = "results.txt";
        file_out.open(file_name);
        if (!file_out.is_open())
        {
          Rcpp::stop("Could not open output file.");
        }

        file_out << "SRD distribution\n";

        for (auto i = 0; i < distribution_c[0].size(); i++)
        file_out << distribution_c[0][i] << ";"<< distribution_c[1][i] / sim_size << std::endl;
        file_out << "xx1: " << Simulation_results.get_xx1(sim_size) << std::endl; //5% threshold
        file_out << "q1: " << Simulation_results.get_q1(sim_size) << std::endl;
        file_out << "med: " << Simulation_results.get_median(sim_size) << std::endl;
        file_out << "q3: " << Simulation_results.get_q3(sim_size) << std::endl;
        file_out << "xx19: " <<  Simulation_results.get_xx19(sim_size) << std::endl; //95% threshold
        file_out << "avg: " << Simulation_results_dyn.get_avg() << std::endl;
        file_out << "std_dev: " << Simulation_results_dyn.get_std_dev() << std::endl;
        file_out.close();

    }

    Rcpp::DataFrame distribution;
    Rcpp::CharacterVector distribution_names = {"SRD_value", "relative_frequency"};

    for (auto i = 0; i < distribution_c[0].size(); i++) {
      //distribution_c[0][i] = 100 * distribution_c[0][i];
      distribution_c[1][i] = distribution_c[1][i] / sim_size;
    }

    for (int i = 0; i < distribution_c.size(); i++){
      std::vector<double> tmp = distribution_c[i];
      Rcpp::NumericVector tmp_r = Rcpp::wrap(tmp);
      Rcpp::String name = distribution_names[i];
      distribution.push_back(tmp_r, name);
    }


    Rcpp::List results;
    results["SRD_Distribution"] = distribution;
    results["xx1"] = Simulation_results.get_xx1(sim_size);
    results["q1"] = Simulation_results.get_q1(sim_size);
    results["median"] = Simulation_results.get_median(sim_size);
    results["q3"] = Simulation_results.get_q3(sim_size);
    results["xx19"] = Simulation_results.get_xx19(sim_size);
    results["avg"] = Simulation_results_dyn.get_avg();
    results["std_dev"] = Simulation_results_dyn.get_std_dev();

    return results;
}
