#include "Cross_Validation.h"
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "Ranking_Matrix.h"
#include "Manhattan_Dist.h"
#include "Distribution_Stats.h"
#include <Rcpp.h>
#include <fstream>
#include <sstream>


std::vector<int> Cross_Validation::Wilcoxon(std::vector<double>& test_statistics, std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold, bool print_it, int precision)
{

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()



    //Reinstate after simulation
    k_fold = std::min(int(data_matrix.size()), k_fold);
    if (int(data_matrix.size()) < k_fold)
    {
        Rcpp::Rcout << "Applying " << data_matrix.size() << "-fold Wilcoxon cross validation." << std::endl;
    }

    int leave_m_out = ceil(double(data_matrix.size()) / double(k_fold));

    

    std::vector <std::vector<int>> random_row_indices(k_fold, std::vector<int>(leave_m_out, 0));
    

    for (int i = 0; i < k_fold; ++i) //Generating random row indices that will be removed during the cross validation
    {

        std::vector<int> ranking;
        for (int j = 0; j < data_matrix.size(); j++)
        {
            ranking.push_back(j);
        }

        std::shuffle(ranking.begin(), ranking.end(), gen); //Random shuffling of the vector elements


        for (int j = 0; j < leave_m_out; ++j)
        {
            random_row_indices[i][j] = ranking[j];
        }

        std::sort(random_row_indices[i].begin(), random_row_indices[i].end(), std::greater<int>()); //Since we will remove rows one-by-one we need to sort the elements

        
        for (int j = 0; j < i; ++j) //If the same configuration of row indices has been used before, we generate a new one
        {
            int abs_differences = 0;

            for (int k = 0; k < leave_m_out; ++k)
            {
                abs_differences += abs(random_row_indices[j][k] - random_row_indices[i][k]);
            }

            if (abs_differences == 0) 
            {
                --i;
            }
        }
    }

    //std::vector<std::vector<double>>cv_srd; //Resulting SRD will be stored in this matrix


    for (int i = 0; i < k_fold; ++i)
    {
        std::vector<std::vector<double>>trimmed_data_matrix = data_matrix;

        for (int j = 0; j < leave_m_out; ++j)
        {
            trimmed_data_matrix.erase(trimmed_data_matrix.begin() + random_row_indices[i][j]);
        }


        std::vector<std::vector<double>>trimmed_ranking_matrix = Ranking_Matrix::transform_input_matrix(trimmed_data_matrix, precision);

        std::vector<double> reference_column;
        for (int j = 0; j < trimmed_ranking_matrix.size(); ++j)
        {
            reference_column.push_back(trimmed_ranking_matrix[j][trimmed_ranking_matrix[0].size() - 1]);
        }

        std::vector<double> tmp_row;
        for (int j = 0; j < trimmed_ranking_matrix[0].size() - 1; ++j)
        {
            std::vector<double> actual_column;
            for (int k = 0; k < trimmed_ranking_matrix.size(); ++k)
            {
                actual_column.push_back(trimmed_ranking_matrix[k][j]);
            }


            double distance = Manhattan_Dist::manhattan(actual_column, reference_column);

            tmp_row.push_back(distance);

            //std::cout<<"\ndist: "<< distance << std::endl;
        }
        cv_srd.push_back(tmp_row);
    }



    std::vector<std::vector<double>> boxplot_stats(7,std::vector<double>(cv_srd[0].size(),0));
    

    for (int j = 0; j < cv_srd[0].size(); ++j)
    {
        Distribution_Stats Column_stat;
    
        for (int i = 0; i < cv_srd.size(); ++i)
        {
            Column_stat.push_data(cv_srd[i][j]);
        }

        boxplot_stats[0][j] = Column_stat.get_min();
        boxplot_stats[1][j] = Column_stat.get_xx1(double(cv_srd.size()));
        boxplot_stats[2][j] = Column_stat.get_q1(double(cv_srd.size()));
        boxplot_stats[3][j] = Column_stat.get_median(double(cv_srd.size()));
        boxplot_stats[4][j] = Column_stat.get_q3(double(cv_srd.size()));
        boxplot_stats[5][j] = Column_stat.get_xx19(double(cv_srd.size()));
        boxplot_stats[6][j] = Column_stat.get_max();

    }
    
   
    std::vector<int> columns_ordered;
    for (int i = 0; i < cv_srd[0].size(); ++i)
    {
        columns_ordered.push_back(i+1);
    }

    bool ordering_completed = false;
    while (!ordering_completed)  //Ordering the columns by the median, if those are the same then by the Q1 values, then the Q3 values then the min then the max.
    {
        ordering_completed = true;

        for (int j = 0; j < cv_srd[0].size()-1; ++j)
        {
            if (boxplot_stats[3][j] > boxplot_stats[3][j + 1])
            {
                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                for (int i = 0; i < cv_srd.size(); ++i)
                {
                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                }
                for (int i = 0; i < boxplot_stats.size(); ++i)
                {
                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                }
                ordering_completed = false;
            }
            else if (boxplot_stats[3][j] == boxplot_stats[3][j + 1])
            {
                if (boxplot_stats[2][j] > boxplot_stats[2][j + 1])
                {
                    std::swap(columns_ordered[j], columns_ordered[j + 1]);
                    for (int i = 0; i < cv_srd.size(); ++i)
                    {
                        std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                    }
                    for (int i = 0; i < boxplot_stats.size(); ++i)
                    {
                        std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                    }
                    ordering_completed = false;
                }
                else if (boxplot_stats[2][j] == boxplot_stats[2][j + 1])
                {
                    if (boxplot_stats[4][j] > boxplot_stats[4][j + 1])
                    {
                        std::swap(columns_ordered[j], columns_ordered[j + 1]);
                        for (int i = 0; i < cv_srd.size(); ++i)
                            {
                                std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                            }
                        for (int i = 0; i < boxplot_stats.size(); ++i)
                        {
                            std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                        }
                        ordering_completed = false;
                    }
                    else if (boxplot_stats[4][j] == boxplot_stats[4][j + 1])
                    {
                        if (boxplot_stats[0][j] > boxplot_stats[0][j + 1])
                        {
                            std::swap(columns_ordered[j], columns_ordered[j + 1]);
                            for (int i = 0; i < cv_srd.size(); ++i)
                                {
                                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                }
                            for (int i = 0; i < boxplot_stats.size(); ++i)
                            {
                                std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                            }
                            ordering_completed = false;
                        }
                        else if (boxplot_stats[0][j] == boxplot_stats[0][j + 1])
                        {
                            if (boxplot_stats[6][j] > boxplot_stats[6][j + 1])
                            {
                                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                                for (int i = 0; i < cv_srd.size(); ++i)
                                {
                                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                }
                                for (int i = 0; i < boxplot_stats.size(); ++i)
                                {
                                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                                }
                                ordering_completed = false;
                            }    
                        }    
                    }
                }
            }
        }
    }
    
    boxplot_values = boxplot_stats;

    
    //std::vector<double>test_statistics;

    for (int j = 0; j < cv_srd[0].size() - 1; ++j)
    {
        std::vector<std::vector<double>>signed_ranks_input;

        for (int i = 0; i < cv_srd.size(); ++i)
        {
            std::vector<double> diff;
            diff.push_back(abs(cv_srd[i][j] - cv_srd[i][j + 1]));
            if (diff[0] != 0)
            {
                signed_ranks_input.push_back(diff);
            }
        }

        if (signed_ranks_input.size() == 0)
        {
            test_statistics.push_back(0);
            continue;
        }


        std::vector<std::vector<double>>signed_ranks = Ranking_Matrix::transform_input_matrix(signed_ranks_input, precision);



        int ind = 0;
        double w_plus = 0; //Ranks of the positive values
        double w_minus = 0; //Ranks of the negative values
        for (int i = 0; i < cv_srd.size(); ++i)
        {
            if (cv_srd[i][j] - cv_srd[i][j + 1] > 0)
            {
                w_plus += signed_ranks[ind][0];
                ++ind;
            }
            else if (cv_srd[i][j] - cv_srd[i][j + 1] < 0)
            {
                w_minus += signed_ranks[ind][0];
                ++ind;
            }
        }

        test_statistics.push_back(abs(w_plus - w_minus));

    }

    if (print_it == true)
    {   

        std::ofstream file_out;
        std::string file_name = "CV_Wilcoxon.csv";
        file_out.open(file_name);
        if (!file_out.is_open())
        {
            Rcpp::stop("Could not open output file.");
        }        

        file_out << "Wilcoxon test" << std::endl << std::endl;

        if (int(data_matrix.size()) < k_fold)
        {
            file_out << "Applying " << data_matrix.size() << "-fold Wilcoxon cross validation." << std::endl;
        }
        else
        {
        file_out << "k fold: " << k_fold << std::endl;
        }

        file_out << "leave m out: " << leave_m_out << std::endl;


        file_out << "\nIndices of the removed rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < leave_m_out; ++j)
            {
                file_out << random_row_indices[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nColumn IDs after reordering\n";
        for (int i = 0; i < columns_ordered.size(); ++i)
        {
            file_out << columns_ordered[i] << "; ";
        }
        file_out << std::endl;

        file_out << "\nSRD values of different folds\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < data_matrix[0].size() - 1; ++j)
            {
                file_out << cv_srd[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nBoxplot\n";
        for (int i = 0; i < 7; ++i)
        {
            for (int k = 0; k < cv_srd[0].size(); ++k)
            {
                if (i == 0 && k == 0)
                    file_out << "min; ";
                if (i == 1 && k == 0)
                    file_out << "xx1; ";
                if (i == 2 && k == 0)
                    file_out << "q1; ";
                if (i == 3 && k == 0)
                    file_out << "median; ";
                if (i == 4 && k == 0)
                    file_out << "q3; ";
                if (i == 5 && k == 0)
                    file_out << "xx19; ";
                if (i == 6 && k == 0)
                    file_out << "max; ";

                file_out << boxplot_stats[i][k] << "; ";
            }
            file_out << std::endl;
        }

        file_out << "\nTest statistics\n";
        for (int i = 0; i < test_statistics.size(); ++i)
        {
            file_out << test_statistics[i] << "; ";
        }
		
		file_out << "\n\nStatistical significance\n";

        std::vector<double> p_val_005 = {15,21,24,30,35,39};
        std::vector<double> p_val_010 = {15,17,22,26,29,35};

        for(int i = 0; i < test_statistics.size(); ++i)
        {
            if (k_fold < 5 || k_fold>10)
            {
                file_out <<  "n.a.; ";
            }
            else if (test_statistics[i] >= p_val_005[k_fold - 5])
            {
                file_out << "(p<0.05*); ";
            }
            else if (test_statistics[i] >= p_val_010[k_fold - 5])
            {
                file_out << "(p<0.1); ";
            }
            else
            {
                file_out << "n.s.; ";
            }

        }

    }
    return columns_ordered;
}





std::vector<int> Cross_Validation::Alpaydin(std::vector<double>& test_statistics, std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold, bool print_it, int precision)
{

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()


    if (int(data_matrix.size()*(data_matrix.size()-1)/2) < k_fold) // Triggers if the number of folds exceeds the possible number of 2-coloring the rows
    {
        k_fold = int(data_matrix.size() * (data_matrix.size() - 1) / 2);
        Rcpp::Rcout << "Applying " << int(data_matrix.size() * (data_matrix.size() - 1) / 2) << "-fold Alpaydin cross validation." << std::endl;
    }

    int leave_m_out = std::floor(data_matrix.size() / 2);

    std::vector <std::vector<int>> random_row_indices(k_fold, std::vector<int>(leave_m_out, 0));



    for (int i = 0; i < k_fold; ++i) //Generating random row indices that will be removed during the cross validation
    {

        std::vector<int> ranking;
        for (int j = 0; j < data_matrix.size(); j++)
        {
            ranking.push_back(j);
        }

        std::shuffle(ranking.begin(), ranking.end(), gen); //Random shuffling of the vector elements


        for (int j = 0; j < leave_m_out; ++j)
        {
            random_row_indices[i][j] = ranking[j];
        }


        std::sort(random_row_indices[i].begin(), random_row_indices[i].end(), std::greater<int>()); //Since we will remove rows one-by-one we need to sort the elements

        for (int j = 0; j < i; ++j) //If the same configuration of row indices has been used before, we generate a new one
        {
            int abs_differences = 0;

            for (int k = 0; k < leave_m_out; ++k)
            {
                abs_differences += abs(random_row_indices[j][k] - random_row_indices[i][k]);
            }

            if (abs_differences == 0)
            {
                //std::cout << "r" << std::endl;
                --i;
            }
        }
    }


    //Resulting SRD will be stored in three matrices (blue = selected rows, orange = left out rows, merged = blue + orange)
    std::vector<std::vector<double>>cv_srd_blue; 
    std::vector<std::vector<double>>cv_srd_orange; 
    //std::vector<std::vector<double>>cv_srd;

    for (int i = 0; i < k_fold; ++i)
    {
        std::vector<std::vector<double>>blue_data_matrix;
        std::vector<std::vector<double>>orange_data_matrix = data_matrix;

        for (int j = 0; j < leave_m_out; ++j)
        {
            blue_data_matrix.push_back(data_matrix[random_row_indices[i][j]]);
        }



        for (int j = 0; j < leave_m_out; ++j)
        {
            orange_data_matrix.erase(orange_data_matrix.begin() + random_row_indices[i][j]);
        }


        std::vector<std::vector<double>>blue_ranking_matrix = Ranking_Matrix::transform_input_matrix(blue_data_matrix, precision);
        std::vector<std::vector<double>>orange_ranking_matrix = Ranking_Matrix::transform_input_matrix(orange_data_matrix, precision);



        std::vector<double> blue_reference_column;
        for (int j = 0; j < blue_ranking_matrix.size(); ++j)
        {
            blue_reference_column.push_back(blue_ranking_matrix[j][blue_ranking_matrix[0].size() - 1]);
        }

        std::vector<double> orange_reference_column;
        for (int j = 0; j < orange_ranking_matrix.size(); ++j)
        {
            orange_reference_column.push_back(orange_ranking_matrix[j][orange_ranking_matrix[0].size() - 1]);
        }



        std::vector<double> blue_row;
        for (int j = 0; j < blue_ranking_matrix[0].size() - 1; ++j)
        {
            std::vector<double> actual_column;
            for (int k = 0; k < blue_ranking_matrix.size(); ++k)
            {
                actual_column.push_back(blue_ranking_matrix[k][j]);
            }


            double distance = Manhattan_Dist::manhattan(actual_column, blue_reference_column);

            blue_row.push_back(distance);

        }
        cv_srd_blue.push_back(blue_row);
        cv_srd.push_back(blue_row);



        std::vector<double> orange_row;
        for (int j = 0; j < orange_ranking_matrix[0].size() - 1; ++j)
        {
            std::vector<double> actual_column;
            for (int k = 0; k < orange_ranking_matrix.size(); ++k)
            {
                actual_column.push_back(orange_ranking_matrix[k][j]);
            }



            double distance = Manhattan_Dist::manhattan(actual_column, orange_reference_column);

            orange_row.push_back(distance);

        }
        cv_srd_orange.push_back(orange_row);
        cv_srd.push_back(orange_row);
    }




    std::vector<std::vector<double>> boxplot_stats(7, std::vector<double>(cv_srd[0].size(), 0));


    for (int j = 0; j < cv_srd[0].size(); ++j)
    {
        Distribution_Stats Column_stat;

        for (int i = 0; i < cv_srd.size(); ++i)
        {
            Column_stat.push_data(cv_srd[i][j]);
        }

        boxplot_stats[0][j] = Column_stat.get_min();
        boxplot_stats[1][j] = Column_stat.get_xx1(double(cv_srd.size()));
        boxplot_stats[2][j] = Column_stat.get_q1(double(cv_srd.size()));
        boxplot_stats[3][j] = Column_stat.get_median(double(cv_srd.size()));
        boxplot_stats[4][j] = Column_stat.get_q3(double(cv_srd.size()));
        boxplot_stats[5][j] = Column_stat.get_xx19(double(cv_srd.size()));
        boxplot_stats[6][j] = Column_stat.get_max();


    }
    
   

    std::vector<int> columns_ordered;
    for (int i = 0; i < cv_srd[0].size(); ++i)
    {
        columns_ordered.push_back(i + 1);
    }

    bool ordering_completed = false;
    while (!ordering_completed)  //Ordering the columns by the median, if those are the same then by the Q1 values, then the Q3 values then the min then the max.
    {
        ordering_completed = true;

        for (int j = 0; j < cv_srd[0].size() - 1; ++j)
        {
            if (boxplot_stats[3][j] > boxplot_stats[3][j + 1])
            {
                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                for (int i = 0; i < cv_srd.size(); ++i)
                {
                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                    if (i < k_fold)
                    {
                        std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                        std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                    }
                }
                for (int i = 0; i < boxplot_stats.size(); ++i)
                {
                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                }
                ordering_completed = false;
            }
            else if (boxplot_stats[3][j] == boxplot_stats[3][j + 1])
            {
                if (boxplot_stats[2][j] > boxplot_stats[2][j + 1])
                {
                    std::swap(columns_ordered[j], columns_ordered[j + 1]);
                    for (int i = 0; i < cv_srd.size(); ++i)
                    {
                        std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                        if (i < k_fold)
                        {
                            std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                            std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                        }
                    }
                    for (int i = 0; i < boxplot_stats.size(); ++i)
                    {
                        std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                    }
                    ordering_completed = false;
                }
                else if (boxplot_stats[2][j] == boxplot_stats[2][j + 1])
                {
                    if (boxplot_stats[4][j] > boxplot_stats[4][j + 1])
                    {
                        std::swap(columns_ordered[j], columns_ordered[j + 1]);
                        for (int i = 0; i < cv_srd.size(); ++i)
                        {
                            std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                            if (i < k_fold)
                            {
                                std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                            }
                        }
                        for (int i = 0; i < boxplot_stats.size(); ++i)
                        {
                            std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                        }
                        ordering_completed = false;
                    }
                    else if (boxplot_stats[4][j] == boxplot_stats[4][j + 1])
                    {
                        if (boxplot_stats[0][j] > boxplot_stats[0][j + 1])
                        {
                            std::swap(columns_ordered[j], columns_ordered[j + 1]);
                            for (int i = 0; i < cv_srd.size(); ++i)
                            {
                                std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                if (i < k_fold)
                                {
                                    std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                    std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                                }
                            }
                            for (int i = 0; i < boxplot_stats.size(); ++i)
                            {
                                std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                            }
                            ordering_completed = false;
                        }
                        else if (boxplot_stats[0][j] == boxplot_stats[0][j + 1])
                        {
                            if (boxplot_stats[6][j] > boxplot_stats[6][j + 1])
                            {
                                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                                for (int i = 0; i < cv_srd.size(); ++i)
                                {
                                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                    if (i < k_fold)
                                    {
                                        std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                        std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                                    }
                                }
                                for (int i = 0; i < boxplot_stats.size(); ++i)
                                {
                                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                                }
                                ordering_completed = false;
                            }
                        }
                    }
                }
            }
        }
    }
    
    boxplot_values = boxplot_stats;

    for (int j = 0; j < cv_srd[0].size() - 1; ++j)
    {
        std::vector<double>alpaydin_p_orange; //the p parameter in the Alpaydin formula
        std::vector<double>alpaydin_p_blue; //the p parameter in the Alpaydin formula
        std::vector<double>alpaydin_s; //the s parameter in the Alpaydin formula

        for (int i = 0; i < k_fold; ++i)
        {
            alpaydin_p_orange.push_back(cv_srd_orange[i][j] - cv_srd_orange[i][j + 1]);
            alpaydin_p_blue.push_back(cv_srd_blue[i][j] - cv_srd_blue[i][j + 1]);

            double color_avg = (alpaydin_p_orange[i] + alpaydin_p_blue[i]) / 2;
            alpaydin_s.push_back(std::pow(alpaydin_p_orange[i] - color_avg, 2) + std::pow(alpaydin_p_blue[i] - color_avg, 2));
        }

        double alpaydin_numerator = 0;
        double alpaydin_denominator = 0;

        for(int i = 0; i < k_fold; ++i)
        {
            alpaydin_numerator += std::pow(alpaydin_p_blue[i], 2);
            alpaydin_numerator += std::pow(alpaydin_p_orange[i], 2);
            alpaydin_denominator += alpaydin_s[i];
        }

		double alpaydin_formula = 0;

		if (alpaydin_denominator!=0)
		{
			alpaydin_formula = alpaydin_numerator / (2*alpaydin_denominator);
		}

        test_statistics.push_back(alpaydin_formula);
    }

    if (print_it == true)
    {

        std::ofstream file_out;
        std::string file_name = "CV_Alpaydin.csv";
        file_out.open(file_name);
        if (!file_out.is_open())
        {
            Rcpp::stop("Could not open output file.");
        }    

        file_out << "Alpaydin test" << std::endl << std::endl;

        file_out << "k fold: " << k_fold << std::endl;
        file_out << "leave m out: " << leave_m_out << std::endl;


        file_out << "\nIndices of the removed rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < leave_m_out; ++j)
            {
                file_out << random_row_indices[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nColumn IDs after reordering\n";
        for (int i = 0; i < columns_ordered.size(); ++i)
        {
            file_out << columns_ordered[i] << "; ";
        }
        file_out << std::endl;

        file_out << "\nSRD values of selected rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < data_matrix[0].size() - 1; ++j)
            {
                file_out << cv_srd_blue[i][j] << "; ";
            }

            file_out << std::endl;
        }


        file_out << "\nSRD values of left out rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < data_matrix[0].size() - 1; ++j)
            {
                file_out << cv_srd_orange[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nBoxplot\n";
        for (int i = 0; i < 7; ++i)
        {
            for (int k = 0; k < cv_srd[0].size(); ++k)
            {
                if (i == 0 && k == 0)
                    file_out << "min; ";
                if (i == 1 && k == 0)
                    file_out << "xx1; ";
                if (i == 2 && k == 0)
                    file_out << "q1; ";
                if (i == 3 && k == 0)
                    file_out << "median; ";
                if (i == 4 && k == 0)
                    file_out << "q3; ";
                if (i == 5 && k == 0)
                    file_out << "xx19; ";
                if (i == 6 && k == 0)
                    file_out << "max; ";

                file_out << boxplot_stats[i][k] << "; ";
            }
            file_out << std::endl;
        }


        file_out << "\nTest statistics\n";
        for (int i = 0; i < test_statistics.size(); ++i)
        {
            file_out << test_statistics[i] << "; ";
        }
		
		file_out << "\n\nStatistical significance\n";

        std::vector<double> p_val_005 = {4.735,4.000,3.529,3.202,2.960,2.774};
        std::vector<double> p_val_010 = {3.297,2.905,2.643,2.455,2.312,2.201};

        for (int i = 0; i < test_statistics.size(); ++i)
        {
            if (k_fold < 5 || k_fold>10)
            {
                file_out << "n.a.; ";
            }
            else if (test_statistics[i] >= p_val_005[k_fold - 5])
            {
                file_out << "(p<0.05*); ";
            }
            else if (test_statistics[i] >= p_val_010[k_fold - 5])
            {
                file_out << "(p<0.1); ";
            }
            else
            {
                file_out << "n.s.; ";
            }

        }

    }

    return columns_ordered;
}











std::vector<int> Cross_Validation::Dietterich(std::vector<double>& test_statistics, std::vector<std::vector<double>>& cv_srd, std::vector<std::vector<double>>& boxplot_values, std::vector<std::vector<double>>data_matrix, int k_fold, bool print_it, int precision)
{

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    if (int(data_matrix.size() * (data_matrix.size() - 1) / 2) < k_fold) // Triggers if the number of folds exceeds the possible number of 2-coloring the rows
    {
      k_fold = int(data_matrix.size() * (data_matrix.size() - 1) / 2);
      Rcpp::Rcout << "Applying " << int(data_matrix.size() * (data_matrix.size() - 1) / 2) << "-fold Dietterich cross validation." << std::endl;
    }
    
    int leave_m_out = std::floor(data_matrix.size() / 2);

    std::vector <std::vector<int>> random_row_indices(k_fold, std::vector<int>(leave_m_out, 0));


    for (int i = 0; i < k_fold; ++i) //Generating random row indices that will be removed during the cross validation
    {

        std::vector<int> ranking;
        for (int j = 0; j < data_matrix.size(); j++)
        {
            ranking.push_back(j);
        }

        std::shuffle(ranking.begin(), ranking.end(), gen); //Random shuffling of the vector elements


        for (int j = 0; j < leave_m_out; ++j)
        {
            random_row_indices[i][j] = ranking[j];
        }

        std::sort(random_row_indices[i].begin(), random_row_indices[i].end(), std::greater<int>()); //Since we will remove rows one-by-one we need to sort the elements

        for (int j = 0; j < i; ++j) //If the same configuration of row indices has been used before, we generate a new one
        {
            int abs_differences = 0;

            for (int k = 0; k < leave_m_out; ++k)
            {
                abs_differences += abs(random_row_indices[j][k] - random_row_indices[i][k]);
            }

            if (abs_differences == 0)
            {
                //std::cout << "r" << std::endl;
                --i;
            }
        }
    }

    //Resulting SRD will be stored in three matrices (blue = selected rows, orange = left out rows, merged = blue + orange)
    std::vector<std::vector<double>>cv_srd_blue;
    std::vector<std::vector<double>>cv_srd_orange;
    //std::vector<std::vector<double>>cv_srd;

    for (int i = 0; i < k_fold; ++i)
    {
        std::vector<std::vector<double>>blue_data_matrix;
        std::vector<std::vector<double>>orange_data_matrix = data_matrix;

        for (int j = 0; j < leave_m_out; ++j)
        {
            std::vector<double> tmp_row;
            for (int k = 0; k < data_matrix[0].size(); ++k)
            {
                tmp_row.push_back(data_matrix[random_row_indices[i][j]][k]);
            }

            blue_data_matrix.push_back(tmp_row);
        }

        for (int j = 0; j < leave_m_out; ++j)
        {
            orange_data_matrix.erase(orange_data_matrix.begin() + random_row_indices[i][j]);
        }

        std::vector<std::vector<double>>blue_ranking_matrix = Ranking_Matrix::transform_input_matrix(blue_data_matrix, precision);
        std::vector<std::vector<double>>orange_ranking_matrix = Ranking_Matrix::transform_input_matrix(orange_data_matrix, precision);

        std::vector<double> blue_reference_column;
        for (int j = 0; j < blue_ranking_matrix.size(); ++j)
        {
            blue_reference_column.push_back(blue_ranking_matrix[j][blue_ranking_matrix[0].size() - 1]);
        }

        std::vector<double> orange_reference_column;
        for (int j = 0; j < orange_ranking_matrix.size(); ++j)
        {
            orange_reference_column.push_back(orange_ranking_matrix[j][orange_ranking_matrix[0].size() - 1]);
        }



        std::vector<double> blue_row;
        for (int j = 0; j < blue_ranking_matrix[0].size() - 1; ++j)
        {
            std::vector<double> actual_column;
            for (int k = 0; k < blue_ranking_matrix.size(); ++k)
            {
                actual_column.push_back(blue_ranking_matrix[k][j]);
            }

            double distance = Manhattan_Dist::manhattan(actual_column, blue_reference_column);

            blue_row.push_back(distance);
            
        }
        cv_srd_blue.push_back(blue_row);
        cv_srd.push_back(blue_row);



        std::vector<double> orange_row;
        for (int j = 0; j < orange_ranking_matrix[0].size() - 1; ++j)
        {
            std::vector<double> actual_column;
            for (int k = 0; k < orange_ranking_matrix.size(); ++k)
            {
                actual_column.push_back(orange_ranking_matrix[k][j]);
            }

            double distance = Manhattan_Dist::manhattan(actual_column, orange_reference_column);

            orange_row.push_back(distance);

            //std::cout<<"\ndist: "<< distance << std::endl;
        }
        cv_srd_orange.push_back(orange_row);
        cv_srd.push_back(orange_row);
    }

    std::vector<std::vector<double>> boxplot_stats(7, std::vector<double>(cv_srd[0].size(), 0));


    for (int j = 0; j < cv_srd[0].size(); ++j)
    {
        Distribution_Stats Column_stat;

        for (int i = 0; i < cv_srd.size(); ++i)
        {
            Column_stat.push_data(cv_srd[i][j]);
        }

        boxplot_stats[0][j] = Column_stat.get_min();
        boxplot_stats[1][j] = Column_stat.get_xx1(double(cv_srd.size()));
        boxplot_stats[2][j] = Column_stat.get_q1(double(cv_srd.size()));
        boxplot_stats[3][j] = Column_stat.get_median(double(cv_srd.size()));
        boxplot_stats[4][j] = Column_stat.get_q3(double(cv_srd.size()));
        boxplot_stats[5][j] = Column_stat.get_xx19(double(cv_srd.size()));
        boxplot_stats[6][j] = Column_stat.get_max();

    }

    std::vector<int> columns_ordered;
    for (int i = 0; i < cv_srd[0].size(); ++i)
    {
        columns_ordered.push_back(i + 1);
    }

    bool ordering_completed = false;
    while (!ordering_completed)  //Ordering the columns by the median, if those are the same then by the Q1 values, then the Q3 values then the min then the max.
    {
        ordering_completed = true;

        for (int j = 0; j < cv_srd[0].size() - 1; ++j)
        {
            if (boxplot_stats[3][j] > boxplot_stats[3][j + 1])
            {
                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                for (int i = 0; i < cv_srd.size(); ++i)
                {
                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                    if (i < k_fold)
                    {
                        std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                        std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                    }
                }
                for (int i = 0; i < boxplot_stats.size(); ++i)
                {
                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                }
                ordering_completed = false;
            }
            else if (boxplot_stats[3][j] == boxplot_stats[3][j + 1])
            {
                if (boxplot_stats[2][j] > boxplot_stats[2][j + 1])
                {
                    std::swap(columns_ordered[j], columns_ordered[j + 1]);
                    for (int i = 0; i < cv_srd.size(); ++i)
                    {
                        std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                        if (i < k_fold)
                        {
                            std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                            std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                        }
                    }
                    for (int i = 0; i < boxplot_stats.size(); ++i)
                    {
                        std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                    }
                    ordering_completed = false;
                }
                else if (boxplot_stats[2][j] == boxplot_stats[2][j + 1])
                {
                    if (boxplot_stats[4][j] > boxplot_stats[4][j + 1])
                    {
                        std::swap(columns_ordered[j], columns_ordered[j + 1]);
                        for (int i = 0; i < cv_srd.size(); ++i)
                        {
                            std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                            if (i < k_fold)
                            {
                                std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                            }
                        }
                        for (int i = 0; i < boxplot_stats.size(); ++i)
                        {
                            std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                        }
                        ordering_completed = false;
                    }
                    else if (boxplot_stats[4][j] == boxplot_stats[4][j + 1])
                    {
                        if (boxplot_stats[0][j] > boxplot_stats[0][j + 1])
                        {
                            std::swap(columns_ordered[j], columns_ordered[j + 1]);
                            for (int i = 0; i < cv_srd.size(); ++i)
                            {
                                std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                if (i < k_fold)
                                {
                                    std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                    std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                                }
                            }
                            for (int i = 0; i < boxplot_stats.size(); ++i)
                            {
                                std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                            }
                            ordering_completed = false;
                        }
                        else if (boxplot_stats[0][j] == boxplot_stats[0][j + 1])
                        {
                            if (boxplot_stats[6][j] > boxplot_stats[6][j + 1])
                            {
                                std::swap(columns_ordered[j], columns_ordered[j + 1]);
                                for (int i = 0; i < cv_srd.size(); ++i)
                                {
                                    std::swap(cv_srd[i][j], cv_srd[i][j + 1]);
                                    if (i < k_fold)
                                    {
                                        std::swap(cv_srd_blue[i][j], cv_srd_blue[i][j + 1]);
                                        std::swap(cv_srd_orange[i][j], cv_srd_orange[i][j + 1]);
                                    }
                                }
                                for (int i = 0; i < boxplot_stats.size(); ++i)
                                {
                                    std::swap(boxplot_stats[i][j], boxplot_stats[i][j + 1]);
                                }
                                ordering_completed = false;
                            }
                        }
                    }
                }
            }
        }
    }
    
    boxplot_values = boxplot_stats;

    for (int j = 0; j < cv_srd[0].size() - 1; ++j)
    {
        std::vector<double>dietterich_p_orange; //the p parameter in the Alpaydin formula
        std::vector<double>dietterich_p_blue; //the p parameter in the Alpaydin formula
        std::vector<double>dietterich_s; //the s parameter in the Alpaydin formula

        for (int i = 0; i < k_fold; ++i)
        {
            dietterich_p_orange.push_back(cv_srd_orange[i][j] - cv_srd_orange[i][j + 1]);
            dietterich_p_blue.push_back(cv_srd_blue[i][j] - cv_srd_blue[i][j + 1]);

            double color_avg = (dietterich_p_orange[i] + dietterich_p_blue[i]) / 2;
            dietterich_s.push_back(std::pow(dietterich_p_orange[i] - color_avg, 2) + std::pow(dietterich_p_blue[i] - color_avg, 2));
        }

        //******************** ORIGINAL DIETTERICH FORMULA
        double dietterich_numerator = dietterich_p_blue[0];
        //********************
        //******************** NEW FORMULA
        //double dietterich_numerator = 0;

        //std::uniform_int_distribution<int> kgen(0, k_fold-1);
        //int random_roll = kgen(gen);
        //dietterich_numerator += dietterich_p_blue[random_roll];

        //********************


        double dietterich_denominator = 0;

        for (int i = 0; i < k_fold; ++i)
        {
            dietterich_denominator += dietterich_s[i];
               
        }

        double dietterich_formula = 0;

		if (dietterich_denominator != 0)
		{
			dietterich_formula = dietterich_numerator / std::sqrt(dietterich_denominator/k_fold);
		}

        test_statistics.push_back(dietterich_formula);
    }


    if (print_it == true)
    {

        std::ofstream file_out;
        std::string file_name = "CV_Dietterich.csv";
        file_out.open(file_name);
        if (!file_out.is_open())
        {
            Rcpp::stop("Could not open output file.");
        }    

        file_out << "Dietterich test" << std::endl << std::endl;

        file_out << "k fold: " << k_fold << std::endl;
        file_out << "leave m out: " << leave_m_out << std::endl;


        file_out << "\nIndices of the removed rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < leave_m_out; ++j)
            {
                file_out << random_row_indices[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nColumn IDs after reordering\n";
        for (int i = 0; i < columns_ordered.size(); ++i)
        {
            file_out << columns_ordered[i] << "; ";
        }
        file_out << std::endl;

        file_out << "\nSRD values of selected rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < data_matrix[0].size() - 1; ++j)
            {
                file_out << cv_srd_blue[i][j] << "; ";
            }

            file_out << std::endl;
        }


        file_out << "\nSRD values of left out rows\n";
        for (int i = 0; i < k_fold; ++i)
        {
            for (int j = 0; j < data_matrix[0].size() - 1; ++j)
            {
                file_out << cv_srd_orange[i][j] << "; ";
            }

            file_out << std::endl;
        }

        file_out << "\nBoxplot\n";
        for (int i = 0; i < 7; ++i)
        {
            for (int k = 0; k < cv_srd[0].size(); ++k)
            {
                if (i == 0 && k == 0)
                    file_out << "min; ";
                if (i == 1 && k == 0)
                    file_out << "xx1; ";
                if (i == 2 && k == 0)
                    file_out << "q1; ";
                if (i == 3 && k == 0)
                    file_out << "median; ";
                if (i == 4 && k == 0)
                    file_out << "q3; ";
                if (i == 5 && k == 0)
                    file_out << "xx19; ";
                if (i == 6 && k == 0)
                    file_out << "max; ";

                file_out << boxplot_stats[i][k] << "; ";
            }
            file_out << std::endl;
        }


        file_out << "\nTest statistics\n";
        for (int i = 0; i < test_statistics.size(); ++i)
        {
            file_out << test_statistics[i] << "; ";
        }
		
		file_out << "\n\nStatistical significance\n";

        std::vector<double> p_val_005 = {2.571,2.447,2.365,2.306,2.262,2.228};
        std::vector<double> p_val_010 = {2.015,1.943,1.895,1.860,1.833,1.812};

        for (int i = 0; i < test_statistics.size(); ++i)
        {
            if (k_fold < 5 || k_fold>10)
            {
                file_out << "n.a.; ";
            }
            else if (abs(test_statistics[i]) >= p_val_005[k_fold - 5])
            {
                file_out << "(p<0.05*); ";
            }
            else if (abs(test_statistics[i]) >= p_val_010[k_fold - 5])
            {
                file_out << "(p<0.1); ";
            }
            else
            {
                file_out << "n.s.; ";
            }

        }

    }

    return columns_ordered;
}
