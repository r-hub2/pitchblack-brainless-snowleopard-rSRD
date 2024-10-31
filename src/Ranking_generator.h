#ifndef RANKING_GENERATOR_H
#define RANKING_GENERATOR_H

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>



class Ranking_generator
{

public:


	static std::vector<double> randomRanking(int size_of_ranking, double tie_probability);
	static std::vector<double> randomRanking_noties(int size_of_ranking);

	static std::vector<double> perturbation_simp_inv(std::vector<double> ranking, int number_of_inversions);
	static std::vector<double> perturbation_first_half(std::vector<double> ranking, int number_of_inversions);
	static std::vector<double> perturbation_second_half(std::vector<double> ranking, int number_of_inversions);
	static std::vector<double> perturbation_med_size(std::vector<double> ranking, int number_of_inversions);
	static std::vector<double> perturbation_underdog(std::vector<double> ranking);

};

#endif /* RANKING_GENERATOR_H */
