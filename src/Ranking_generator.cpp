#include "Ranking_generator.h"
#include <algorithm>


std::vector<double> Ranking_generator::randomRanking(int size_of_ranking, double tie_probability)
{

    std::vector<double>ranking; //We initialize the ranking without any randomness or ties
    for (int i = 1; i <= size_of_ranking; i++)
    {
        ranking.push_back(i);
    }

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> rtie(0, 1); //Generating random booleans (0/1)

    std::vector<int>tied_ranks; //Shows which ranks are tied. Contains 'size_of_ranking-1' boolean (0/1) values

    for (int i = 1; i <= size_of_ranking - 1; i++)
    {
        if (rtie(gen) <= tie_probability)
        {
            tied_ranks.push_back(1);
        }
        else
        {
            tied_ranks.push_back(0);
        }
    }

    tied_ranks.push_back(0); //We need to end the sequence with a 0 as the last element cannot be tied with the next (as there is no next element)

    int tie_size = 0; //Shows the size of a tie block. If 'tied_ranks' contains more than one consequtive 1's it means that there is a threeway or even larger tie block. 
    for (int i = 0; i < size_of_ranking; i++)
    {

        if (tied_ranks[i] == 1)
        {
            tie_size++; //With each 1 we increase the size of the tie block
        }
        else if (tie_size != 0) //At the end of a tie block we adjust the ranking. E.g. if 'tied_ranks' starts with two 1's followed by a 0, that means that the first three elements of 'ranking' are tied.
        {
            double new_rank = (ranking[i] + ranking[i - tie_size]) / 2;
            for (int j = i - tie_size; j <= i; j++)
            {
                ranking[j] = new_rank; //Tied values are registered
            }
            tie_size = 0;
        }

    }


    std::shuffle(ranking.begin(), ranking.end(), gen); //Random shuffling of the vector elements

    return ranking;
}


std::vector<double> Ranking_generator::randomRanking_noties(int size_of_ranking)
{

    std::vector<double> ranking; //We initialize the ranking without any randomness or ties
    for (int i = 1; i <= size_of_ranking; i++)
    {
        ranking.push_back(i);
    }

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::shuffle(ranking.begin(), ranking.end(), gen); //Random shuffling of the vector elements

    return ranking;
}



std::vector<double> Ranking_generator::perturbation_simp_inv(std::vector<double> ranking, int number_of_inversions)
{


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::uniform_int_distribution<> deg(0, ranking.size()-2);

    for (int i = 0; i < number_of_inversions; ++i)
    {
        int index = deg(gen);
        std::swap(ranking[index], ranking[index + 1]);
    }


    return ranking;
}


std::vector<double> Ranking_generator::perturbation_first_half(std::vector<double> ranking, int number_of_inversions)
{


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::uniform_int_distribution<> deg(0, std::floor(ranking.size()/2) - 2);

    for (int i = 0; i < number_of_inversions; ++i)
    {
        int index = deg(gen);
        std::swap(ranking[index], ranking[index + 1]);
    }


    return ranking;
}


std::vector<double> Ranking_generator::perturbation_second_half(std::vector<double> ranking, int number_of_inversions)
{


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::uniform_int_distribution<> deg(std::floor(ranking.size() / 2), ranking.size() - 2);

    for (int i = 0; i < number_of_inversions; ++i)
    {
        int index = deg(gen);
        std::swap(ranking[index], ranking[index + 1]);
    }


    return ranking;
}

std::vector<double> Ranking_generator::perturbation_med_size(std::vector<double> ranking, int number_of_inversions)
{


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::uniform_int_distribution<> deg(0, ranking.size() - std::ceil(double(ranking.size()) / 4)-1);

    for (int i = 0; i < number_of_inversions; ++i)
    {
        int index = deg(gen);
        std::swap(ranking[index], ranking[index + ceil(double(ranking.size()) / 4)]);
    }


    return ranking;
}


std::vector<double> Ranking_generator::perturbation_underdog(std::vector<double> ranking)
{


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::uniform_int_distribution<> deg(std::floor(ranking.size() / 2), ranking.size() - 1);

    int index = deg(gen);
    std::swap(ranking[0], ranking[index]);

    return ranking;
}
