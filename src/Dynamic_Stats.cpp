#include "Dynamic_Stats.h"
#include <cmath>

Dynamic_Stats::Dynamic_Stats()
{
	round = 1;
	mean_new = 0;
	std_dev_new = 0;
	mean_old = 0;
	std_dev_old = 0;
	data = 0;
}


// This dynamic mean and std deviation computation is superior to storing the values and applying a formula. 
// It goes back to a 1962 paper by B. P. Welford and is presented in Donald Knuth's Art of Computer Programming, 
// Vol 2, page 232, 3rd edition. https://www.johndcook.com/blog/standard_deviation/

void Dynamic_Stats::push_data(double data) 
{
	

	if (round == 1)
	{
		mean_new = data;
		std_dev_new = 0;
		mean_old = mean_new;
		std_dev_old = std_dev_new;
	}
	else
	{
		mean_new = mean_old + (data - mean_old) / round;
		std_dev_new = std_dev_old + (data - mean_old) * (data - mean_new);
		mean_old = mean_new;
		std_dev_old = std_dev_new;
	}

	++round;
}

double Dynamic_Stats::get_avg()
{
	return mean_new;
}

double Dynamic_Stats::get_std_dev()
{
	if(round - 2 == 0)
	{
		return 0;
	}
	else
	{
	return std::sqrt((std_dev_new) / (round - 2));
	}
	
}
