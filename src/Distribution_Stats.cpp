#include "Distribution_Stats.h"
#include <cmath>

Distribution_Stats::Distribution_Stats()
{
	SRD_distribution = { {0}, {0} };
}


void Distribution_Stats::push_data(double data)
{
	data = double(round(10000 * data)) / double(10000); // Round data to 2 decimal point precision
	bool new_data_inserted = false;
	std::vector<double>::iterator iter0 = SRD_distribution[0].begin();
	std::vector<double>::iterator iter1 = SRD_distribution[1].begin();

	for (int i = 0; i < SRD_distribution[0].size(); ++i)
	{
		if (data == SRD_distribution[0][i])
		{
			++SRD_distribution[1][i];
			new_data_inserted = true;
			break;
		}
		else if (data < SRD_distribution[0][i])
		{
			SRD_distribution[0].insert(iter0 + i, data);
			SRD_distribution[1].insert(iter1 + i, 1);
			new_data_inserted = true;
			break;
		}
	}



	if(new_data_inserted == false)
	{
		SRD_distribution[0].push_back(data);
		SRD_distribution[1].push_back(1);
	}

}



double Distribution_Stats::get_xx1(double sim_size)
{
	double cumulated_frequency = 0;
	for (int i = 0; i < SRD_distribution[1].size(); ++i)
	{
		cumulated_frequency += SRD_distribution[1][i]/sim_size;
		if (cumulated_frequency >= 0.05)
		{
			return SRD_distribution[0][i];
		}
	}

	return -1;
}


double Distribution_Stats::get_xx19(double sim_size)
{
	double cumulated_frequency = 0;
	for (int i = 0; i < SRD_distribution[1].size(); ++i)
	{
		cumulated_frequency += SRD_distribution[1][i] / sim_size;
		if (cumulated_frequency >= 0.95)
		{
			return SRD_distribution[0][i];
		}
	}

	return -1;
}


double Distribution_Stats::get_q1(double sim_size)
{
	double cumulated_frequency = 0;
	for (int i = 0; i < SRD_distribution[1].size(); ++i)
	{
		cumulated_frequency += SRD_distribution[1][i] / sim_size;
		if (cumulated_frequency >= 0.25)
		{
			return SRD_distribution[0][i];
		}
	}

	return -1;
}


double Distribution_Stats::get_q3(double sim_size)
{
	double cumulated_frequency = 0;
	for (int i = 0; i < SRD_distribution[1].size(); ++i)
	{
		cumulated_frequency += SRD_distribution[1][i] / sim_size;
		if (cumulated_frequency >= 0.75)
		{
			return SRD_distribution[0][i];
		}
	}

	return -1;
}



double Distribution_Stats::get_median(double sim_size)
{
	double cumulated_frequency = 0;
	for (int i = 0; i < SRD_distribution[1].size(); ++i)
	{
		cumulated_frequency += SRD_distribution[1][i] / sim_size;
		if (cumulated_frequency >= 0.5)
		{
			return SRD_distribution[0][i];
		}
	}

	return -1;
}

double Distribution_Stats::get_min()
{
	if (SRD_distribution[1][0] > 0)
	{
		return SRD_distribution[0][0];
	}
	else
	{
		return SRD_distribution[0][1];
	}
	
	return -1;
}


double Distribution_Stats::get_max()
{
	return SRD_distribution[0][SRD_distribution[0].size() - 1];
}


std::vector<std::vector<double>>Distribution_Stats::get_distribution()
{
	return SRD_distribution;
}
