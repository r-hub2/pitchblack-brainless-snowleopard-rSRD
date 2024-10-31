#include "Manhattan_Dist.h"
#include <stdlib.h>
#include <cmath>

double Manhattan_Dist::manhattan(std::vector<double> first, std::vector<double> second)
{

    double nSRD=0;
    for (int i = 0; i < first.size(); i++)
    {
        nSRD += std::abs(first[i] - second[i]);
    }

    if (first.size() % 2 == 0) //We normalize by the max SRD value which depends on the number of rows (there are two formulas depending on the parity)
    {
        nSRD /= 2 * std::pow((first.size() / 2), 2);
    }
    else
    {
        nSRD /= (first.size()-1) * (first.size() + 1)/2;
    }
	return nSRD;
}
