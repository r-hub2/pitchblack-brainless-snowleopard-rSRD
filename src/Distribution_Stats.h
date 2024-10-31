#ifndef DISTRIBUTION_STATS_H
#define DISTRIBUTION_STATS_H

#include <iostream>
#include <vector>


class Distribution_Stats {
public:
    Distribution_Stats();
    void push_data(double data);
    std::vector<std::vector<double>>get_distribution();
    double get_xx1(double sim_size);
    double get_xx19(double sim_size);
    double get_q1(double sim_size);
    double get_q3(double sim_size);
    double get_median(double sim_size);
    double get_min();
    double get_max();

private:
    std::vector<std::vector<double>>SRD_distribution;
};


#endif /* DISTRIBUTION_STATS_H */
