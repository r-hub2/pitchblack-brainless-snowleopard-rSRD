#ifndef DYNAMIC_STATS_H
#define DYNAMIC_STATS_H

#include <iostream>

class Dynamic_Stats {
public:
    Dynamic_Stats();
    void push_data(double data);
    double get_avg();
    double get_std_dev();

private:
    int round;
    double mean_old, mean_new, std_dev_old, std_dev_new, data;
};

#endif /* DYNAMIC_STATS_H */
