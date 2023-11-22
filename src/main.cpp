#include <KernelsLoader.hpp>
#include <ThreadPool.hpp>

#include <iostream>
#include <span>
#include <thread>
#include <vector>
// #include io-aerospace
// #include matplotlib-cpp

const uint days_in_seconds = 86400;

const uint step = 1;  // Step in days

int main() {
    // Load kernel
    IO::Astrodynamics::Kernels::KernelsLoader::Load("Data/SolarSystem");
    // Define vector of departure days, 1 day step for 365 days
    std::vector<uint> departure_days(365);
    // Define vector of arrival days, from 3 months to 4 years from departure, 1
    // day step
    std::vector<uint> arrival_days(365 * 4 - 90);

    std::span dep_span(departure_days);
    std::span trav_span(arrival_days);

// Fill departure days vector starting from 1 September 2030
std:
    string date = "2030-09-01";
    for (size_t i = 0; i < departure_days.size(); ++i) {
        departure_days[i] = IO::Time::DateToSeconds(date);
        date = IO::Time::SecondsToDate(departure_days[i] + days_in_seconds);
    }
    Janu

        // Fill travel times vector starting from 3 months after 1 September
        // 2030
        date = "2030-12-01";
    for (size_t i = 0; i < arrival_days.size(); ++i) {
        arrival_days[i] = IO::Time::DateToSeconds(date);
        date = IO::Time::SecondsToDate(arrival_days[i] + days_in_seconds);
    }

    // i is departure day, j is travel time
    double res_mat_v[departure_days.size()][arrival_days.size()];
    double res_mat_m[departure_days.size()][arrival_days.size()];

    // Create thread pool
    ThreadPool pool(6);

    // Add jobs to thread pool
    for (size_t i = 0; i < departure_days.size(); ++i) {
        for (size_t j = 0; j < arrival_days.size(); ++j) {
            pool.addjob(computeDV, i, j, dep_span, trav_span, res_mat_v,
                        res_mat_m);
        }
    }
    // Wait for jobs to finish
    pool.wait();
    std::cout << "Done!" << std::endl;


    

   
    // join threads here

    // Plot results with matplotlib-cpp
    // std::vector<double> x_axis;
    // std::vector<double> y_axis;
    // for (size_t i = 0; i < departure_days.size(); ++i) {
    //     for (size_t j = 0; j < arrival_days.size(); ++j) {
    //         x_axis.push_back(departure_days[i]);
    //         y_axis.push_back(arrival_days[j]);
    // add dv as gradient of color

    //     }
    // }
    // plt::plot(x_axis, y_axis, "o");
    // plt::show();

    return 0;
}
