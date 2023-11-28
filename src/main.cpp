#include <SpiceUsr.h>
#include <matplot/matplot.h>

#include <Lambert.hpp>
#include <ThreadPool.hpp>
#include <Utils.hpp>
#include <iostream>
#include <span>
#include <string>
#include <thread>
#include <vector>

const uint days_in_seconds = 86400;

const uint step = 1;  // Step in days
const char cspice_lsk[] = "kernel/lsk/naif0012.tls";
const char cspice_pck[] = "kernel/pck/pck00010.tpc";
const char cspice_spk[] = "kernel/spk/de432s.bsp";

int main() {
    // Load kernel
    furnsh_c(cspice_lsk);  // Leap seconds kernel
    furnsh_c(cspice_spk);  // Planetary ephemeris kernel
    furnsh_c(cspice_pck);  // Leap seconds kernel

    // Define vector of departure days, 1 day step for 365 days
    std::vector<SpiceDouble> departure_days(365 * 4);
    // Define vector of arrival days, from 3 months to 4 years from departure, 1
    // day step
    std::vector<SpiceDouble> arrival_days(365 * 4 - 90);
    std::span<SpiceDouble> dep_span(departure_days);
    std::span<SpiceDouble> trav_span(arrival_days);

    // Fill departure days vector starting from 1 September 2030

    std::string dep_date = "2030-09-01";
    SpiceDouble et_dep;
    str2et_c(dep_date.c_str(), &et_dep);
    for (size_t i = 0; i < departure_days.size(); ++i) {
        departure_days[i] = et_dep + i * days_in_seconds;
    }

    // Fill travel times vector starting from 3 months after 1 September
    // 2030
    std::string arr_date = "2030-12-01";
    SpiceDouble et_arr;
    str2et_c(arr_date.c_str(), &et_arr);

    for (size_t i = 0; i < arrival_days.size(); ++i) {
        arrival_days[i] = et_arr + i * days_in_seconds;
    }
    double** res_mat_v = new double*[departure_days.size()];
    double** res_mat_m = new double*[departure_days.size()];

    for (size_t i = 0; i < departure_days.size(); ++i) {
        res_mat_v[i] = new double[arrival_days.size()];
        res_mat_m[i] = new double[arrival_days.size()];
    }
    // i is departure day, j is travel time
    using namespace thread;
    // Create thread pool
    ThreadPool<int(size_t, size_t, std::span<SpiceDouble>,
                   std::span<SpiceDouble>, double&, double&)>
        pool(computeDV);

    // Add jobs to thread pool
    for (size_t i = 0; i < departure_days.size(); ++i) {
        for (size_t j = 0; j < arrival_days.size(); ++j) {
            pool.addjob(i, j, dep_span, trav_span, res_mat_v[i][j],
                        res_mat_m[i][j]);
        }
    }
    // Wait for jobs to finish
    pool.wait();
    std::cout << "Done!" << std::endl;

    // // Plot results
    // matplot::figure();
    // matplot::contour(res_mat_v, departure_days.size(), arrival_days.size());
    // matplot::title("DeltaV to Venus");
    // matplot::xlabel("Departure day");
    // matplot::ylabel("Travel time");
    // matplot::save("deltaV_v.png");

    // matplot::figure();
    // matplot::contour(res_mat_m, departure_days.size(), arrival_days.size());
    // matplot::title("DeltaV to Mars");
    // matplot::xlabel("Departure day");
    // matplot::ylabel("Travel time");
    // matplot::save("deltaV_m.png");

    return 0;
}
