#include <iostream>
#include <span>
#include <thread>
#include <vector>
// #include io-aerospace
// #include matplotlib-cpp

const uint days_in_seconds = 86400;

const uint step = 1;  // Step in days

int main() {
    // Define vector of departure days, 1 day step for 400 days
    std::vector<uint> departure_days(400);
    // Define vector of travel times, from 3 months to 4 years, 1 day step
    std::vector<uint> travel_times(365 * 4 - 90);

    // for each departure day
    //    for each travel time
    //       State of Earth(r1)
    //       State of Venus(r2)
    //       State of Mars(r3)

    return 0;
}
