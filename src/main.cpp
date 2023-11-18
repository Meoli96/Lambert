#include <iostream>
#include <span>
#include <thread>
#include <vector>
#include <IO>
// #include io-aerospace
// #include matplotlib-cpp

const uint days_in_seconds = 86400;

const uint step = 1;  // Step in days



int main() {

    // Load kernel
    IO::Astrodynamics::Kernels::KernelsLoader::Load("Data/SolarSystem");
    // Define vector of departure days, 1 day step for 400 days
    std::vector<uint> departure_days(400);
    // Define vector of travel times, from 3 months to 4 years, 1 day step
    std::vector<uint> travel_times(365 * 4 - 90);

    std::span dep_span(departure_days);
    std::span trav_span(travel_times);


    // i is departure day, j is travel time
    double res_mat[departure_days.size()][travel_times.size()];

    for (size_t i = 0; i < departure_days.size(); ++i) {
        // Get state of Earth at departure

        /* 
            Posso fare che ogni thread si prende un suo giorno di partenza e si itera
            in parallelo su tutti i tempi di viaggio -> gestione della memoria senza lock
            
            Oppure fare tipo un threadpool che prende una coppia i,j e va per se -> servono lock?
            */
        for (size_t j = 0; j < travel_times.size(); ++j) {
            // Get state of Venus at departure + travel time

            // Get state of Mars at departure + travel time

            // Lambert - Get deltaV
            // Store deltaV in res_mat[i][j]

        }
    }
    return 0;
}
