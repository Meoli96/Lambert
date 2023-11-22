#include <IO>
#include <cmath>
#include <span>

#include "Utils.hpp"

double LambertSolve(double r1[3], double r2[3], double dt, bool prograde = 1) {
    double alpha, beta;

    double r1_norm = sqrt(pow(r1[0], 2) + pow(r1[1], 2) +
                          pow(r1[2], 2));  // Distance between r1 and origin
    double r2_norm = sqrt(pow(r2[0], 2) + pow(r2[1], 2) +
                          pow(r2[2], 2));  // Distance between r2 and origin

    double r_cross[3];
    cross(r1, r2, r_cross);  // Cross product between r1 and r2
    double r_cross_z = r_cross[2];
    // Compute dtheta based on prograde/retrograde and r_cross[2] sign
    double dtheta;
    if (prograde) {
        // Prograde
        if (r_cross_z >= 0) {
            dtheta = acos(dot(r1, r2) / (r1_norm * r2_norm));
        } else {
            // r_cross_z < 0
            dtheta = 2 * M_PI - acos(dot(r1, r2) / (r1_norm * r2_norm));
        }

    } else {
        // Retrograde
        if (r_cross_z >= 0) {
            dtheta = 2 * M_PI - acos(dot(r1, r2) / (r1_norm * r2_norm));
        } else {
            // r_cross_z < 0
            dtheta = acos(dot(r1, r2) / (r1_norm * r2_norm));
        }
    }
    double dtheta_norm =
        normalize(dtheta, 0, 2 * M_PI);  // Normalize dtheta between 0 and 2pi

    double A =
        sin(dtheta_norm) * sqrt(r1_norm * r2_norm / (1 - cos(dtheta_norm)));

    // Determine aproximatively where F(z,t) changes sign
    double z0 = -100;
    while (F(z0, dt, A, r1_norm, r2_norm) < 0) {
        z0 += 0.1;
    }
    double tol = 1e-8;
    double args[4] = {A, r1_norm, r2_norm, dt};

    double z_n = newton_iter(F, dFdz, z0, tol, 1000, args);

    if (z_n > 0) {
        // Elliptic orbit

        // Compute Lagrange coefficients
        double f = 1 - y(z_n, A, r1_norm, r2_norm) / r1_norm;
        double g = A * sqrt(y(z_n, A, r1_norm, r2_norm) / mu_s);
        double g_dot = 1 - y(z_n, A, r1_norm, r2_norm) / r2_norm;

        // Compute velocity at r1 and r2
        double v1[3], v2[3];
        for (size_t i = 0; i < 3; ++i) {
            v1[i] = (r2[i] - f * r1[i]) / g;
            v2[i] = (g_dot * r2[i] - r1[i]) / g;
        }

        // Or i could return abs(v2-v1) as deltaV
        double dV = sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) +
                         pow(v1[2] - v2[2], 2));
        return dV;

    } else {
        // Hyperbolic or parabolic orbit, return 0
        return 0;
    }
}

int computeDV(size_t i, size_t j, std::span<uint> dep_days,
              std::span<uint> arr_days, double** res_v, double** res_m) {
    // Compute state of Earth at departure
    IO::Astrodynamics::State earth_state = IO::Astrodynamics::GetState(
        "Earth", dep_days[i], IO::Astrodynamics::ReferenceFrame::J2000);

    const double r_earth[3] = {state_earth.position[0], state_earth.position[1],
                               state_earth.position[2]};

    // Compute state of Venus at arrival

    IO::Astrodynamics::State state_venus = IO::Astrodynamics::GetState(
        "Venus", arr_days[j], IO::Astrodynamics::ReferenceFrame::J2000);

    IO::Astrodynamics::State state_mars = IO::Astrodynamics::GetState(
        "Mars", arr_days[j], IO::Astrodynamics::ReferenceFrame::J2000);

    double dt = arr_days[j] - dep_days[i];

    // Lambert - Get deltaV

    double r_venus[3] = {state_venus.position[0], state_venus.position[1],
                         state_venus.position[2]};

    double r_mars[3] = {state_mars.position[0], state_mars.position[1],
                        state_mars.position[2]};

    double deltaV_v = LambertSolve(r_earth, r_venus, dt);
    double deltaV_m = LambertSolve(r_earth, r_mars, dt);

    // Compute deltaV

    // Store deltaV in res_mat[i][j] - 7.5 is the max deltaV
    if (deltaV_v > 7.5) {
        deltaV_v = 0;
    }
    if (deltaV_m > 7.5) {
        deltaV_m = 0;
    }
    res_mat_v[i][j] = deltaV_v;
    res_mat_m[i][j] = deltaV_m;
}