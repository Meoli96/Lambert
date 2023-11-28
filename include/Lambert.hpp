

#include <cmath>
#include <span>

#include "Utils.hpp"

double LambertSolve(double r1[6], double r2[6], double dt, bool prograde = 1) {
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
            v1[i] = (r2[i] - f * r1[i]) / g;  // Remove Earth velociy
            v2[i] =
                (g_dot * r2[i] - r1[i]) / g;  // Remove this from Mars velocity

            for (size_t i = 0; i < 3; ++i) {
                v1[i] -= r1[i + 3];
                v2[i] = r2[i + 3] - v2[i];
            }
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

int computeDV(size_t i, size_t j, std::span<SpiceDouble> dep_days,
              std::span<SpiceDouble> arr_days, double& res_v, double& res_m) {
    // Compute state of Earth at departure
    SpiceDouble state_earth[6];
    spkezr_c("EARTH", arr_days[j], "J2000", "NONE", "SUN", state_earth,
             nullptr);

    // Compute state of Venus at arrival

    SpiceDouble state_venus[6];
    spkezr_c("VENUS", arr_days[j], "J2000", "NONE", "SUN", state_venus,
             nullptr);

    // Compute state of Mars at arrival
    SpiceDouble state_mars[6];
    spkezr_c("MARS", arr_days[j], "J2000", "NONE", "SUN", state_mars, nullptr);

    double dt = arr_days[j] - dep_days[i];

    // Lambert - Get deltaV

    double deltaV_v = LambertSolve(state_earth, state_venus, dt);
    double deltaV_m = LambertSolve(state_earth, state_mars, dt);

    // Compute deltaV

    // Store deltaV in res_mat[i][j] - 7.5 is the max deltaV
    if (deltaV_v > 7.5) {
        deltaV_v = 0;
    }
    if (deltaV_m > 7.5) {
        deltaV_m = 0;
    }
    res_v = deltaV_v;
    res_m = deltaV_m;
    

    return 0;
}