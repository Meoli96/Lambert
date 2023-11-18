#include <cmath>

#include "Utils.hpp"

void LambertSolve(double r1[3], double r2[3], double dt, bool prograde = 1) {
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

    double c =
        sqrt(pow(r2[0] - r1[0], 2) + pow(r2[1] - r1[1], 2) +
             pow(r2[2] - r1[2], 2));  // Distance between r1 and r2 - chord
    double sp = (r1_norm + r2_norm + c) / 2;  // Semiperimeter

    double t_hyp = (sqrt(2 / mu_s) / 3) *
                   (pow(sp, 3 / 2) * sign(sin(dtheta)) *
                    pow(sp - c, 3 / 2));  // Parabolic time of flight

    if (dt < t_hyp) {
        // Hyperbolic case, return
        return;
    } else {
        // Elliptic case
        // Compute beta_m, t_m of minimal energy semimajor axis
        double sma_m = sp / 2;
        double beta_m;
        if (dtheta_norm < M_PI) {
            beta_m = acos(1 - (sp - c) / sma_m);
        } else if (dtheta_norm > M_PI) {
            beta_m = -acos(1 - (sp - c) / sma_m);
        } else {
            // Throw numerical error

            beta_m = 0;
        }
        double t_m =
            sqrt(pow(sma_m, 3) / (8 * mu_s)) * (M_PI - beta_m + sin(beta_m));
    }
    // compute M(  )
}