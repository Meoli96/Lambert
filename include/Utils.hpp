#pragma once
#include <cassert>
#include <cmath>
#include <functional>
double mu_s =
    1.372 * pow(10, 11);  // gravitational parameter of the Sun - km^3/s^2

double dot(double a[3], double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double sign(double x) {
    if (x == 0)
        return 0;
    else {
        return (x > 0) - (x < 0);
    }
}
double normalize(const double value, const double start, const double end) {
    const double width = end - start;          //
    const double offsetValue = value - start;  // value relative to 0

    return (offsetValue - (floor(offsetValue / width) * width)) + start;
    // + start to reset back to start of original range
}

double* cross(double a[3], double b[3], double* res) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    return res;
}

double newton_iter(
    std::function<double(double, double, double, double, double)> f,
    std::function<double(double, double, double, double, double)> df, double x0,
    double tol = 1e-6, int max_iter = 1000, double* args = nullptr) {
    if (args != nullptr) {
        double x = x0;

      
        double t = args[0];
        double A = args[1];
        double r1_mod = args[2];
        double r2_mod = args[3];

        double fx = f(x, t, A, r1_mod, r2_mod);
        int i = 0;
        while (abs(fx) > tol && i < max_iter) {
            x = x - fx / df(x, t, A, r1_mod, r2_mod);
            fx = f(x, t, A, r1_mod, r2_mod);
            i++;
        }
        return x;
    }
    return 0;
}

double S(double z) {
    if (z > 0) {
        return (sqrt(z) - sin(sqrt(z))) / pow(z, 1.5);
    } else if (z < 0) {
        return (sinh(sqrt(-z)) - sqrt(-z)) / pow(-z, 1.5);
    } else {
        return 1.0 / 6.0;
    }
}

double C(double z) {
    if (z > 0) {
        return (1 - cos(sqrt(z))) / z;
    } else if (z < 0) {
        return (cosh(sqrt(-z)) - 1) / (-z);
    } else {
        return 0.5;
    }
}

double y(double z, double A, double r1_mod, double r2_mod) {
    return r1_mod + r2_mod + A * (z * S(z) - 1) / sqrt(C(z));
}

double F(double z, double t, double A, double r1_mod, double r2_mod) {
    return pow(y(z, A, r1_mod, r2_mod) / C(z), 1.5) * S(z) +
           A * sqrt(y(z, A, r1_mod, r2_mod)) - sqrt(mu_s) * t;
}

double dFdz(double z, double t, double A, double r1_mod, double r2_mod) {
    if (z == 0) {
        return sqrt(2) / 40 * pow(y(0, A, r1_mod, r2_mod), 1.5) +
               A / 8 * sqrt(y(0, A, r1_mod, r2_mod)) +
               A * sqrt(0.5 / y(0, A, r1_mod, r2_mod));
    } else {
        return pow(y(z, A, r1_mod, r2_mod) / C(z), 1.5) *
                   ((0.5 / z) * (C(z) - (3 * S(z)) / (2 * C(z))) +
                    (3 * S(z) * S(z)) / (4 * C(z))) +
               A / 8 *
                   ((3 * S(z) / C(z)) * sqrt(y(z, A, r1_mod, r2_mod)) +
                    A * sqrt(C(z) / y(z, A, r1_mod, r2_mod)));
    }
}

double f(double z, double t, double A, double r1_mod, double r2_mod) {
    return 1 - (y(z, A, r1_mod, r2_mod) / r1_mod);
}
double g(double z, double t, double A, double r1_mod, double r2_mod) {
    return A * sqrt(y(z, A, r1_mod, r2_mod) / mu_s);
}
double f_dot(double z, double t, double A, double r1_mod, double r2_mod) {
    return sqrt(y(z, A, r1_mod, r2_mod) * mu_s / C(z)) *
           ((z * S(z) - 1) / (r1_mod * r2_mod));
}
double g_dot(double z, double t, double A, double r1_mod, double r2_mod) {
    return 1 - (y(z, A, r1_mod, r2_mod) / r2_mod);
}