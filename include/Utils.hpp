#pragma once
#include <cmath>
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