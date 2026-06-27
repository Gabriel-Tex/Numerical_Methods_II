#include "../include/function.h"

double F(double x) {
    return pow(x, 3);
}

double F_NC(double x) {
    return pow(sin(2.0 * x) + 4.0 * x * x + 3.0 * x, 2);
}

double F_sing1(double x) {
    return pow(x, -1.0 / 3.0);
}

double F_sing2(double x) {
    return pow(x, -0.5) * exp(-x);
}

double F_surf(double x, double y) {
    return 0.2 * (x * x - y * y);
}

double dF_surf_dx(double x, double y) {
    (void)y;
    return 0.4 * x;
}

double dF_surf_dy(double x, double y) {
    (void)x;
    return -0.4 * y;
}