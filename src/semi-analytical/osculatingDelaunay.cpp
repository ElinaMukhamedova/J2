#include "osculatingDelaunay.hpp"
#include <cmath>


double l_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                        double f_prime, double alpha, double mu) {
    double const eta_prime = G_prime / L_prime;
    double const e_prime = std::sqrt(1.0 - G_prime * G_prime / L_prime / L_prime);
    double const s2 = 1.0 - H_prime * H_prime / G_prime / G_prime;
    double const factor = alpha * alpha * mu * mu / 8.0 / G_prime / G_prime / G_prime;
    double const partial_e_prime_L_prime = eta_prime * eta_prime / L_prime / e_prime;
    double const partial_eta_prime_L_prime = -eta_prime/L_prime;

    double const first = 2.0 * std::sin(f_prime) * (3.0 * s2 - 2.0) * partial_e_prime_L_prime;
    double const second = -3.0 * s2 * std::sin(f_prime + 2.0 * g_prime) * partial_e_prime_L_prime;
    double const third = -s2 * std::sin(3.0 * f_prime + 2.0 * g_prime) * partial_e_prime_L_prime;
    double const fourth = -s2 * std::sin(2.0 * g_prime) * (2.0 * eta_prime * eta_prime + 4.0 * eta_prime) / (eta_prime + 1.0) / (eta_prime + 1.0) * partial_eta_prime_L_prime;

    return factor * (first + second + third + fourth);
}

double g_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                    double f_prime, double alpha, double mu) {
    double const eta_prime = G_prime / L_prime;
    double const e_prime = std::sqrt(1.0 - G_prime * G_prime / L_prime / L_prime);
    double const s2 = 1.0 - H_prime * H_prime / G_prime / G_prime;

    double const first = 2.0 * (f_prime - l_prime) * alpha * alpha * mu * mu *
                            (-3.0 / 8.0 / G_prime / G_prime / G_prime / G_prime +
                            15.0 * H_prime * H_prime / 8.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime);
    double const second = -alpha * alpha * mu * mu * (3.0*s2 - 2.0) * std::sin(f_prime) / 4.0 / G_prime / G_prime / L_prime / L_prime / e_prime
                        -3.0 * alpha * alpha * mu * mu * (3.0*s2 - 2.0) * e_prime * std::sin(f_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime
                        +3.0 * H_prime + H_prime * alpha * alpha * mu * mu * e_prime * std::sin(f_prime) / 2.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const third = 9.0 * alpha * alpha * mu * mu * s2 * std::sin(2.0*f_prime+2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime
                        -3.0 * H_prime * H_prime * alpha * alpha * mu * mu * std::sin(2.0*f_prime+2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const fourth = 3.0 * alpha * alpha * mu * mu * s2 * std::sin(f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / L_prime / L_prime / e_prime
                        +9.0 * alpha * alpha * mu * mu * s2 * e_prime * std::sin(f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime
                        -3.0 * H_prime * H_prime * alpha * alpha * mu * mu * e_prime * std::sin(f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const fifth = alpha * alpha * mu * mu * s2 * std::sin(3.0*f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / L_prime / L_prime / e_prime
                        +3.0 * alpha * alpha * mu * mu * s2 * e_prime * std::sin(3.0*f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime
                        -H_prime * H_prime * alpha * alpha * mu * mu * e_prime * std::sin(3.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const sixth = alpha * alpha * mu * mu * s2 * (1.0 - eta_prime) * std::sin(2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / L_prime / (1.0 + eta_prime)
                        -alpha * alpha * mu * mu * s2 * (1.0 - eta_prime) * (1.0 + 2.0*eta_prime) * std::sin(2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / L_prime / (1.0 + eta_prime) / (1.0 + eta_prime)
                        -alpha * alpha * mu * mu * s2 * (1.0 + 2.0*eta_prime) * std::sin(2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / L_prime / (1.0 + eta_prime)
                        -3.0 * alpha * alpha * mu * mu * s2 * (1.0 - eta_prime) * (1.0 + 2.0*eta_prime) * std::sin(2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime / (1.0 + eta_prime)
                        +H_prime * H_prime * alpha * alpha * mu * mu * (1.0 - eta_prime) * (1.0 + 2.0*eta_prime) * std::sin(2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime / G_prime / (1.0 + eta_prime);
    return first+second+third+fourth+fifth+sixth;
}

double h_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                    double f_prime, double alpha, double mu) {
    double const eta_prime = G_prime / L_prime;
    double const e_prime = std::sqrt(1.0 - G_prime * G_prime / L_prime / L_prime);
    double const s2 = 1.0 - H_prime * H_prime / G_prime / G_prime;

    double const first = -3.0 * H_prime * alpha * alpha * mu * mu * (f_prime - l_prime) / 2.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const second = -3.0 * H_prime * alpha * alpha * mu * mu * e_prime * std::sin(f_prime) / 2.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const third = 3.0 * H_prime * alpha * alpha * mu * mu * std::sin(2.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const fourth = 3.0 * H_prime * alpha * alpha * mu * mu * e_prime * std::sin(f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const fifth = H_prime * alpha * alpha * mu * mu * e_prime * std::sin(3.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const sixth = -H_prime * alpha * alpha * mu * mu * (1.0 - eta_prime) * (1.0 + 2.0*eta_prime) * std::sin(2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime / (1.0 + eta_prime);

    return first+second+third+fourth+fifth+sixth;
}

double L_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                    double f_prime, double alpha, double mu) {
    double const eta_prime = G_prime / L_prime;
    double const e_prime = std::sqrt(1.0 - G_prime * G_prime / L_prime / L_prime);
    double const s2 = 1.0 - H_prime * H_prime / G_prime / G_prime;

    double const dW1_dl = -alpha * alpha * mu * mu * (G_prime * G_prime - 3.0 * H_prime * H_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const df_dl = L_prime * L_prime * L_prime * (1.0 + e_prime * std::cos(f_prime)) * (1.0 + e_prime * std::cos(f_prime)) / G_prime / G_prime / G_prime;
    double const dfirst_df = alpha * alpha * mu * mu * (G_prime * G_prime - 3.0 * H_prime * H_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const dsecond_df = alpha * alpha * mu * mu * e_prime * (G_prime * G_prime - 3.0 * H_prime * H_prime) * std::cos(f_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const dthird_df = -3.0 * alpha * alpha * mu * mu * (G_prime * G_prime - H_prime * H_prime) * std::cos(2.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const dfourth_df = -3.0 * alpha * alpha * mu * mu * e_prime * (G_prime * G_prime - H_prime * H_prime) * std::cos(f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime / G_prime;
    double const dfifth_df = -3.0 * alpha * alpha * mu * mu * e_prime * (G_prime * G_prime - H_prime * H_prime) * std::cos(3.0*f_prime + 2.0*g_prime) / 8.0 / G_prime / G_prime / G_prime / G_prime / G_prime;

    double const dW1_df_df_dl = (dfirst_df+dsecond_df+dthird_df+dfourth_df+dfifth_df) * df_dl;

    return -(dW1_dl + dW1_df_df_dl);
}

double G_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                    double f_prime, double alpha, double mu) {
    double const eta_prime = G_prime / L_prime;
    double const e_prime = std::sqrt(1.0 - G_prime * G_prime / L_prime / L_prime);
    double const s2 = 1.0 - H_prime * H_prime / G_prime / G_prime;

    double const third = 3.0 * alpha * alpha * mu * mu * s2 * std::cos(2.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime;
    double const fourth = 3.0 * alpha * alpha * mu * mu * s2 * e_prime * std::cos(f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime;
    double const fifth = alpha * alpha * mu * mu * s2 * e_prime * std::cos(3.0*f_prime + 2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime;
    double const sixth = -alpha * alpha * mu * mu * s2 * (1.0 - eta_prime) * (1.0 + 2.0*eta_prime) * std::cos(2.0*g_prime) / 4.0 / G_prime / G_prime / G_prime / (1.0 + eta_prime);

    return third+fourth+fifth+sixth;
}

double H_prime_1(double l_prime, double g_prime, double h_prime,
                    double L_prime, double G_prime, double H_prime,
                    double f_prime, double alpha, double mu) {
    return 0;
}
