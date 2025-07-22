#include "transformations.hpp"
#include "osculatingDelaunay.hpp"
#include "second_order.hpp"


DelaunayElements MeanToOsculating(DelaunayElements const& mean,double const f_prime, double const alpha, double const mu) {
    double const l_prime = mean.l;
    double const g_prime = mean.g;
    double const h_prime = mean.h;
    double const L_prime = mean.L;
    double const G_prime = mean.G;
    double const H_prime = mean.H;

    double const J2 = 0.001082636;

    double const l = l_prime
                    + J2 * l_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * l_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);
    double const g = g_prime
                    + J2 * g_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * g_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);
    double const h = h_prime
                    + J2 * h_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * h_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);
    double const L = L_prime
                    + J2 * L_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * L_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);
    double const G = G_prime
                    + J2 * G_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * G_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);
    double const H = H_prime
                    + J2 * H_prime_1(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu)
                    + 0.5 * J2 * J2 * H_2(l_prime, g_prime, h_prime, L_prime, G_prime, H_prime, f_prime, alpha, mu);

    return DelaunayElements{l, g, h, L, G, H};
}

DelaunayElements OsculatingToMean(DelaunayElements const& osc, double const f, double const alpha, double const mu) {
    double const l = osc.l;
    double const g = osc.g;
    double const h = osc.h;
    double const L = osc.L;
    double const G = osc.G;
    double const H = osc.H;

    double const J2 = 0.001082636;

    double const l_prime = l
                        - J2 * l_prime_1(l, g, h, L, G, H, f, alpha, mu);
    double const g_prime = g
                        - J2 * g_prime_1(l, g, h, L, G, H, f, alpha, mu);
    double const h_prime = h
                        - J2 * h_prime_1(l, g, h, L, G, H, f, alpha, mu);
    double const L_prime = L
                        - J2 * L_prime_1(l, g, h, L, G, H, f, alpha, mu);
    double const G_prime = G
                        - J2 * G_prime_1(l, g, h, L, G, H, f, alpha, mu);
    double const H_prime = H
                        - J2 * H_prime_1(l, g, h, L, G, H, f, alpha, mu);

    return DelaunayElements{l_prime, g_prime, h_prime, L_prime, G_prime, H_prime};
}