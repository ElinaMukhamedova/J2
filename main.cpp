#include "Kozai/transformations_YoshihideKozai.hpp"


int main() {
    double const mu = 398600441500000;
    double J2 = 0.001082636;
    double const R = 6371000;
    J2 *= R * R;

    double const l = -2.60400297084762;
    double const g = 0.2186689458739421;
    double const h = 1.10714871779409;
    double const L = 52088323885.51392;
    double const G = 50967166149.63148;
    double const H = 46103900839.90248;

    double const lp = -2.609002970847619;
    double const gp = 0.2186689458739421;
    double const hp = 1.10714871779409;
    double const Lp = 52062299240.42924;
    double const Gp = 50943964064.51405;
    double const Hp = 46092157963.13176;

    DelaunayElements const osculating{l, g, h, L, G, H};
    DelaunayElements const mean{lp, gp, hp, Lp, Gp, Hp};

    deltas(mean, osculating, mu, J2);

}