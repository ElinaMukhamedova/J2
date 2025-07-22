#pragma once

#include <cmath>

#include "celestial_mechanics/orbital_elements/AnomalyConverter.hpp"
#include "celestial_mechanics/orbital_elements/Elements.hpp"

double L_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
                double fp, double e, double mu, double J2) {
    double const eta = std::sqrt(1 - e * e);
    double const theta = Hp / Gp;
    double const gs = gp - (3 * mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * (1 - 5 * theta * theta) * (fp - lp);
    double const ap_rp = (1 + e * std::cos(fp)) / eta / eta;
    double const B20 = -(1 - 3 * theta * theta) / 4;
    double const B22 = 3 * (1 - theta * theta) / 4;

    double const delta_L1 = (mu * mu * J2 / Lp / Lp / Lp) * (B20 * (ap_rp * ap_rp * ap_rp - 1 / eta / eta / eta)
                                            + B22 * ap_rp * ap_rp * ap_rp * std::cos(2 * (fp + gs)));
    double const delta_L2 = (3 * mu * mu * mu * mu * J2 * J2 / 128 / Gp / Gp / Gp / Gp / Gp / Gp / Gp) *
                                        (8 * theta * theta * (1 - 5 * theta * theta)
                                        + e * e * (5 - 18 * theta * theta + 5 * theta * theta * theta * theta)
                                        - 2 * e * e * (1 - theta * theta) * (1 - 15 * theta * theta) * std::cos(2 * gp))
        + (3 * mu * mu * mu * mu * J2 * J2 / 512 / Gp / Gp / Gp / Gp / Gp / Gp / Lp) * (ap_rp * ap_rp * ap_rp) *
                                        (8 * (9 - 26 * theta * theta + 49 * theta * theta * theta * theta)
                                        + 4 * e * e * (37 - 98 * theta * theta + 37 * theta * theta * theta * theta)
                                        + 16 * eta * eta * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)
                                        + 2 * e * (16 * (1 - 3 * theta * theta) * (1 - 3 * theta * theta) *
                                                   (1 - eta * eta * eta) / e / e
                                                   + 4 * (19 - 30 * theta * theta + 35 * theta * theta * theta * theta)
                                                   + e * e * (73 - 234 * theta * theta
                                                              + 121 * theta * theta * theta * theta)) * std::cos(fp)
                                        + 4 * e * e * (4 * (1 - 3 * theta * theta) * (1 - 3 * theta * theta) *
                                                       (1 - eta * eta * eta) / e / e + 29 - 66 * theta * theta
                                                       + 45 * theta * theta * theta * theta) * std::cos(2*fp)
                                        + 2 * e * e * e * (11 - 30 * theta * theta
                                                           + 27 * theta * theta * theta * theta) * std::cos(3 * fp)
                                        + 4 * (1 - theta * theta) * (- 3 * e * e * e * (1 - 3 * theta * theta) *
                                                                     std::cos(fp - 2 * gp)
                                                                     - 2 * e * e * (1 - 3 * theta * theta) *
                                                                     ((1 - eta * eta * eta) / e / e + 8) * std::cos(2 * gp)
                                                                     + e * (4 * (1 - 3 * theta * theta) *
                                                                            (1 - eta * eta * eta) / e / e - 32 - e * e
                                                                            * (17 - 147 * theta * theta)) *
                                                                                                    std::cos(fp + 2 * gp)
                                                                     - 4 * (13 - 27 * theta * theta
                                                                            + 2 * e * e * (1 - 9 * theta * theta)
                                                                            + 3 * eta * eta * eta *
                                                                            (1 - 3 * theta * theta)) *
                                                                                                std::cos(2 * (fp + gp))
                                                                     - e * (28 * (1 - 3 * theta * theta) *
                                                                            (1 - eta * eta * eta) / e / e
                                                                            + 32 * (1 - 4 * theta * theta)
                                                                            - e * e * (15 - 77 * theta * theta)) *
                                                                                    std::cos(3 * fp + 2 * gp)
                                                                     - 2 * e * e * (1 - 3 * theta * theta) *
                                                                            (5 * (1 - eta * eta * eta) / e / e + 4) *
                                                                                    std::cos(4 * fp + 2 * gp)
                                                                     - 3 * e * e * e * (1 - 3 * theta * theta) *
                                                                                    std::cos(5 * fp + 2 * gp))
                                        + (1 - theta * theta) * (1 - theta * theta) * (9 * e * e * e *
                                                                                            std::cos(fp + 4 * gp)
                                                                                       + 54 * e * e *
                                                                                            std::cos(2 * fp + 4 * gp)
                                                                                       + e * (148 - 13 * e * e) *
                                                                                            std::cos(3 * fp + 4 * gp)
                                                                                       + 20 * (2 + 7 * e * e) *
                                                                                            std::cos(4 * fp + 4 * gp)
                                                                                       + 3 * e * (28 + 17 * e * e) *
                                                                                            std::cos(5 * fp + 4 * gp)
                                                                                       + 54 * e * e *
                                                                                            std::cos(6 * fp + 4 * gp)
                                                                                       + 9 * e * e * e *
                                                                                            std::cos(7 * fp + 4 * gp))
                                        );
    return Lp + delta_L1 + delta_L2;
}

double G_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
                double fp, double e, double mu, double J2) {
     double const eta = std::sqrt(1 - e * e);
     double const theta = Hp / Gp;
     double const gs = gp - (3 * mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * (1 - 5 * theta * theta) * (fp - lp);
     double const B22 = 3 * (1 - theta * theta) / 4;

     double const delta_G1 = mu * mu * J2 / 3 / Gp / Gp / Gp * B22 * (3 * e * std::cos(fp + 2 * gs)
                                                                      + 3 * std::cos(2 * fp + 2 * gs)
                                                                      + e * std::cos(3 * fp + 2 * gs));
     double const delta_G2 = mu * mu * mu * mu * J2 * J2 / 128 / Gp / Gp / Gp / Gp / Gp / Gp / Gp * (1 - theta * theta) *
          (-12 * e * e * (1 - 15 * theta * theta) * (fp - lp) * std::sin(2 * gp) - 4 * (7 - 25 * theta * theta)
               - 8 * e * e * (7 - 17 * theta * theta) - 96 * e * (1 - 3 * theta * theta) * std::cos(fp)
               - 24 * e * e * (1 - 3 * theta * theta) * cos(2 * fp) + 6 * e * (1 - 3 * theta * theta) * (1 - eta) * std::cos(fp - 2 * gp)
               + e * e * (20 * (1 - 3 * theta * theta) * (1 - eta) / e / e + (239 - 1581 * theta * theta) / 4 - 6 * eta * (1 - 15 * theta * theta)) * std::cos(2 * gp)
               + 6 * e * (4 * (1 - 3 * theta * theta) * (1 - eta) / e / e + 3 * (3 - 41 * theta * theta) - 5 * eta * (1 - 3 * theta * theta)) * std::cos(fp + 2 * gp)
               - 36 * (2 * (1 - theta * theta) + e * e * (1 - 3 * theta * theta)) * std::cos(2 * fp + 2 * gp)
               - 2 * e * (28 * (1 - 3 * theta * theta) * (1 - eta) / e / e + 43 - 161 * theta * theta + eta * (1 - 3 * theta * theta)) * std::cos(3 * fp + 2 * gp)
               - 3 * e * e * (12 * (1 - 3 * theta * theta) * (1 - eta) / e / e + 7 - 33 * theta * theta) * std::cos(4 * fp + 2 * gp)
               - 6 * e * (1 - 3 * theta * theta) * (1 - eta) * std::cos(5 * fp + 2 * gp)
               + 3 * (1 - theta * theta) * (5 * e * e * std::cos(2 * fp + 4 * gp)
                                             + 4 * e * std::cos(3 * fp + 4 * gp)
                                             - (4 - e * e) * std::cos(4 * fp + 4 * gp)
                                             - 4 * e * std::cos(5 * fp + 4 * gp)
                                             - e * e * std::cos(6 * fp + 4 * gp))
          );
     return Gp + delta_G1 + delta_G2;
}

double H_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
                double fp, double e, double mu, double J2) {
     return Hp;
}

double l_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
                double fp, double e, double mu, double J2) {
     double const eta = std::sqrt(1 - e * e);
     double const theta = Hp / Gp;
     double const gs = gp - (3 * mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * (1 - 5 * theta * theta) * (fp - lp);
     double const B20 = -(1 - 3 * theta * theta) / 4;
     double const B22 = 3 * (1 - theta * theta) / 4;

     double const delta_l1 = - mu * mu * J2 / 24 / Lp / Gp / Gp / Gp *
          (6 * B20 * (3 / e * (4 - e * e) * std::sin(fp) + 6 * std::sin(2 * fp) + e * std::sin(3 * fp))
               + B22 * (3 * e * std::sin(fp - 2 * gs) - 18 * std::sin(2 * gs) - 3 / e * (4 + 5 * e * e) * std::sin(fp + 2 * gs) + 1 / e * (28 - e * e) * std::sin(3 * fp + 2 * gs)
                         + 18 * std::sin(4 * fp + 2 * gs) + 3 * e * std::sin(5 * fp + 2 * gs)
                        )
          );
     double const delta_l2 = mu * mu * mu * mu * J2 * J2 * Lp / 4096 / Gp / Gp / Gp / Gp / Gp / Gp / Gp / Gp / Gp *
          (192 * eta * eta * eta * eta * (5 - 18 * theta * theta + 5 * theta * theta * theta * theta
                                             - 2 * (1 - theta * theta) * (1 - 15 * theta * theta) * std::cos(2 * gp)) * (fp - lp)
          + 48 / e * (2 * (1 - 3 * theta * theta) * (1 - 3 * theta * theta) * (8 - eta * eta * eta * (8 + 8 * e * e - 3 * e * e * e * e)) / e / e
                      + 4 * (29 - 74 * theta * theta - 11 * theta * theta * theta * theta)
                      + 2 * e * e * (7 - 22 * theta * theta + 155 * theta * theta * theta * theta)
                      - e * e * e * e * (47 - 166 * theta * theta + 167 * theta * theta * theta * theta)) * std::sin(fp)
          + 6 / e / e * (32 * (9 - 26 * theta * theta + 25 * theta * theta * theta * theta)
                         + 16 * e * e * (103 - 318 * theta * theta + 175 * theta * theta * theta * theta)
                         - 2 * e * e * e * e * (483 - 1614 * theta * theta + 563 * theta * theta * theta * theta)
                         - e * e * e * e * e * e * (13 - 82 * theta * theta + 125 * theta * theta * theta * theta)
                         + 128 * eta * eta * eta * eta * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)) * std::sin(2 * fp)
          + 8 / e * (8 * (47 - 102 * theta * theta + 63 * theta * theta * theta * theta)
                     + 40 * e * e * (10 - 39 * theta * theta + 30 * theta * theta * theta * theta)
                     - e * e * e * e * (325 - 1146 * theta * theta + 597 * theta * theta * theta * theta)
                     + 2 * eta * eta * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta) * (64 - 11 * e * e)) * std::sin(3 * fp)
          + 24 * (2 * (39 - 86 * theta * theta + 55 * theta * theta * theta * theta)
                  - 2 * e * e * (5 + 14 * theta * theta - 27 * theta * theta * theta * theta)
                  - e * e * e * e * (13 - 50 * theta * theta + 29 * theta * theta * theta * theta)
                  + 16 * eta * eta * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)) * std::sin(4 * fp)
          + 24 * e * (4 * (5 - 12 * theta * theta + 9 * theta * theta * theta * theta)
                      - 9 * e * e * (1 - theta * theta) * (1 - theta * theta)
                      + 2 * eta * eta * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)) * std::sin(5 * fp)
          + 2 * e * e * (11 - 30 * theta * theta + 27 * theta * theta * theta * theta) * (2 - e * e) * std::sin(6 * fp)
          + 4 * (1 - theta * theta) * (- 3 * e * e * (1 - 3 * theta * theta) * (2 - e * e) * std::sin(4 * fp - 2 * gp)
                                       - 6 * e * (1 - 3 * theta * theta) * (11 - 5 * e * e + eta * eta * eta) * std::sin(3 * fp - 2 * gp)
                                       - 12 * (2 * (1 - 3 * theta * theta) * (12 + eta * eta * eta)
                                               - 6 * e * e * (1 + theta * theta) - 3 * e * e * e * e * (1 - 11 * theta * theta)) * std::sin(2 * fp - 2 * gp)
                                       - 6 / e * ((1 - 3 * theta * theta) * (96 - 13 * e * e * eta * eta * eta) + e * e * (25 - 387 * theta * theta)
                                                 - 39 * e * e * e * e * (1 - 11 * theta * theta)) * std::sin(fp - 2 * gp)
                                       + 1 / e / e * (96 * (1 - 3 * theta * theta) * (4 - eta * eta * eta) + 2 * e * e * (113 - 1347 * theta * theta)
                                                      + 2 * e * e * e * e * (221 - 327 * theta * theta)
                                                      - e * e * e * e * e * e * (269 - 2151 * theta * theta)
                                                      - 16 * e * e * eta * eta * eta * (11 + 3 * theta * theta)
                                                      - 72 * e * e * e * e * eta * eta * eta * (1 - 15 * theta * theta)) * std::sin(2 * gp)
                                       + 12 / e * ((1 - 3 * theta * theta) * (8 - eta * eta * eta * (8 - 24 * e * e + 21 * e * e * e * e)) / e / e
                                                   - 48 * (1 - 4 * theta * theta) + e * e * (53 - 315 * theta * theta) + e * e * e * e * (41 - 15 * theta * theta)) * std::sin(fp + 2 * gp)
                                       - 48 * ((1 - 3 * theta * theta) * (16 - eta * eta * eta * (16 + 5 * e * e)) / e / e + 17 - 123 * theta * theta
                                               - 4 * e * e * (7 - 39 * theta * theta) - 5 * e * e * e * e * (1 - 3 * theta * theta)) * std::sin(2 * fp + 2 * gp)
                                       - 4 / e * ((1 - 3 * theta * theta) * (56 - eta * eta * eta * (56 + 152 * e * e + 11 * e * e * e * e))
                                                  + 16 * (12 - 59 * theta * theta) + e * e * (51 - 173 * theta * theta)
                                                  - e * e * e * e * (137 - 799 * theta * theta)) * std::sin(3 * fp + 2 * gp)
                                       - 3 / e / e * (8 * (1 - 3 * theta * theta) * (17 * e * e + 4 * eta * eta * eta * (7 - 3 * e * e)) + 2 * e * e * e * e * (77 - 279 * theta * theta)
                                                      - e * e * e * e * e * e * (29 - 183 * theta * theta)) * std::sin(4 * fp + 2 * gp)
                                       - 6 / e * (e * e * (37 - 39 * theta * theta) + 9 * e * e * e * e * (5 - 23 * theta * theta)
                                                  + eta * eta * eta * (1 - 3 * theta * theta) * (128 - 9 * e * e)) * std::sin(5 * fp + 2 * gp)
                                       - 4 * (2 * (1 - 3 * theta * theta) * (10 + 33 * eta * eta * eta) + 2 * e * e * (7 - 9 * theta * theta)
                                              + e * e * e * e * (11 - 57 * theta * theta)) * std::sin(6 * fp + 2 * gp)
                                       - 6 * e * (1 - 3 * theta * theta) * (7 - e * e + 5 * eta * eta * eta) * std::sin(7 * fp + 2 * gp)
                                       - 3 * e * e * (1 - 3 * theta * theta) * (2 - e * e) * std::sin(8 * fp + 2 * gp)
                                       )
          + (1 - theta * theta) * (1 - theta * theta) * (9 * e * e * (2 - e * e) * std::sin(2 * fp - 4 * gp)
                                                         + 108 * e * (2 - e * e) * std::sin(fp - 4 * gp)
                                                         - 12 * (66 - 34 * e * e + 13 * e * e * e * e) * std::sin(4 * gp)
                                                         - 12 / e * (72 + 2 * e * e + 49 * e * e * e * e) * std::sin(fp + 4 * gp)
                                                         - 3 / e / e * (96 + 368 * e * e + 158 * e * e * e * e + 161 * e * e * e * e * e * e) * std::sin(2 * fp + 4 * gp)
                                                         - 216 / e * (4 + 5 * e * e * e * e) * std::sin(3 * fp + 4 * gp)
                                                         + 144 * (17 - 16 * e * e - e * e * e * e) * std::sin(4 * fp + 4 * gp)
                                                         + 24 / e * (212 - 88 * e * e - 43 * e * e * e * e) * std::sin(5 * fp + 4 * gp)
                                                         + 1 / e / e * (1568 + 4112 * e * e - 3238 * e * e * e * e - 93 * e * e * e * e * e * e) * std::sin(6 * fp + 4 * gp)
                                                         + 12 / e * (168 + 50 * e * e - 95 * e * e * e * e) * std::sin(7 * fp + 4 * gp)
                                                         + 12 * (82 - 26 * e * e - 11 * e * e * e * e) * std::sin(8 * fp + 4 * gp)
                                                         + 108 * e * (2 - e * e) * std::sin(9 * fp + 4 * gp) + 9 * e * e * (2 - e * e) * std::sin(10 * fp + 4 * gp)
                                                         )
          );
     return lp + delta_l1 + delta_l2;
}

double g_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
                double fp, double e, double mu, double J2) {
     double const eta = std::sqrt(1 - e * e);
     double const theta = Hp / Gp;
     double const gs = gp - (3 * mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * (1 - 5 * theta * theta) * (fp - lp);
     double const B20 = -(1 - 3 * theta * theta) / 4;
     double const B22 = 3 * (1 - theta * theta) / 4;

     double const delta_g1 = mu * mu * J2 / 24 / Gp / Gp / Gp / Gp *
          (6 * B20 * (3 / e * (4 - e * e) * std::sin(fp) + 6 * std::sin(2 * fp) + e * std::sin(3 * fp))
           + B22 * (3 * e * std::sin(fp - 2 * gs) - 18 * std::sin(2 * gs)
                    - 3 / e * (4 + 5 * e * e) * std::sin(fp + 2 * gs)
                    + 1 / e * (28 - e * e) * std::sin(3 * fp + 2 * gs) + 18 * std::sin(4 * fp + 2 * gs) + 3 * e * std::sin(5 * fp + 2 * gs))
           - 18 * (1 - 5 * theta * theta) * (fp - lp + e * std::sin(fp))
           + 3 * (3 - 5 * theta * theta) * (3 * e * std::sin(fp + 2 * gs) + 3 * std::sin(2 * fp + 2 * gs) + e * std::sin(3 * fp + 2 * gs))
           );
     double const delta_g2 = - mu * mu * mu * mu * J2 * J2 / 2048 / Gp / Gp / Gp / Gp / Gp / Gp / Gp / Gp *
          (48 * (2 * (5 + 18 * theta * theta - 215 * theta * theta * theta * theta) + e * e * (25 - 126 * theta * theta + 45 * theta * theta * theta * theta)) * (fp - lp)
           - 96 * (2 * (1 - theta * theta) * (1 - 15  * theta * theta) + e * e * (5 - 112  * theta * theta + 135 * theta * theta * theta * theta)) * (fp - lp) * std::cos(2 * gp)
           + 48 / e * (8 * (1 - 3 * theta * theta) * (1 - 3 * theta * theta) * (1 - eta) / e / e
                       + 2 * (17 + 22 * theta * theta - 191 * theta * theta * theta * theta)
                       + 4 * e * e * (29 - 105 * theta * theta - 33 * theta * theta * theta * theta)
                       + eta * (1 - 3 * theta * theta) * (24 * (1 - 5 * theta * theta) - e * e * (5 - 27 * theta * theta))) * std::sin(fp)
           + 6 / e / e * (16 * (9 - 26 * theta * theta + 25 * theta * theta * theta * theta)
                          + 32 * e * e *(23 - 50 * theta * theta - 13 * theta * theta * theta * theta)
                          + e * e * e * e * (173 - 738 * theta * theta + 301 * theta * theta * theta * theta)
                          + 64 * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta + e * e * (1 - 6 * theta * theta))) * std::sin(2 * fp)
           + 8 / e * (4 * (47 - 102 * theta * theta + 63 * theta * theta * theta * theta)
                      + 6 * e * e * (41 - 112 * theta * theta + 37 * theta * theta * theta * theta)
                      + eta * (1 - 3 * theta * theta) *
                              (64 * (1 - 3 * theta * theta) + e * e * (5 - 39 * theta * theta))) * std::sin(3 * fp)
           + 12 * (2 * (39 - 86 * theta * theta + 55 * theta * theta * theta * theta)
                   + e * e * (21 - 66 * theta * theta + 37 * theta * theta * theta * theta)
                   + 16 * eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)) * std::sin(4 * fp)
           + 24 * e * (2 * (5 - 12 * theta * theta + 9 * theta * theta * theta * theta)
                       + eta * (1 - 3 * theta * theta) * (1 - 3 * theta * theta)) * std::sin(5 * fp)
           + 2 * e * e * (11 - 30 * theta * theta + 27 * theta * theta * theta * theta) * std::sin(6 * fp)
           - 12 * e * e * (1 - theta * theta) * (1 - 3 * theta * theta) * std::sin(4 * fp - 2 * gp)
           - 12 * e * (1 - theta * theta) * (1 - 3 * theta * theta) * (11 + eta) * std::sin(3 * fp - 2 * gp)
           - 24 * (1 - theta * theta) * (2 * (1 - 3 * theta * theta) * (12 + eta) + e * e * (5 - 39 * theta * theta)) * std::sin(2 * fp - 2 * gp)
           - 12 / e * (96 * (1 - theta * theta) * (1 - 3 * theta * theta) + e * e * (57 - 532 * theta * theta + 459 * theta * theta * theta * theta)
                       + e * e * eta * (1 - 3 * theta * theta) * (3 - 11 * theta * theta)) * std::sin(fp - 2 * gp)
           + 2 / e / e * (96 * (1 - theta * theta) * (1 - 3 * theta * theta) * (4 - eta)
                          + 2 * e * e * (169 - 1796 * theta * theta + 1467 * theta * theta * theta * theta)
                          - e * e * e * e * (709 - 6764 * theta * theta + 8739 * theta * theta * theta * theta)
                          - 16 * eta * (e * e * (3 + 10 * theta * theta - 33 * theta * theta * theta * theta)
                                        - 6 * e * e * e * e * (1 - 24 * theta * theta + 30 * theta * theta * theta * theta))) * std::sin(2 * gp)
           + 24 / e * (8 * (1 - theta * theta) * (1 - 3 * theta * theta) * (1 - eta) / e / e
                       - 8 * (9 - 44 * theta * theta + 39 * theta * theta * theta * theta)
                       - e * e * (139 - 2020 * theta * theta + 2337 * theta * theta * theta * theta)
                       + eta * (1 - 3 * theta * theta) * (8 * (7 - 9 * theta * theta) + e * e * (19 - 39 * theta * theta))) * std::sin(fp + 2 * gp)
           - 96 * ((1 - theta * theta) * (1 - 3 * theta * theta) * (16 - eta * (16 + 5 * e * e)) / e / e
                   + (19 - 256 * theta * theta + 237 * theta * theta * theta * theta)
                   - e * e * (1 + 44 * theta * theta - 33 * theta * theta * theta * theta)) * std::sin(2 * fp + 2 * gp)
           - 8 / e * (56 * (1 - theta * theta) * (1 - 3 * theta * theta) * (1 - eta) / e / e
                      + 8 * (3 - 44 * theta * theta + 13 * theta * theta * theta * theta)
                      - e * e * (69 + 324 * theta * theta + 79 * theta * theta * theta * theta)
                      + eta * (1 - 3 * theta * theta) * (8 * (9 - 23 * theta * theta) - e * e * (19 - 23 * theta * theta))) * std::sin(3 * fp + 2 * gp)
           + 12 / e / e * (12 * e * e * (7 - 32 * theta * theta + 33 * theta * theta * theta * theta)
                           + e * e * e * e * (21 - 76 * theta * theta + 159 * theta * theta * theta * theta)
                           - 16 * eta * (1 - 3 * theta * theta) * (7 * (1 - theta * theta) + 3 * e * e * (1 - 2 * theta * theta))) * std::sin(4 * fp + 2 * gp)
           + 12 / e * (e * e * (43 - 252 * theta * theta + 225 * theta * theta * theta * theta)
                       - eta * (1 - 3 * theta * theta) * (128 * (1 - theta * theta) + e * e * (7 - 15 * theta * theta))) * std::sin(5 * fp + 2 * gp)
           - 8 * (1 - theta * theta) * (2 * (1 - 3 * theta * theta) * (10 + 33 * eta) - 3 * e * e * (3 - 17 * theta * theta)) * std::sin(6 * fp + 2 * gp)
           - 12 * e * (1 - theta * theta) * (1 - 3 * theta * theta) * (7 + 5 * eta) * std::sin(7 * fp + 2 * gp)
           - 12 * e * e * (1 - theta * theta) * (1 - 3 * theta * theta) * std::sin(8 * fp + 2 * gp)
           + 9 * e * e * (1 - theta * theta) * (1 - theta * theta) * std::sin(2 * fp - 4 * gp)
           + 108 * e * (1 - theta * theta) * (1 - theta * theta) * std::sin(fp - 4 * gp)
           - 6 * (1 - theta * theta) * (1 - theta * theta) * (66 - e * e) * std::sin(4 * gp)
           - 12 / e * (1 - theta * theta) * (1 - theta * theta) * (36 + 19 * e * e) * std::sin(fp + 4 * gp)
           - 3 / e / e * (16 * (1 - theta * theta) * (1 - theta * theta) * (3 + 13 * e * e) + e * e * e * e * (7 + 194 * theta * theta - 9 * theta * theta * theta * theta)) * std::sin(2 * fp + 4 * gp)
           - 24 / e * (18 * (1 - theta * theta) * (1 - theta * theta) - e * e * (19 - 102 * theta * theta + 35 * theta * theta * theta * theta)) * std::sin(3 * fp + 4 * gp)
           + 12 * (2 * (83 - 210 * theta * theta + 103 * theta * theta * theta * theta) + e * e * (11 - 66 * theta * theta + 23 * theta * theta * theta * theta)) * std::sin(4 * fp + 4 * gp)
           + 24 / e * (106 * (1 - theta * theta) * (1 - theta * theta) + e * e * (53 - 138 * theta * theta + 69 * theta * theta * theta * theta)) * std::sin(5 * fp + 4 * gp)
           + 1 / e / e * (16 * (1 - theta * theta) * (1 - theta * theta) * (49 + 153 * e * e)
                          + e * e * e * e * (197 - 538 * theta * theta + 277 * theta * theta * theta * theta)) * std::sin(6 * fp + 4 * gp)
           + 12 / e * (1 - theta * theta) * (1 - theta * theta) * (84 + 67 * e * e) * std::sin(7 * fp + 4 * gp)
           + 6 * (1 - theta * theta) * (1 - theta * theta) * (82 + 15 * e * e) * std::sin(8 * fp + 4 * gp)
           + 108 * e * (1 - theta * theta) * (1 - theta * theta) * std::sin(9 * fp + 4 * gp)
           + 9 * e * e * (1 - theta * theta) * (1 - theta * theta) * std::sin(10 * fp + 4 * gp)
          );
    return gp + delta_g1 + delta_g2;
}

double h_osc(double Lp, double Gp, double Hp, double lp, double gp, double hp,
             double fp, double e, double mu, double J2) {
    double const eta = std::sqrt(1 - e * e);
    double const theta = Hp / Gp;
    double const gs = gp - (3 * mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * (1 - 5 * theta * theta) * (fp - lp);
    double const B20 = -(1 - 3 * theta * theta) / 4;
    double const B22 = 3 * (1 - theta * theta) / 4;

    double const delta_h1 = (mu * mu * J2 / 4 / Gp / Gp / Gp / Gp) * theta *
                (-6 * (fp - lp + e * std::sin(fp)) + 3 * e * std::sin(fp + 2 * gs)
                 + 3 * std::sin(2 * fp + 2 * gs) + e * std::sin(3 * fp + 2 * gs)
                );
    double const delta_h2 = mu * mu * mu * mu * J2 * J2 / 128 / Gp / Gp / Gp / Gp / Gp / Gp / Gp / Gp * theta *
        (12 * (4 * (1 - 10 * theta * theta) - e * e * (9 - 5 * theta * theta) + 2 * e * e * (8 - 15 * theta * theta)) * (fp - lp)
         + 12 * e * (3 * (1 - 3 * theta * theta) * (4 - eta * (4 - e * e)) / e / e - (17 + 21 * theta * theta)) * std::sin(fp)
         + 6 * (12 * (1 - 3 * theta * theta) * (1 - eta) / e / e - (9 - 5 * theta * theta)) * std::sin(2 * fp)
         + 12 * e * (1 - 3 * theta * theta) * (1 - eta) * std::sin(3 * fp)
         - 6 * e * (1 - 3 * theta * theta) * (1 - eta) * std::sin(fp - 2 * gp)
         - e * e * (4 * (7 + 3 * theta * theta) * (1 - eta) / e / e - 191 + 1053 * theta * theta / 2 + 12 * eta * (8 - 15 * theta * theta)) * std::sin(2 * gp)
         + 6 * e * ((1 - 3 * theta * theta) * (4 - eta * (4 + 5 * e * e)) / e / e + 3 * (35 - 73 * theta * theta)) * std::sin(fp + 2 * gp)
         + 24 * (6 * (1 - theta * theta) + e * e * (4 - theta * theta)) * std::sin(2 * fp + 2 * gp)
         - 2 * e * ((1 - 3 * theta * theta) * (28 - eta * (28 - e * e)) / e / e - (21 + 97 * theta * theta)) * std::sin(3 * fp + 2 * gp)
         - 6 * e * e * (6 * (1 - 3 * theta * theta) * (1 - eta) / e / e - 13 * e * e * theta * theta) * std::sin(4 * fp + 2 * gp)
         - 6 * e * (1 - 3 * theta * theta) * (1 - eta) * std::sin(5 * fp + 2 * gp)
         - 3 * e * e * (7 + 5 * theta * theta) * std::sin(2 * fp + 4 * gp)
         - 12 * e * (5 + theta * theta) * std::sin(3 * fp + 4 * gp)
         - 3 * (4 * (4 - theta * theta) + e * e * (7 + theta * theta)) * std::sin(4 * fp + 4 * gp)
         - 12 * e * (3 - theta * theta) * std::sin(5 * fp + 4 * gp)
         - e * e * (7 - 3 * theta * theta) * std::sin(6 * fp + 4 * gp)
        );
    return hp + delta_h1 + delta_h2;
}

DelaunayElements MeanToOsculating(DelaunayElements const& mean, double mu, double J2) {
    double const e = std::sqrt(1 - mean.G * mean.G / mean.L / mean.L);
    double const meanAnomaly = mean.l;
    double const eccentricAnomaly = Orbit::MeanToEccentric(e, meanAnomaly);
    double const trueAnomaly = Orbit::EccentricToTrue(e, eccentricAnomaly);
    double const l = l_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    double const g = g_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    double const h = h_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    double const L = L_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    double const G = G_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    double const H = H_osc(mean.L, mean.G, mean.H, mean.l, mean.g, mean.h, trueAnomaly, e, mu, J2);
    return DelaunayElements{l, g, h, L, G, H};
}

DelaunayElements OsculatingToMean(DelaunayElements const& osculating, double mu, double J2) {
    double const e = std::sqrt(1 - osculating.G * osculating.G / osculating.L / osculating.L);
    double const meanAnomaly = osculating.l;
    double const eccentricAnomaly = Orbit::MeanToEccentric(e, meanAnomaly);
    double const trueAnomaly = Orbit::EccentricToTrue(e, eccentricAnomaly);

    double const l1 = l_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);
    double const g1 = g_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);
    double const h1 = h_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);
    double const L1 = L_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);
    double const G1 = G_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);
    double const H1 = H_osc(osculating.L, osculating.G, osculating.H, osculating.l, osculating.g, osculating.h, trueAnomaly, e, mu, J2);

    double const l2 = l_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);
    double const g2 = g_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);
    double const h2 = h_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);
    double const L2 = L_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);
    double const G2 = G_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);
    double const H2 = H_osc(L1, G1, H1, l1, g1, h1, trueAnomaly, e, mu, J2);

    double const l3 = l_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);
    double const g3 = g_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);
    double const h3 = h_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);
    double const L3 = L_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);
    double const G3 = G_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);
    double const H3 = H_osc(L2, G2, H2, l2, g2, h2, trueAnomaly, e, mu, J2);

    return DelaunayElements{l3, g3, h3, L3, G3, H3};
}