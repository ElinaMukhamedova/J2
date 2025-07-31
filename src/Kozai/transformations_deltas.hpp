#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <array>
#include "celestial_mechanics/orbital_elements/AnomalyConverter.hpp"
#include "celestial_mechanics/orbital_elements/Elements.hpp"


inline double wrapToPi(double angle) {
    angle = std::fmod(angle + std::numbers::pi, 2 * std::numbers::pi);
    if (angle < 0)
        angle += 2 * std::numbers::pi;
    return angle - std::numbers::pi;
}

std::array<double, 6> deltas(DelaunayElements const& mean, DelaunayElements const& osculating, double mu, double J2) {
    double const L = osculating.L;
    double const G = osculating.G;
    double const H = osculating.H;
    double const l = osculating.l;
    double const g = osculating.g;
    double const h = osculating.h;

    double const Lp = mean.L;
    double const Gp = mean.G;
    double const Hp = mean.H;
    double const lp = mean.l;
    double const gp = mean.g;
    double const hp = mean.h;

    double const e2 = 1 - G * G / L / L;
    double const e = std::sqrt(e2);
    double const e3 = e * e2;
    double const e4 = e2 * e2;
    double const e6 = e2 * e2 * e2;

    double const ep2 = 1 - Gp * Gp / Lp / Lp;
    double const ep = std::sqrt(ep2);
    double const ep3 = ep * ep2;
    double const ep4 = ep2 * ep2;
    double const ep6 = ep2 * ep2 * ep2;

    double const eta2 = 1 - e2;
    double const eta = std::sqrt(eta2);
    double const eta3 = eta2 * eta;
    double const eta4 = eta2 * eta2;
    double const eta5 = eta3 * eta2;

    double const etap2 = 1 - ep2;
    double const etap = std::sqrt(etap2);
    double const etap3 = etap2 * etap;
    double const etap4 = etap2 * etap2;
    double const etap5 = etap3 * etap2;

    double const theta = H / G;
    double const theta2 = theta * theta;
    double const theta4 = theta * theta * theta * theta;

    double const thetap = Hp / Gp;
    double const thetap2 = thetap * thetap;
    double const thetap4 = thetap * thetap * thetap * thetap;

    double const mu2 = mu * mu;
    double const mu4 = mu * mu * mu * mu;
    double const J22 = J2 * J2;
    double const L3 = L * L * L;
    double const Lp3 = Lp * Lp * Lp;
    double const G4 = G * G * G * G;
    double const Gp4 = Gp * Gp * Gp * Gp;
    double const G3 = G * G * G;
    double const Gp3 = Gp * Gp * Gp;
    double const G6 = G3 * G3;
    double const Gp6 = Gp3 * Gp3;
    double const G7 = G3 * G4;
    double const Gp7 = Gp3 * Gp4;
    double const G8 = G * G7;
    double const Gp8 = Gp * Gp7;
    double const G9 = G3 * G6;
    double const Gp9 = Gp3 * Gp6;

    double const E = Orbit::MeanToEccentric(e, l);
    double const f = wrapToPi(Orbit::EccentricToTrue(e, E));

    double const Ep = Orbit::MeanToEccentric(ep, lp);
    double const fp = wrapToPi(Orbit::EccentricToTrue(ep, Ep));

    double const gs = gp - (3 * mu2 * J2 / 4 / G4) * (1 - 5 * theta2) * (f - l);

    double const cos_f = std::cos(f);
    double const cos_2f = std::cos(2 * f);
    double const cos_3f = std::cos(3 * f);
    double const cos_f_2gs = std::cos(f + 2 * gs);
    double const cos_2f_2gs = std::cos(2 * (f + gs));
    double const cos_3f_2gs = std::cos(3 * f + 2 * gs);
    double const cos_2gp = std::cos(2 * gp);
    double const cos_f_2g = std::cos(f + 2 * g);
    double const cos_f_2gp = std::cos(f + 2 * gp);
    double const cos_2f_2g = std::cos(2 * (f + g));
    double const cos_2f_2gp = std::cos(2 * (f + gp));
    double const cos_3f_2g = std::cos(3 * f + 2 * g);
    double const cos_3f_2gp = std::cos(3 * f + 2 * gp);
    double const cos_4f_2g = std::cos(4 * f + 2 * g);
    double const cos_4f_2gp = std::cos(4 * f + 2 * gp);
    double const cos_5f_2g = std::cos(5 * f + 2 * g);
    double const cos_5f_2gp = std::cos(5 * f + 2 * gp);
    double const cos_f_m_2g = std::cos(f - 2 * g);
    double const cos_f_m_2gp = std::cos(f - 2 * gp);
    double const cos_f_4gp = std::cos(f + 4 * gp);
    double const cos_2f_4gp = std::cos(2 * f + 4 * gp);
    double const cos_3f_4gp = std::cos(3 * f + 4 * gp);
    double const cos_4f_4gp = std::cos(4 * f + 4 * gp);
    double const cos_5f_4gp = std::cos(5 * f + 4 * gp);
    double const cos_6f_4gp = std::cos(6 * f + 4 * gp);
     double const cos_7f_4gp = std::cos(7 * f + 4 * gp);
     double const cos_f_4g = std::cos(f + 4 * g);
     double const cos_2f_4g = std::cos(2 * f + 4 * g);
     double const cos_3f_4g = std::cos(3 * f + 4 * g);
     double const cos_4f_4g = std::cos(4 * f + 4 * g);
     double const cos_5f_4g = std::cos(5 * f + 4 * g);
     double const cos_6f_4g = std::cos(6 * f + 4 * g);
     double const cos_7f_4g = std::cos(7 * f + 4 * g);

    double const cos_fp = std::cos(fp);
    double const cos_2fp = std::cos(2 * fp);
    double const cos_3fp = std::cos(3 * fp);
    double const cos_fp_2gs = std::cos(fp + 2 * gs);
    double const cos_2fp_2gs = std::cos(2 * (fp + gs));
    double const cos_3fp_2gs = std::cos(3 * fp + 2 * gs);
    double const cos_2g = std::cos(2 * g);

    double const sin_f = std::sin(f);
    double const sin_2f = std::sin(2 * f);
    double const sin_3f = std::sin(3 * f);
    double const sin_4f = std::sin(4 * f);
    double const sin_5f = std::sin(5 * f);
    double const sin_6f = std::sin(6 * f);
    double const sin_2gs = std::sin(2 * gs);
    double const sin_f_2gs = std::sin(f + 2 * gs);
    double const sin_f_m_2gs = std::sin(f - 2 * gs);
    double const sin_2f_2gs = std::sin(2 * f + 2 * gs);
    double const sin_3f_2gs = std::sin(3 * f + 2 * gs);
    double const sin_4f_2gs = std::sin(4 * f + 2 * gs);
    double const sin_5f_2gs = std::sin(5 * f + 2 * gs);
    double const sin_2g = std::sin(2 * g);
    double const sin_2gp = std::sin(2 * gp);
    double const sin_f_2gp = std::sin(f + 2 * gp);
    double const sin_2f_2gp = std::sin(2 * f + 2 * gp);
    double const sin_3f_2gp = std::sin(3 * f + 2 * gp);
    double const sin_4f_2gp = std::sin(4 * f + 2 * gp);
    double const sin_5f_2gp = std::sin(5 * f + 2 * gp);
    double const sin_6f_2gp = std::sin(6 * f + 2 * gp);
    double const sin_7f_2gp = std::sin(7 * f + 2 * gp);
    double const sin_8f_2gp = std::sin(8 * f + 2 * gp);
    double const sin_f_m_2gp = std::sin(f - 2 * gp);
    double const sin_2f_m_2gp = std::sin(2 * f - 2 * gp);
    double const sin_3f_m_2gp = std::sin(3 * f - 2 * gp);
    double const sin_4f_m_2gp = std::sin(4 * f - 2 * gp);
    double const sin_f_m_4gp = std::sin(f - 4 * gp);
    double const sin_2f_m_4gp = std::sin(2 * f - 4 * gp);
    double const sin_4gp = std::sin(4 * gp);
    double const sin_f_4gp = std::sin(f + 4 * gp);
    double const sin_2f_4gp = std::sin(2 * f + 4 * gp);
    double const sin_3f_4gp = std::sin(3 * f + 4 * gp);
    double const sin_4f_4gp = std::sin(4 * f + 4 * gp);
    double const sin_5f_4gp = std::sin(5 * f + 4 * gp);
    double const sin_6f_4gp = std::sin(6 * f + 4 * gp);
    double const sin_7f_4gp = std::sin(7 * f + 4 * gp);
    double const sin_8f_4gp = std::sin(8 * f + 4 * gp);
    double const sin_9f_4gp = std::sin(9 * f + 4 * gp);
    double const sin_10f_4gp = std::sin(10 * f + 4 * gp);

    double const sin_fp = std::sin(fp);
    double const sin_2fp = std::sin(2 * fp);
    double const sin_3fp = std::sin(3 * fp);
    double const sin_4fp = std::sin(4 * fp);
    double const sin_5fp = std::sin(5 * fp);
    double const sin_6fp = std::sin(6 * fp);
    double const sin_fp_2gs = std::sin(fp + 2 * gs);
    double const sin_fp_m_2gs = std::sin(fp - 2 * gs);
    double const sin_2fp_2gs = std::sin(2 * fp + 2 * gs);
    double const sin_3fp_2gs = std::sin(3 * fp + 2 * gs);
    double const sin_4fp_2gs = std::sin(4 * fp + 2 * gs);
    double const sin_5fp_2gs = std::sin(5 * fp + 2 * gs);
    double const sin_f_2g = std::sin(f + 2 * g);
    double const sin_2f_2g = std::sin(2 * f + 2 * g);
    double const sin_3f_2g = std::sin(3 * f + 2 * g);
    double const sin_4f_2g = std::sin(4 * f + 2 * g);
    double const sin_5f_2g = std::sin(5 * f + 2 * g);
    double const sin_6f_2g = std::sin(6 * f + 2 * g);
    double const sin_7f_2g = std::sin(7 * f + 2 * g);
    double const sin_8f_2g = std::sin(8 * f + 2 * g);
    double const sin_f_m_2g = std::sin(f - 2 * g);
    double const sin_2f_m_2g = std::sin(2 * f - 2 * g);
    double const sin_3f_m_2g = std::sin(3 * f - 2 * g);
    double const sin_4f_m_2g = std::sin(4 * f - 2 * g);
    double const sin_f_m_4g = std::sin(f - 4 * g);
    double const sin_2f_m_4g = std::sin(2 * f - 4 * g);
    double const sin_4g = std::sin(4 * g);
    double const sin_f_4g = std::sin(f + 4 * g);
    double const sin_2f_4g = std::sin(2 * f + 4 * g);
    double const sin_3f_4g = std::sin(3 * f + 4 * g);
    double const sin_4f_4g = std::sin(4 * f + 4 * g);
    double const sin_5f_4g = std::sin(5 * f + 4 * g);
    double const sin_6f_4g = std::sin(6 * f + 4 * g);
    double const sin_7f_4g = std::sin(7 * f + 4 * g);
    double const sin_8f_4g = std::sin(8 * f + 4 * g);
    double const sin_9f_4g = std::sin(9 * f + 4 * g);
    double const sin_10f_4g = std::sin(10 * f + 4 * g);

    double const ap_rp = (1 + ep * cos_fp) / etap2;
    double const ap_rp3 = ap_rp * ap_rp * ap_rp;
    double const a_r = (1 + ep * cos_f) / etap2;
    double const a_r3 = a_r * a_r * a_r;
    double const B20 = -(1 - 3 * thetap2) / 4;
    double const B22 = 3 * (1 - thetap2) / 4;

    double const delta_L1 = (mu2 * J2 / Lp3) * (B20 * (ap_rp3 - 1 / etap3) + B22 * ap_rp3 * cos_2fp_2gs);
    double const delta_L2 = (3 * mu4 * J22 / 128 / G7) *
                                        (8 * theta2 * (1 - 5 * theta2)
                                        + e2 * (5 - 18 * theta2 + 5 * theta4)
                                        - 2 * e2 * (1 - theta2) * (1 - 15 * theta2) * cos_2g)
        + (3 * mu4 * J22 / 512 / G6 / L) * (a_r3) *
                                        (8 * (9 - 26 * theta2 + 49 * theta4)
                                        + 4 * e2 * (37 - 98 * theta2 + 37 * theta4)
                                        + 16 * eta3 * (1 - 3 * theta2) * (1 - 3 * theta2)
                                        + 2 * e * (16 * (1 - 3 * theta2) * (1 - 3 * theta2) *
                                                   (1 - eta3) / e2
                                                   + 4 * (19 - 30 * theta2 + 35 * theta4)
                                                   + e2 * (73 - 234 * theta2
                                                              + 121 * theta4)) * cos_f
                                        + 4 * e2 * (4 * (1 - 3 * theta2) * (1 - 3 * theta2) *
                                                       (1 - eta3) / e2 + 29 - 66 * theta2
                                                       + 45 * theta4) * cos_2f
                                        + 2 * e3 * (11 - 30 * theta2
                                                           + 27 * theta4) * cos_3f
                                        + 4 * (1 - theta2) * (- 3 * e3 * (1 - 3 * theta2) *
                                                                     cos_f_m_2g
                                                                     - 2 * e2 * (1 - 3 * theta2) *
                                                                     ((1 - eta3) / e2 + 8) * cos_2g
                                                                     + e * (4 * (1 - 3 * theta2) *
                                                                            (1 - eta3) / e2 - 32 - e2
                                                                            * (17 - 147 * theta2)) *
                                                                                                    cos_f_2g
                                                                     - 4 * (13 - 27 * theta2
                                                                            + 2 * e2 * (1 - 9 * theta2)
                                                                            + 3 * eta3 *
                                                                            (1 - 3 * theta2)) *
                                                                                                cos_2f_2g
                                                                     - e * (28 * (1 - 3 * theta2) *
                                                                            (1 - eta3) / e2
                                                                            + 32 * (1 - 4 * theta2)
                                                                            - e2 * (15 - 77 * theta2)) *
                                                                                    cos_3f_2g
                                                                     - 2 * e2 * (1 - 3 * theta2) *
                                                                            (5 * (1 - eta3) / e2 + 4) *
                                                                                    cos_4f_2g
                                                                     - 3 * e3 * (1 - 3 * theta2) *
                                                                                    cos_5f_2g)
                                        + (1 - theta2) * (1 - theta2) * (9 * e3 *
                                                                                            cos_f_4g
                                                                                       + 54 * e2 *
                                                                                            cos_2f_4g
                                                                                       + e * (148 - 13 * e2) *
                                                                                            cos_3f_4g
                                                                                       + 20 * (2 + 7 * e2) *
                                                                                            cos_4f_4g
                                                                                       + 3 * e * (28 + 17 * e2) *
                                                                                            cos_5f_4g
                                                                                       + 54 * e2 *
                                                                                            cos_6f_4g
                                                                                       + 9 * e3 *
                                                                                            cos_7f_4g));


    double const delta_G1 = mu2 * J2 / 3 / Gp3 * B22 * (3 * ep * cos_fp_2gs + 3 * cos_2fp_2gs + ep * cos_3fp_2gs);
    double const delta_G2 = mu4 * J22 / 128 / G7 * (1 - theta2) *
           (-12 * e2 * (1 - 15 * theta2) * (f - l) * sin_2g - 4 * (7 - 25 * theta2)
                - 8 * e2 * (7 - 17 * theta2) - 96 * e * (1 - 3 * theta2) * cos_f
                - 24 * e2 * (1 - 3 * theta2) * cos_2f + 6 * e * (1 - 3 * theta2) * (1 - eta) * cos_f_m_2g
                + e2 * (20 * (1 - 3 * theta2) * (1 - eta) / e2 + (239 - 1581 * theta2) / 4 - 6 * eta * (1 - 15 * theta2)) * cos_2g
                + 6 * e * (4 * (1 - 3 * theta2) * (1 - eta) / e2 + 3 * (3 - 41 * theta2) - 5 * eta * (1 - 3 * theta2)) * cos_f_2g
                - 36 * (2 * (1 - theta2) + e2 * (1 - 3 * theta2)) * cos_2f_2g
                - 2 * e * (28 * (1 - 3 * theta2) * (1 - eta) / e2 + 43 - 161 * theta2 + eta * (1 - 3 * theta2)) * cos_3f_2g
                - 3 * e2 * (12 * (1 - 3 * theta2) * (1 - eta) / e2 + 7 - 33 * theta2) * cos_4f_2g
                - 6 * e * (1 - 3 * theta2) * (1 - eta) * cos_5f_2g
                + 3 * (1 - theta2) * (5 * e2 * cos_2f_4g
                                              + 4 * e * cos_3f_4g
                                              - (4 - e2) * cos_4f_4g
                                              - 4 * e * cos_5f_4g
                                              - e2 * cos_6f_4g)
           );

    double const delta_l1 = - mu2 * J2 / 24 / Lp / Gp3 *
           (6 * B20 * (3 / ep * (4 - ep2) * sin_fp + 6 * sin_2fp + ep * sin_3fp)
                + B22 * (3 * ep * sin_fp_m_2gs - 18 * sin_2gs - 3 / ep * (4 + 5 * ep2) * sin_fp_2gs + 1 / ep * (28 - ep2) * sin_3fp_2gs
                          + 18 * sin_4fp_2gs + 3 * ep * sin_5fp_2gs
                         )
           );
    double const delta_l2 = mu4 * J22 * L / 4096 / G9 *
          (192 * eta4 * (5 - 18 * theta2 + 5 * theta4
                                             - 2 * (1 - theta2) * (1 - 15 * theta2) * cos_2g) * (f - l)
          + 48 / e * (2 * (1 - 3 * theta2) * (1 - 3 * theta2) * (8 - eta3 * (8 + 8 * e2 - 3 * e4)) / e2
                      + 4 * (29 - 74 * theta2 - 11 * theta4)
                      + 2 * e2 * (7 - 22 * theta2 + 155 * theta4)
                      - e4 * (47 - 166 * theta2 + 167 * theta4)) * sin_f
          + 6 / e2 * (32 * (9 - 26 * theta2 + 25 * theta4)
                         + 16 * e2 * (103 - 318 * theta2 + 175 * theta4)
                         - 2 * e4 * (483 - 1614 * theta2 + 563 * theta4)
                         - e6 * (13 - 82 * theta2 + 125 * theta4)
                         + 128 * eta5 * (1 - 3 * theta2) * (1 - 3 * theta2)) * sin_2f
          + 8 / e * (8 * (47 - 102 * theta2 + 63 * theta4)
                     + 40 * e2 * (10 - 39 * theta2 + 30 * theta4)
                     - e4 * (325 - 1146 * theta2 + 597 * theta4)
                     + 2 * eta3 * (1 - 3 * theta2) * (1 - 3 * theta2) * (64 - 11 * e2)) * sin_3f
          + 24 * (2 * (39 - 86 * theta2 + 55 * theta4)
                  - 2 * e2 * (5 + 14 * theta2 - 27 * theta4)
                  - e4 * (13 - 50 * theta2 + 29 * theta4)
                  + 16 * eta3 * (1 - 3 * theta2) * (1 - 3 * theta2)) * sin_4f
          + 24 * e * (4 * (5 - 12 * theta2 + 9 * theta4)
                      - 9 * e2 * (1 - theta2) * (1 - theta2)
                      + 2 * eta3 * (1 - 3 * theta2) * (1 - 3 * theta2)) * sin_5f
          + 2 * e2 * (11 - 30 * theta2 + 27 * theta4) * (2 - e2) * sin_6f
          + 4 * (1 - theta2) * (- 3 * e2 * (1 - 3 * theta2) * (2 - e2) * sin_4f_m_2g
                                       - 6 * e * (1 - 3 * theta2) * (11 - 5 * e2 + eta3) * sin_3f_m_2g
                                       - 12 * (2 * (1 - 3 * theta2) * (12 + eta3)
                                               - 6 * e2 * (1 + theta2) - 3 * e4 * (1 - 11 * theta2)) * sin_2f_m_2g
                                       - 6 / e * ((1 - 3 * theta2) * (96 - 13 * e2 * eta3) + e2 * (25 - 387 * theta2)
                                                 - 39 * e4 * (1 - 11 * theta2)) * sin_f_m_2g
                                       + 1 / e2 * (96 * (1 - 3 * theta2) * (4 - eta3) + 2 * e2 * (113 - 1347 * theta2)
                                                      + 2 * e4 * (221 - 327 * theta2)
                                                      - e6 * (269 - 2151 * theta2)
                                                      - 16 * e2 * eta3 * (11 + 3 * theta2)
                                                      - 72 * e4 * eta3 * (1 - 15 * theta2)) * sin_2g
                                       + 12 / e * ((1 - 3 * theta2) * (8 - eta3 * (8 - 24 * e2 + 21 * e4)) / e2
                                                   - 48 * (1 - 4 * theta2) + e2 * (53 - 315 * theta2) + e4 * (41 - 15 * theta2)) * sin_f_2g
                                       - 48 * ((1 - 3 * theta2) * (16 - eta3 * (16 + 5 * e2)) / e2 + 17 - 123 * theta2
                                               - 4 * e2 * (7 - 39 * theta2) - 5 * e4 * (1 - 3 * theta2)) * sin_2f_2g
                                       - 4 / e * ((1 - 3 * theta2) * (56 - eta3 * (56 + 152 * e2 + 11 * e4))
                                                  + 16 * (12 - 59 * theta2) + e2 * (51 - 173 * theta2)
                                                  - e4 * (137 - 799 * theta2)) * sin_3f_2g
                                       - 3 / e2 * (8 * (1 - 3 * theta2) * (17 * e2 + 4 * eta3 * (7 - 3 * e2)) + 2 * e4 * (77 - 279 * theta2)
                                                      - e6 * (29 - 183 * theta2)) * sin_4f_2g
                                       - 6 / e * (e2 * (37 - 39 * theta2) + 9 * e4 * (5 - 23 * theta2)
                                                  + eta3 * (1 - 3 * theta2) * (128 - 9 * e2)) * sin_5f_2g
                                       - 4 * (2 * (1 - 3 * theta2) * (10 + 33 * eta3) + 2 * e2 * (7 - 9 * theta2)
                                              + e4 * (11 - 57 * theta2)) * sin_6f_2g
                                       - 6 * e * (1 - 3 * theta2) * (7 - e2 + 5 * eta3) * sin_7f_2g
                                       - 3 * e2 * (1 - 3 * theta2) * (2 - e2) * sin_8f_2g
                                       )
          + (1 - theta2) * (1 - theta2) * (9 * e2 * (2 - e2) * sin_2f_m_4g
                                                         + 108 * e * (2 - e2) * sin_f_m_4g
                                                         - 12 * (66 - 34 * e2 + 13 * e4) * sin_4g
                                                         - 12 / e * (72 + 2 * e2 + 49 * e4) * sin_f_4g
                                                         - 3 / e2 * (96 + 368 * e2 + 158 * e4 + 161 * e6) * sin_2f_4g
                                                         - 216 / e * (4 + 5 * e4) * sin_3f_4g
                                                         + 144 * (17 - 16 * e2 - e4) * sin_4f_4g
                                                         + 24 / e * (212 - 88 * e2 - 43 * e4) * sin_5f_4g
                                                         + 1 / e2 * (1568 + 4112 * e2 - 3238 * e4 - 93 * e6) * sin_6f_4g
                                                         + 12 / e * (168 + 50 * e2 - 95 * e4) * sin_7f_4g
                                                         + 12 * (82 - 26 * e2 - 11 * e4) * sin_8f_4g
                                                         + 108 * e * (2 - e2) * sin_9f_4g + 9 * e2 * (2 - e2) * sin_10f_4g
                                                         )
          );

    double const delta_g1 = mu2 * J2 / 24 / Gp4 *
          (6 * B20 * (3 / ep * (4 - ep2) * sin_fp + 6 * sin_2fp + ep * sin_3fp)
           + B22 * (3 * ep * sin_fp_m_2gs - 18 * sin_2gs
                    - 3 / ep * (4 + 5 * ep2) * sin_fp_2gs
                    + 1 / ep * (28 - ep2) * sin_3fp_2gs + 18 * sin_4fp_2gs + 3 * ep * sin_5fp_2gs)
           - 18 * (1 - 5 * thetap2) * (fp - lp + ep * sin_fp)
           + 3 * (3 - 5 * thetap2) * (3 * ep * sin_fp_2gs + 3 * sin_2fp_2gs + ep * sin_3fp_2gs)
           );
    double const delta_g2 = - mu4 * J22 / 2048 / G8 *
          (48 * (2 * (5 + 18 * theta2 - 215 * theta4) + e2 * (25 - 126 * theta2 + 45 * theta4)) * (f - l)
           - 96 * (2 * (1 - theta2) * (1 - 15  * theta2) + e2 * (5 - 112  * theta2 + 135 * theta4)) * (f - l) * cos_2g
           + 48 / e * (8 * (1 - 3 * theta2) * (1 - 3 * theta2) * (1 - eta) / e2
                       + 2 * (17 + 22 * theta2 - 191 * theta4)
                       + 4 * e2 * (29 - 105 * theta2 - 33 * theta4)
                       + eta * (1 - 3 * theta2) * (24 * (1 - 5 * theta2) - e2 * (5 - 27 * theta2))) * sin_f
           + 6 / e2 * (16 * (9 - 26 * theta2 + 25 * theta4)
                          + 32 * e2 *(23 - 50 * theta2 - 13 * theta4)
                          + e4 * (173 - 738 * theta2 + 301 * theta4)
                          + 64 * eta * (1 - 3 * theta2) * (1 - 3 * theta2 + e2 * (1 - 6 * theta2))) * sin_2f
           + 8 / e * (4 * (47 - 102 * theta2 + 63 * theta4)
                      + 6 * e2 * (41 - 112 * theta2 + 37 * theta4)
                      + eta * (1 - 3 * theta2) *
                              (64 * (1 - 3 * theta2) + e2 * (5 - 39 * theta2))) * sin_3f
           + 12 * (2 * (39 - 86 * theta2 + 55 * theta4)
                   + e2 * (21 - 66 * theta2 + 37 * theta4)
                   + 16 * eta * (1 - 3 * theta2) * (1 - 3 * theta2)) * sin_4f
           + 24 * e * (2 * (5 - 12 * theta2 + 9 * theta4)
                       + eta * (1 - 3 * theta2) * (1 - 3 * theta2)) * sin_5f
           + 2 * e2 * (11 - 30 * theta2 + 27 * theta4) * sin_6f
           - 12 * e2 * (1 - theta2) * (1 - 3 * theta2) * sin_4f_m_2g
           - 12 * e * (1 - theta2) * (1 - 3 * theta2) * (11 + eta) * sin_3f_m_2g
           - 24 * (1 - theta2) * (2 * (1 - 3 * theta2) * (12 + eta) + e2 * (5 - 39 * theta2)) * sin_2f_m_2g
           - 12 / e * (96 * (1 - theta2) * (1 - 3 * theta2) + e2 * (57 - 532 * theta2 + 459 * theta4)
                       + e2 * eta * (1 - 3 * theta2) * (3 - 11 * theta2)) * sin_f_m_2g
           + 2 / e2 * (96 * (1 - theta2) * (1 - 3 * theta2) * (4 - eta)
                          + 2 * e2 * (169 - 1796 * theta2 + 1467 * theta4)
                          - e4 * (709 - 6764 * theta2 + 8739 * theta4)
                          - 16 * eta * (e2 * (3 + 10 * theta2 - 33 * theta4)
                                        - 6 * e4 * (1 - 24 * theta2 + 30 * theta4))) * sin_2g
           + 24 / e * (8 * (1 - theta2) * (1 - 3 * theta2) * (1 - eta) / e2
                       - 8 * (9 - 44 * theta2 + 39 * theta4)
                       - e2 * (139 - 2020 * theta2 + 2337 * theta4)
                       + eta * (1 - 3 * theta2) * (8 * (7 - 9 * theta2) + e2 * (19 - 39 * theta2))) * sin_f_2g
           - 96 * ((1 - theta2) * (1 - 3 * theta2) * (16 - eta * (16 + 5 * e2)) / e2
                   + (19 - 256 * theta2 + 237 * theta4)
                   - e2 * (1 + 44 * theta2 - 33 * theta4)) * sin_2f_2g
           - 8 / e * (56 * (1 - theta2) * (1 - 3 * theta2) * (1 - eta) / e2
                      + 8 * (3 - 44 * theta2 + 13 * theta4)
                      - e2 * (69 + 324 * theta2 + 79 * theta4)
                      + eta * (1 - 3 * theta2) * (8 * (9 - 23 * theta2) - e2 * (19 - 23 * theta2))) * sin_3f_2g
           + 12 / e2 * (12 * e2 * (7 - 32 * theta2 + 33 * theta4)
                           + e4 * (21 - 76 * theta2 + 159 * theta4)
                           - 16 * eta * (1 - 3 * theta2) * (7 * (1 - theta2) + 3 * e2 * (1 - 2 * theta2))) * sin_4f_2g
           + 12 / e * (e2 * (43 - 252 * theta2 + 225 * theta4)
                       - eta * (1 - 3 * theta2) * (128 * (1 - theta2) + e2 * (7 - 15 * theta2))) * sin_5f_2g
           - 8 * (1 - theta2) * (2 * (1 - 3 * theta2) * (10 + 33 * eta) - 3 * e2 * (3 - 17 * theta2)) * sin_6f_2g
           - 12 * e * (1 - theta2) * (1 - 3 * theta2) * (7 + 5 * eta) * sin_7f_2g
           - 12 * e2 * (1 - theta2) * (1 - 3 * theta2) * sin_8f_2g
           + 9 * e2 * (1 - theta2) * (1 - theta2) * sin_2f_m_4g
           + 108 * e * (1 - theta2) * (1 - theta2) * sin_f_m_4g
           - 6 * (1 - theta2) * (1 - theta2) * (66 - e2) * sin_4g
           - 12 / e * (1 - theta2) * (1 - theta2) * (36 + 19 * e2) * sin_f_4g
           - 3 / e2 * (16 * (1 - theta2) * (1 - theta2) * (3 + 13 * e2) + e4 * (7 + 194 * theta2 - 9 * theta4)) * sin_2f_4g
           - 24 / e * (18 * (1 - theta2) * (1 - theta2) - e2 * (19 - 102 * theta2 + 35 * theta4)) * sin_3f_4g
           + 12 * (2 * (83 - 210 * theta2 + 103 * theta4) + e2 * (11 - 66 * theta2 + 23 * theta4)) * sin_4f_4g
           + 24 / e * (106 * (1 - theta2) * (1 - theta2) + e2 * (53 - 138 * theta2 + 69 * theta4)) * sin_5f_4g
           + 1 / e2 * (16 * (1 - theta2) * (1 - theta2) * (49 + 153 * e2)
                          + e4 * (197 - 538 * theta2 + 277 * theta4)) * sin_6f_4g
           + 12 / e * (1 - theta2) * (1 - theta2) * (84 + 67 * e2) * sin_7f_4g
           + 6 * (1 - theta2) * (1 - theta2) * (82 + 15 * e2) * sin_8f_4g
           + 108 * e * (1 - theta2) * (1 - theta2) * sin_9f_4g
           + 9 * e2 * (1 - theta2) * (1 - theta2) * sin_10f_4g
          );

    double const delta_h1 = (mu2 * J2 / 4 / Gp4) * thetap *
                (-6 * (fp - lp + ep * sin_fp) + 3 * ep * sin_fp_2gs
                 + 3 * sin_2fp_2gs + ep * sin_3fp_2gs
                );
    double const delta_h2 = mu4 * J22 / 128 / G8 * theta *
        (12 * (4 * (1 - 10 * theta2) - e2 * (9 - 5 * theta2) + 2 * e2 * (8 - 15 * theta2) * cos_2g) * (f - l)
         + 12 * e * (3 * (1 - 3 * theta2) * (4 - eta * (4 - e2)) / e2 - (17 + 21 * theta2)) * sin_f
         + 6 * (12 * (1 - 3 * theta2) * (1 - eta) / e2 - (9 - 5 * theta2)) * sin_2f
         + 12 * e * (1 - 3 * theta2) * (1 - eta) * sin_3f
         - 6 * e * (1 - 3 * theta2) * (1 - eta) * sin_f_m_2g
         - e2 * (4 * (7 + 3 * theta2) * (1 - eta) / e2 - 191 + 1053 * theta2 / 2 + 12 * eta * (8 - 15 * theta2)) * sin_2g
         + 6 * e * ((1 - 3 * theta2) * (4 - eta * (4 + 5 * e2)) / e2 + 3 * (35 - 73 * theta2)) * sin_f_2g
         + 24 * (6 * (1 - theta2) + e2 * (4 - theta2)) * sin_2f_2g
         - 2 * e * ((1 - 3 * theta2) * (28 - eta * (28 - e2)) / e2 - (21 + 97 * theta2)) * sin_3f_2g
         - 6 * e2 * (6 * (1 - 3 * theta2) * (1 - eta) / e2 - 13 * e2 * theta2) * sin_4f_2g
         - 6 * e * (1 - 3 * theta2) * (1 - eta) * sin_5f_2g
         - 3 * e2 * (7 + 5 * theta2) * sin_2f_4g
         - 12 * e * (5 + theta2) * sin_3f_4g
         - 3 * (4 * (4 - theta2) + e2 * (7 + theta2)) * sin_4f_4g
         - 12 * e * (3 - theta2) * sin_5f_4g
         - e2 * (7 - 3 * theta2) * sin_6f_4g
        );

    return {delta_l1, delta_g1, delta_h1, delta_L1, delta_G1, 0};
}