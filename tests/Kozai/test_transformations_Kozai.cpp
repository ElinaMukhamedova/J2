#include <gtest/gtest.h>
#include "Kozai/transformations_Kozai.hpp"
#include "celestial_mechanics/orbital_elements/ElementsConverter.hpp"
#include <numbers>


TEST(TransformationsKozaiTest, MeanToOsculatingToMean_InclinationRange_Works) {
    double const mu = 3986004418e5;
    double const J2 = 0.001082636;

    double const R = 6378137;
    double const perigee = 200000 + R;
    double const apogee = 500000 + R;
    double const semiMajor = 0.5 * (perigee + apogee);
    double const eccentricity = (apogee - perigee) / (apogee + perigee);

    double const number = 180;
    double const delta_inc = std::numbers::pi / number;
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;

    for (int i = 0; i < number; ++i) {
        double const inclination = i*delta_inc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semiMajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const mean_origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const osc = MeanToOsculating_Kozai(mean_origin, mu, J2);
        DelaunayElements const mean_result = OsculatingToMean_Kozai(osc, mu, J2);

        std::cout << "eccentricity = " << eccentricity << ", inclination = " << inclination / std::numbers::pi * 180 << " degree(s)" << std::endl;
        std::cout << mean_origin.l << " " << mean_result.l << std::endl;
        std::cout << mean_origin.g << " " << mean_result.g << std::endl;
        std::cout << mean_origin.h << " " << mean_result.h << std::endl;
        std::cout << mean_origin.L << " " << mean_result.L << std::endl;
        std::cout << mean_origin.G << " " << mean_result.G << std::endl;
        std::cout << mean_origin.H << " " << mean_result.H << std::endl;
        std::cout << std::endl;

        ASSERT_NEAR((mean_result.l-mean_origin.l)/mean_origin.l, 0, 1e-14);
        ASSERT_NEAR((mean_result.g-mean_origin.g)/mean_origin.g, 0, 1e-14);
        ASSERT_NEAR((mean_result.h-mean_origin.h)/mean_origin.h, 0, 1e-14);
        ASSERT_NEAR((mean_result.L-mean_origin.L)/mean_origin.L, 0, 1e-14);
        ASSERT_NEAR((mean_result.G-mean_origin.G)/mean_origin.G, 0, 1e-14);
        ASSERT_NEAR((mean_result.H-mean_origin.H)/mean_origin.H, 0, 1e-14);
    }
}

TEST(TransformationsKozaiTest, OsculatingToMeanToOsculating_InclinationRange_Works) {
    double const mu = 3986004418e5;
    double const J2 = 0.001082636;

    double const R = 6378137;
    double const perigee = 200000 + R;
    double const apogee = 500000 + R;
    double const semiMajor = 0.5 * (perigee + apogee);
    double const eccentricity = (apogee - perigee) / (apogee + perigee);

    double const number = 180;
    double const delta_inc = std::numbers::pi / number;
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;

    for (int i = 0; i < number; ++i) {
        double const inclination = i*delta_inc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semiMajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const osculating_origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const mean = OsculatingToMean_Kozai(osculating_origin, mu, J2);
        DelaunayElements const osculating_result = MeanToOsculating_Kozai(mean, mu, J2);

        std::cout << "eccentricity = " << eccentricity << ", inclination = " << inclination / std::numbers::pi * 180 << " degree(s)" << std::endl;
        std::cout << osculating_origin.l << " " << osculating_result.l << std::endl;
        std::cout << osculating_origin.g << " " << osculating_result.g << std::endl;
        std::cout << osculating_origin.h << " " << osculating_result.h << std::endl;
        std::cout << osculating_origin.L << " " << osculating_result.L << std::endl;
        std::cout << osculating_origin.G << " " << osculating_result.G << std::endl;
        std::cout << osculating_origin.H << " " << osculating_result.H << std::endl;
        std::cout << std::endl;

        ASSERT_NEAR((osculating_result.l-osculating_origin.l)/osculating_origin.l, 0, 1e-14);
        ASSERT_NEAR((osculating_result.g-osculating_origin.g)/osculating_origin.g, 0, 1e-14);
        ASSERT_NEAR((osculating_result.h-osculating_origin.h)/osculating_origin.h, 0, 1e-14);
        ASSERT_NEAR((osculating_result.L-osculating_origin.L)/osculating_origin.L, 0, 1e-14);
        ASSERT_NEAR((osculating_result.G-osculating_origin.G)/osculating_origin.G, 0, 1e-14);
        ASSERT_NEAR((osculating_result.H-osculating_origin.H)/osculating_origin.H, 0, 1e-14);
    }
}

TEST(TransformationsKozaiTest, MeanToOsculatingToMean_EccentricityRange_Works) {
    double const mu = 3986004418e5;
    double const J2 = 0.001082636;

    double const number = 9;
    double const delta_ecc = 1 / (number + 1);
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;
    double const inclination = 51.0 / 180.0 * std::numbers::pi;
    double const semiMajor = 40000000;

    for (int i = 1; i < number + 1; ++i) {
        double const eccentricity = i*delta_ecc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semiMajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const mean_origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const osc = MeanToOsculating_Kozai(mean_origin, mu, J2);
        DelaunayElements const mean_result = OsculatingToMean_Kozai(osc, mu, J2);

        std::cout << "eccentricity = " << eccentricity << ", inclination = " << inclination / std::numbers::pi * 180 << " degree(s)" << std::endl;
        std::cout << mean_origin.l << " " << mean_result.l << std::endl;
        std::cout << mean_origin.g << " " << mean_result.g << std::endl;
        std::cout << mean_origin.h << " " << mean_result.h << std::endl;
        std::cout << mean_origin.L << " " << mean_result.L << std::endl;
        std::cout << mean_origin.G << " " << mean_result.G << std::endl;
        std::cout << mean_origin.H << " " << mean_result.H << std::endl;
        std::cout << std::endl;

        ASSERT_NEAR((mean_result.l-mean_origin.l)/mean_origin.l, 0, 1e-15);
        ASSERT_NEAR((mean_result.g-mean_origin.g)/mean_origin.g, 0, 1e-15);
        ASSERT_NEAR((mean_result.h-mean_origin.h)/mean_origin.h, 0, 1e-15);
        ASSERT_NEAR((mean_result.L-mean_origin.L)/mean_origin.L, 0, 1e-15);
        ASSERT_NEAR((mean_result.G-mean_origin.G)/mean_origin.G, 0, 1e-15);
        ASSERT_NEAR((mean_result.H-mean_origin.H)/mean_origin.H, 0, 1e-15);
    }
}

TEST(TransformationsKozaiTest, OsculatingToMeanToOsculating_EccentricityRange_Works) {
    double const mu = 3986004418e5;
    double const J2 = 0.001082636;

    double const number = 9;
    double const delta_ecc = 1 / (number + 1);
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;
    double const inclination = 51.0 / 180.0 * std::numbers::pi;
    double const semiMajor = 40000000;

    for (int i = 1; i < number + 1; ++i) {
        double const eccentricity = i*delta_ecc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semiMajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const osculating_origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const mean = OsculatingToMean_Kozai(osculating_origin, mu, J2);
        DelaunayElements const osculating_result = MeanToOsculating_Kozai(mean, mu, J2);

        std::cout << "eccentricity = " << eccentricity << ", inclination = " << inclination / std::numbers::pi * 180 << " degree(s)" << std::endl;
        std::cout << osculating_origin.l << " " << osculating_result.l << std::endl;
        std::cout << osculating_origin.g << " " << osculating_result.g << std::endl;
        std::cout << osculating_origin.h << " " << osculating_result.h << std::endl;
        std::cout << osculating_origin.L << " " << osculating_result.L << std::endl;
        std::cout << osculating_origin.G << " " << osculating_result.G << std::endl;
        std::cout << osculating_origin.H << " " << osculating_result.H << std::endl;
        std::cout << std::endl;

        ASSERT_NEAR((osculating_result.l-osculating_origin.l)/osculating_origin.l, 0, 1e-15);
        ASSERT_NEAR((osculating_result.g-osculating_origin.g)/osculating_origin.g, 0, 1e-15);
        ASSERT_NEAR((osculating_result.h-osculating_origin.h)/osculating_origin.h, 0, 1e-15);
        ASSERT_NEAR((osculating_result.L-osculating_origin.L)/osculating_origin.L, 0, 1e-15);
        ASSERT_NEAR((osculating_result.G-osculating_origin.G)/osculating_origin.G, 0, 1e-15);
        ASSERT_NEAR((osculating_result.H-osculating_origin.H)/osculating_origin.H, 0, 1e-15);
    }
}