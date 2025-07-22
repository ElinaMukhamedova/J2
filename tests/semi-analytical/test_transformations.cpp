#include <gtest/gtest.h>
#include "semi-analytical/transformations.hpp"
#include "celestial_mechanics/orbital_elements/AnomalyConverter.hpp"
#include "celestial_mechanics/orbital_elements/ElementsConverter.hpp"
#include <iostream>
#include <numbers>
#include <algorithm>


TEST(TranformationsTest, MeanToOsculatingToMeanWorks) {
    double const alpha = 6378137;
    double const mu = 3986004418e5;

    DelaunayElements const mean_origin{0.0199672, 0.628319, 0.523599,7.93022e+11, 3.4567e+11,1.72835e+11};
    double const mean_meanAnomaly = mean_origin.l;
    double const mean_G = mean_origin.G;
    double const mean_L = mean_origin.L;
    double const mean_eccentricity = std::sqrt(1 - mean_G*mean_G/mean_L/mean_L);
    double const mean_eccentricAnomaly = Orbit::MeanToEccentric(mean_eccentricity, mean_meanAnomaly);
    double const mean_trueAnomaly = Orbit::EccentricToTrue(mean_eccentricity, mean_eccentricAnomaly);

    DelaunayElements const osc = MeanToOsculating(mean_origin, mean_trueAnomaly, alpha, mu);
    double const osc_meanAnomaly = osc.l;
    double const osc_G = osc.G;
    double const osc_L = osc.L;
    double const osc_eccentricity = std::sqrt(1 - osc_G*osc_G/osc_L/osc_L);
    double const osc_eccentricAnomaly = Orbit::MeanToEccentric(osc_eccentricity, osc_meanAnomaly);
    double const osc_trueAnomaly = Orbit::EccentricToTrue(osc_eccentricity, osc_eccentricAnomaly);

    DelaunayElements const mean_result = OsculatingToMean(osc, osc_trueAnomaly, alpha, mu);

    std::cout << mean_origin.l << " " << mean_result.l << std::endl;
    std::cout << mean_origin.g << " " << mean_result.g << std::endl;
    std::cout << mean_origin.h << " " << mean_result.h << std::endl;
    std::cout << mean_origin.L << " " << mean_result.L << std::endl;
    std::cout << mean_origin.G << " " << mean_result.G << std::endl;
    std::cout << mean_origin.H << " " << mean_result.H << std::endl;

    ASSERT_NEAR((mean_result.l-mean_origin.l)/mean_origin.l, 0, 1e-5);
    ASSERT_NEAR((mean_result.g-mean_origin.g)/mean_origin.g, 0, 1e-5);
    ASSERT_NEAR((mean_result.h-mean_origin.h)/mean_origin.h, 0, 1e-5);
    ASSERT_NEAR((mean_result.L-mean_origin.L)/mean_origin.L, 0, 1e-5);
    ASSERT_NEAR((mean_result.G-mean_origin.G)/mean_origin.G, 0, 1e-5);
    ASSERT_NEAR((mean_result.H-mean_origin.H)/mean_origin.H, 0, 1e-5);
}

TEST(TransformationsTest, OsculatingToMeanToOsculatingWorks) {
    double const alpha = 6378137;
    double const mu = 3986004418e5;

    DelaunayElements const osc_origin{0.0199672, 0.628319, 0.523599,7.93022e+11, 3.4567e+11,1.72835e+11};
    double const osc_meanAnomaly = osc_origin.l;
    double const osc_G = osc_origin.G;
    double const osc_L = osc_origin.L;
    double const osc_eccentricity = std::sqrt(1 - osc_G*osc_G/osc_L/osc_L);
    double const osc_eccentricAnomaly = Orbit::MeanToEccentric(osc_eccentricity, osc_meanAnomaly);
    double const osc_trueAnomaly = Orbit::EccentricToTrue(osc_eccentricity, osc_eccentricAnomaly);

    DelaunayElements const mean = OsculatingToMean(osc_origin, osc_trueAnomaly, alpha, mu);
    double const mean_meanAnomaly = mean.l;
    double const mean_G = mean.G;
    double const mean_L = mean.L;
    double const mean_eccentricity = std::sqrt(1 - mean_G*mean_G/mean_L/mean_L);
    double const mean_eccentricAnomaly = Orbit::MeanToEccentric(mean_eccentricity, mean_meanAnomaly);
    double const mean_trueAnomaly = Orbit::EccentricToTrue(mean_eccentricity, mean_eccentricAnomaly);

    DelaunayElements const osc_result = MeanToOsculating(mean, mean_trueAnomaly, alpha, mu);

    std::cout << osc_origin.l << " " << osc_result.l << std::endl;
    std::cout << osc_origin.g << " " << osc_result.g << std::endl;
    std::cout << osc_origin.h << " " << osc_result.h << std::endl;
    std::cout << osc_origin.L << " " << osc_result.L << std::endl;
    std::cout << osc_origin.G << " " << osc_result.G << std::endl;
    std::cout << osc_origin.H << " " << osc_result.H << std::endl;

    ASSERT_NEAR((osc_result.l-osc_origin.l)/osc_origin.l, 0, 1e-6);
    ASSERT_NEAR((osc_result.g-osc_origin.g)/osc_origin.g, 0, 1e-6);
    ASSERT_NEAR((osc_result.h-osc_origin.h)/osc_origin.h, 0, 1e-6);
    ASSERT_NEAR((osc_result.L-osc_origin.L)/osc_origin.L, 0, 1e-6);
    ASSERT_NEAR((osc_result.G-osc_origin.G)/osc_origin.G, 0, 1e-6);
    ASSERT_NEAR((osc_result.H-osc_origin.H)/osc_origin.H, 0, 1e-6);
}

double relError(double a, double b) {
    return (std::abs(b) > 1e-12) ? (a - b) / b : a - b;
}
double angle_diff(double result, double origin) { return std::fmod(std::abs(result - origin), 2 * std::numbers::pi); }
TEST(TransformationsTest, InclinationRange_200kmPerigeeHight_500kmApogeeHight_Works) {
    double const alpha = 6378137;
    double const mu = 3986004418e5;

    double const perigee = 200000 + alpha;
    double const apogee = 500000 + alpha;
    double const semiMajor = 0.5 * (perigee + apogee);
    double const eccentricity = (apogee - perigee) / (apogee + perigee);

    double const number = 5;
    double const delta_inc = std::numbers::pi / number;
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;

    for (int i = 1; i < number; ++i) {
        double const inclination = i*delta_inc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semiMajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const osc = MeanToOsculating(origin, trueAnomaly, alpha, mu);
        double const osc_meanAnomaly = osc.l;
        double const osc_G = osc.G;
        double const osc_L = osc.L;
        double const e2 = 1 - osc_G*osc_G/osc_L/osc_L;
        double const osc_eccentricity = std::sqrt(std::max(1e-13, e2));
        double const osc_eccentricAnomaly = Orbit::MeanToEccentric(osc_eccentricity, osc_meanAnomaly);
        double const osc_trueAnomaly = Orbit::EccentricToTrue(osc_eccentricity, osc_eccentricAnomaly);
        DelaunayElements const mean_result = OsculatingToMean(osc, osc_trueAnomaly, alpha, mu);

        ASSERT_LT(angle_diff(mean_result.l, origin.l), 3e-2);
        ASSERT_LT(angle_diff(mean_result.g, origin.g), 3e-2);
        ASSERT_LT(angle_diff(mean_result.h, origin.h), 3e-2);
        ASSERT_NEAR(relError(mean_result.L, origin.L), 0, 1e-3);
        ASSERT_NEAR(relError(mean_result.G, origin.G), 0, 1e-3);
        ASSERT_NEAR(relError(mean_result.H, origin.H), 0, 1e-3);

        // DelaunayElements const mean = OsculatingToMean(origin, trueAnomaly, alpha, mu);
        // double const mean_meanAnomaly = mean.l;
        // double const mean_G = mean.G;
        // double const mean_L = mean.L;
        // double const mean_e2 = 1 - mean_G*mean_G/mean_L/mean_L;
        // double const mean_eccentricity = std::sqrt(std::max(1e-13, mean_e2));
        // double const mean_eccentricAnomaly = Orbit::MeanToEccentric(mean_eccentricity, mean_meanAnomaly);
        // double const mean_trueAnomaly = Orbit::EccentricToTrue(mean_eccentricity, mean_eccentricAnomaly);
        // DelaunayElements const osc_result = MeanToOsculating(mean, mean_trueAnomaly, alpha, mu);
        //
        // ASSERT_LT(angle_diff(osc_result.l, origin.l), 3e-2);
        // ASSERT_LT(angle_diff(osc_result.g, origin.g), 3e-2);
        // ASSERT_LT(angle_diff(osc_result.h, origin.h), 3e-2);
        // ASSERT_NEAR(relError(osc_result.L, origin.L), 0, 1e-3);
        // ASSERT_NEAR(relError(osc_result.G, origin.G), 0, 1e-3);
        // ASSERT_NEAR(relError(osc_result.H, origin.H), 0, 1e-3);
    }
}

TEST(TransformationsTest, InclinationRangeWorks) {
    double const alpha = 6378137;
    double const mu = 3986004418e5;

    double const eccentricity = 0.9;
    double const perigee = 200000;
    double const semimajor = perigee / (1 - eccentricity);

    double const number = 5;
    double const delta_inc = std::numbers::pi / number;
    double const ascendingNode = std::numbers::pi/8;
    double const argumentPeriapsis = std::numbers::pi/7;
    double const trueAnomaly = std::numbers::pi/6;

    for (int i = 1; i < number; ++i) {
        double const inclination = i * delta_inc;
        KeplerianElements const originKeplerian{inclination, ascendingNode, semimajor, eccentricity, argumentPeriapsis, trueAnomaly};
        DelaunayElements const origin = Orbit::convertKeplerianToDelaunay(originKeplerian, mu);

        DelaunayElements const osc = MeanToOsculating(origin, trueAnomaly, alpha, mu);
        double const osc_meanAnomaly = osc.l;
        double const osc_G = osc.G;
        double const osc_L = osc.L;
        double const e2 = 1 - osc_G*osc_G/osc_L/osc_L;
        double const osc_eccentricity = std::sqrt(std::max(1e-13, e2));
        double const osc_eccentricAnomaly = Orbit::MeanToEccentric(osc_eccentricity, osc_meanAnomaly);
        double const osc_trueAnomaly = Orbit::EccentricToTrue(osc_eccentricity, osc_eccentricAnomaly);
        DelaunayElements const mean_result = OsculatingToMean(osc, osc_trueAnomaly, alpha, mu);

        ASSERT_LT(angle_diff(mean_result.l, origin.l), 0.5);
        ASSERT_LT(angle_diff(mean_result.g, origin.g), 0.5);
        ASSERT_LT(angle_diff(mean_result.h, origin.h), 0.5);
        ASSERT_NEAR(relError(mean_result.L, origin.L), 0, 0.5);
        ASSERT_NEAR(relError(mean_result.G, origin.G), 0, 0.5);
        ASSERT_NEAR(relError(mean_result.H, origin.H), 0, 0.5);
    }
}