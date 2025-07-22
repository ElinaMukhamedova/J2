#include <gtest/gtest.h>
#include "celestial_mechanics/orbital_elements/Elements.hpp"
#include "celestial_mechanics/orbital_elements/ElementsConverter.hpp"
#include <Core>
#include <numbers>
//#include <iostream>


TEST(ElementsConverterTest, KeplerianToCartesianToKeplerianWorks) {
    KeplerianElements const origin{std::numbers::pi/3, std::numbers::pi/6, 6378137, 0.9, std::numbers::pi/5, std::numbers::pi/4};
    CartesianElements const transfer = Orbit::convertKeplerianToCartesian(origin, 3.986e14);
    KeplerianElements const result = Orbit::convertCartesianToKeplerian(transfer, 3.986e14);
    ASSERT_NEAR(std::abs(origin.inclination-result.inclination)/origin.inclination, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.ascendingNode-result.ascendingNode)/origin.ascendingNode, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.semimajor-result.semimajor)/origin.semimajor, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.eccentricity-result.eccentricity)/origin.eccentricity, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.argumentPeriapsis-result.argumentPeriapsis)/origin.argumentPeriapsis, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.trueAnomaly-result.trueAnomaly)/origin.trueAnomaly, 0, 1e-14);
}

TEST(ElementsConverterTest, CartesianToKeplerianToCartesianWorks) {
    Eigen::Vector3d position{6400000, 6500000, 6600000};
    Eigen::Vector3d velocity{100, -5000, -1000};
    CartesianElements const origin{position, velocity};
    KeplerianElements const transfer = Orbit::convertCartesianToKeplerian(origin, 3.986e14);
    CartesianElements const result = Orbit::convertKeplerianToCartesian(transfer, 3.986e14);
    ASSERT_NEAR(std::abs(origin.position(0)-result.position(0))/origin.position(0), 0, 1e-13);  
    ASSERT_NEAR(std::abs(origin.position(1)-result.position(1))/origin.position(1), 0, 1e-13);
    ASSERT_NEAR(std::abs(origin.position(2)-result.position(2))/origin.position(2), 0, 1e-13);
    ASSERT_NEAR(std::abs(origin.velocity(0)-result.velocity(0))/origin.velocity(0), 0, 1e-13);
    ASSERT_NEAR(std::abs(origin.velocity(1)-result.velocity(1))/origin.velocity(1), 0, 1e-13);
    ASSERT_NEAR(std::abs(origin.velocity(2)-result.velocity(2))/origin.velocity(2), 0, 1e-13);
}

TEST(ElementsConverterTest, KeplerianToDelaunayToKeplerianWorks) {
    KeplerianElements const origin{std::numbers::pi/3, std::numbers::pi/6, 6378137, 0.9, std::numbers::pi/5, std::numbers::pi/4};
    DelaunayElements const transfer = Orbit::convertKeplerianToDelaunay(origin, 3.986e14);
    KeplerianElements const result = Orbit::convertDelaunayToKeplerian(transfer, 3.986e14);
    //std::cout << transfer.l << '\n';
    //std::cout << transfer.L << '\n';
    //std::cout << transfer.g << '\n';
    //std::cout << transfer.G << '\n';
    //std::cout << transfer.h << '\n';
    //std::cout << transfer.H << '\n';
    ASSERT_NEAR(std::abs(origin.inclination-result.inclination)/origin.inclination, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.ascendingNode-result.ascendingNode)/origin.ascendingNode, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.semimajor-result.semimajor)/origin.semimajor, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.eccentricity-result.eccentricity)/origin.eccentricity, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.argumentPeriapsis-result.argumentPeriapsis)/origin.argumentPeriapsis, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.trueAnomaly-result.trueAnomaly)/origin.trueAnomaly, 0, 1e-14);
}

TEST(ElementsConverterTest, DelaunayToKeplerianToDelaunayWorks) {
    DelaunayElements const origin{0.0199672, 0.628319, 0.523599,7.93022e+11, 3.4567e+11,1.72835e+11};
    KeplerianElements const transfer = Orbit::convertDelaunayToKeplerian(origin, 3.986e14);
    DelaunayElements const result = Orbit::convertKeplerianToDelaunay(transfer, 3.986e14);
    ASSERT_NEAR(std::abs(origin.l-result.l)/origin.l, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.L-result.L)/origin.L, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.g-result.g)/origin.g, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.G-result.G)/origin.G, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.h-result.h)/origin.h, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.H-result.H)/origin.H, 0, 1e-14);
}

TEST(ElementsConverterTest, DelaunayToEquinoctialToDelaunayWorks) {
    DelaunayElements const origin{0.0199672, 0.628319, 0.523599,7.93022e+11, 3.4567e+11,1.72835e+11};
    EquinoctialElements const transfer = Orbit::convertDelaunayToEquinoctial(origin, 3.986e14);
    //std::cout << transfer.I << ", "<< transfer.a << ", " << transfer.h << ", " << transfer.k << ", " << transfer.p << ", " << transfer.q << ", " << transfer.lambda << '\n';
    DelaunayElements const result = Orbit::convertEquinoctialToDelaunay(transfer, 3.986e14);
    ASSERT_NEAR(std::abs(origin.l-result.l)/origin.l, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.L-result.L)/origin.L, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.g-result.g)/origin.g, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.G-result.G)/origin.G, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.h-result.h)/origin.h, 0, 1e-14);
    ASSERT_NEAR(std::abs(origin.H-result.H)/origin.H, 0, 1e-14);
}

TEST(ElementsConverterTest, EquinoctialToDelaunayToEquinoctialWorks) {
    EquinoctialElements origin{1, 6.37813e+06, 0.822191, 0.366062, 0.288675, 0.5, 1.17189};
    DelaunayElements transfer = Orbit::convertEquinoctialToDelaunay(origin, 3.986e14);
    EquinoctialElements result = Orbit::convertDelaunayToEquinoctial(transfer, 3.986e14);
    ASSERT_NEAR(std::abs(result.I-origin.I)/origin.I, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.a-origin.a)/origin.a, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.h-origin.h)/origin.h, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.k-origin.k)/origin.k, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.p-origin.p)/origin.p, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.q-origin.q)/origin.q, 0, 1e-14);
    ASSERT_NEAR(std::abs(result.lambda-origin.lambda)/origin.lambda, 0, 1e-14);
}