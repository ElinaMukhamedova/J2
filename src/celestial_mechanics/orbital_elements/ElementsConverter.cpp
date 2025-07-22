#include "ElementsConverter.hpp"
#include "AnomalyConverter.hpp"
#include <cmath>
#include <numbers>
#include <Core>
#include <iostream>


Eigen::Matrix<double, 3, 3> R1(double const& phi) {
    double const cos_phi = std::cos(phi);
    double const sin_phi = std::sin(phi);
    return Eigen::Matrix<double, 3, 3>{{1, 0, 0},
                                        {0, cos_phi, sin_phi},
                                        {0, -sin_phi, cos_phi}};
}

Eigen::Matrix<double, 3, 3> R3(double const& phi) {
    double const cos_phi = std::cos(phi);
    double const sin_phi = std::sin(phi);
    return Eigen::Matrix<double, 3, 3>{{cos_phi, sin_phi, 0},
                                        {-sin_phi, cos_phi, 0},
                                        {0, 0, 1}};
}

double N(double const& x) {return x>=0 ? x : x + 2*std::numbers::pi;}

namespace Orbit {

CartesianElements convertKeplerianToCartesian(KeplerianElements const& el, double const mu) {
    double const inclination = el.inclination;
    double const ascendingNode = el.ascendingNode;
    double const semimajor = el.semimajor;
    double const eccentricity = el.eccentricity;
    double const argumentPeriapsis = el.argumentPeriapsis;
    double const trueAnomaly = el.trueAnomaly;

    double const cos_trueAnomaly = std::cos(trueAnomaly);
    double const sin_trueAnomaly = std::sin(trueAnomaly);

    double const p = semimajor * (1 - eccentricity*eccentricity);
    double const r = p / (1 + eccentricity * cos_trueAnomaly);
    double const coef = std::sqrt(mu / p);

    Eigen::Vector3d r_elli{r * cos_trueAnomaly, r * sin_trueAnomaly, 0};
    Eigen::Vector3d v_elli{-coef * sin_trueAnomaly, coef * (eccentricity + cos_trueAnomaly), 0};
    Eigen::Matrix<double, 3, 3> Q = R3(-ascendingNode) * R1(-inclination) * R3(-argumentPeriapsis);
    
    Eigen::Vector3d position = Q * r_elli;
    Eigen::Vector3d velocity = Q * v_elli;

    return CartesianElements{position, velocity};
}

KeplerianElements convertCartesianToKeplerian(CartesianElements const& el, double const mu) {
    const Eigen::Vector3d& position = el.position;
    const Eigen::Vector3d& velocity = el.velocity;
    Eigen::Vector3d const x{1, 0, 0};
    Eigen::Vector3d const y{0, 1, 0};
    Eigen::Vector3d const z{0, 0, 1};
    double const r = position.norm();
    double const v = velocity.norm();

    Eigen::Vector3d const L = position.cross(velocity);
    Eigen::Vector3d const l1 = z.cross(L);
    Eigen::Vector3d const l2 = l1.norm()==0 ? x : l1.normalized();
    Eigen::Vector3d const l = l2.cross(z);
    double const inclination = std::atan2(L.dot(l), L.dot(z));
    double const ascendingNode = N(std::atan2(l2.dot(y), l2.dot(x)));
    
    double const semimajor = mu / (2*mu/r - v*v);
    Eigen::Vector3d const e = (1/mu) * ((v*v - mu/r) * position - position.dot(velocity) * velocity);
    double const eccentricity = e.norm();

    Eigen::Vector3d const e1 = e.norm()==0 ? l2 : e.normalized();
    Eigen::Vector3d const e2 = L.cross(e1).normalized();
    Eigen::Vector3d const g1 = l2;
    Eigen::Vector3d const g2 = L.cross(g1).normalized();
    double const argumentPeriapsis = N(std::atan2(e1.dot(g2), e1.dot(g1)));
    double const trueAnomaly = N(atan2(position.dot(e2), position.dot(e1)));

    return KeplerianElements{inclination, ascendingNode, semimajor,
                                eccentricity, argumentPeriapsis, trueAnomaly};
}

DelaunayElements convertKeplerianToDelaunay(KeplerianElements const& el, double const mu) {
    double const E = TrueToEccentric(el.eccentricity, el.trueAnomaly);
    double const l = EccentricToMean(el.eccentricity, E);
    double const g = el.argumentPeriapsis;
    double const h = el.ascendingNode;
    double const L = std::sqrt(mu*el.semimajor);
    double const G = L*std::sqrt(1-el.eccentricity*el.eccentricity);
    double const H = G*std::cos(el.inclination);

    return DelaunayElements{l, g, h, L, G, H};
}

KeplerianElements convertDelaunayToKeplerian(DelaunayElements const& el, double const mu) {
    double const inclination = std::acos(el.H/el.G);
    double const ascendingNode = el.h;
    double const semimajor = el.L*el.L / mu;
    double const eccentricity = std::sqrt(1-(el.G/el.L)*(el.G/el.L));
    double const argumentPeriapsis = el.g;
    double const E = MeanToEccentric(eccentricity, el.l);
    double const trueAnomaly = EccentricToTrue(eccentricity, E);

    return KeplerianElements{inclination, ascendingNode, semimajor, eccentricity, argumentPeriapsis, trueAnomaly};
}

EquinoctialElements convertDelaunayToEquinoctial(DelaunayElements const& el, double const mu) {
    KeplerianElements const Keplerian = convertDelaunayToKeplerian(el, mu);
    double const I = Keplerian.inclination < std::numbers::pi_v<double>/2 ? 1.0 : -1.0;
    double const a = Keplerian.semimajor;
    double const h = Keplerian.inclination < std::numbers::pi_v<double>/2 ? Keplerian.eccentricity*std::sin(el.g+el.h) : Keplerian.eccentricity*std::sin(el.g-el.h);
    double const k = Keplerian.inclination < std::numbers::pi_v<double>/2 ? Keplerian.eccentricity*std::cos(el.g+el.h) : Keplerian.eccentricity*std::cos(el.g-el.h);
    double const p = Keplerian.inclination < std::numbers::pi_v<double>/2 ? std::tan(Keplerian.inclination/2)*std::sin(el.h) : 1/std::tan(Keplerian.inclination/2)*std::sin(el.h);
    double const q = Keplerian.inclination < std::numbers::pi_v<double>/2 ? std::tan(Keplerian.inclination/2)*std::cos(el.h) : 1/std::tan(Keplerian.inclination/2)*std::cos(el.h);
    double const lambda = Keplerian.inclination < std::numbers::pi_v<double>/2 ? el.l+el.g+el.h : el.l+el.g-el.h;
    
    return EquinoctialElements{I, a, h, k, p ,q, lambda};
}

KeplerianElements convertEquinoctialToKeplerian(EquinoctialElements const& el, double const mu) {
    double const semimajor = el.a;
    double const eccentricity = std::sqrt(el.h*el.h+el.k*el.k);
    double const inclination = std::numbers::pi_v<double>*(1-el.I)/2+2*el.I*std::atan(std::sqrt(el.p*el.p+el.q*el.q));
    double const ascendingNode = std::atan2(el.p, el.q);
    double const xi = std::atan2(el.h, el.k);
    double const argumentPeriapsis = xi - el.I*ascendingNode;
    double const M = el.lambda - xi;
    double const E = MeanToEccentric(eccentricity, M);
    double const trueAnomaly = EccentricToTrue(eccentricity, E);

    return KeplerianElements{inclination, ascendingNode, semimajor, eccentricity, argumentPeriapsis, trueAnomaly};
}

DelaunayElements convertEquinoctialToDelaunay(EquinoctialElements const& el, double const mu) {
    return convertKeplerianToDelaunay(convertEquinoctialToKeplerian(el, mu), mu);
}

};