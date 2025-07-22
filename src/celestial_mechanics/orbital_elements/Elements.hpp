#pragma once

#include <Core>


struct KeplerianElements {

    double inclination;
    double ascendingNode;

    double semimajor;
    double eccentricity;

    double argumentPeriapsis;
    double trueAnomaly;
};

struct CartesianElements {

    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
};

struct DelaunayElements {

    double l; double g; double h;
    double L; double G; double H;

};

struct EquinoctialElements {

    double I;

    double a;
    double h;
    double k;
    double p;
    double q;
    double lambda;
};