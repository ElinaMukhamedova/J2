#pragma once

#include "Elements.hpp"


namespace Orbit {

CartesianElements convertKeplerianToCartesian(KeplerianElements const& el, double const mu);

KeplerianElements convertCartesianToKeplerian(CartesianElements const& el, double const mu);

DelaunayElements convertKeplerianToDelaunay(KeplerianElements const& el, double const mu);

KeplerianElements convertDelaunayToKeplerian(DelaunayElements const& el, double const mu);

EquinoctialElements convertDelaunayToEquinoctial(DelaunayElements const& el, double const mu);

KeplerianElements convertEquinoctialToKeplerian(EquinoctialElements const& el, double const mu);

DelaunayElements convertEquinoctialToDelaunay(EquinoctialElements const& el, double const mu);

};