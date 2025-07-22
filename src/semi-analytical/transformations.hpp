#pragma once

#include "celestial_mechanics/orbital_elements/Elements.hpp"


DelaunayElements OsculatingToMean(DelaunayElements const& osc, double const f_prime, double const alpha, double const mu);
DelaunayElements MeanToOsculating(DelaunayElements const& mean, double const f_prime, double const alpha, double const mu);