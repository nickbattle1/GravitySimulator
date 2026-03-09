#pragma once

#include "GravitySimulator.hpp"

namespace Scenarios {
    void loadSolarSystem(GravitySimulator& sim);
    void loadTwoBody(GravitySimulator& sim);
    void loadBinaryStars(GravitySimulator& sim);
}