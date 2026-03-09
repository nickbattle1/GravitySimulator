#include "Scenarios.hpp"

// All positions in metres, velocities in m/s — real SI units.
// Display radii are in pixels and only affect rendering.

namespace Scenarios {

void loadSolarSystem(GravitySimulator& sim) {
    sim.clear();

    // Sun
    sim.addBody(std::make_shared<CelestialBody>(
        "Sun", 1.989e30, Vec2(0, 0), Vec2(0, 0),
        sf::Color::Yellow, 20.0f));

    // Mercury
    sim.addBody(std::make_shared<CelestialBody>(
        "Mercury", 3.3011e23, Vec2(57.9e9, 0), Vec2(0, 47.4e3),
        sf::Color(169, 169, 169), 3.5f));

    // Venus
    sim.addBody(std::make_shared<CelestialBody>(
        "Venus", 4.8675e24, Vec2(108.2e9, 0), Vec2(0, 35.0e3),
        sf::Color(255, 198, 73), 6.0f));

    // Earth
    sim.addBody(std::make_shared<CelestialBody>(
        "Earth", 5.972e24, Vec2(149.6e9, 0), Vec2(0, 29.8e3),
        sf::Color(100, 149, 237), 7.0f));

    // Moon
    sim.addBody(std::make_shared<CelestialBody>(
        "Moon", 7.342e22, Vec2(149.6e9 + 0.384e9, 0), Vec2(0, 29.8e3 + 1.022e3),
        sf::Color(200, 200, 200), 2.5f));

    // Mars
    sim.addBody(std::make_shared<CelestialBody>(
        "Mars", 6.4171e23, Vec2(227.9e9, 0), Vec2(0, 24.1e3),
        sf::Color(188, 39, 50), 5.0f));

    // Jupiter
    sim.addBody(std::make_shared<CelestialBody>(
        "Jupiter", 1.8982e27, Vec2(778.5e9, 0), Vec2(0, 13.1e3),
        sf::Color(255, 140, 0), 14.0f));

    // Saturn
    sim.addBody(std::make_shared<CelestialBody>(
        "Saturn", 5.6834e26, Vec2(1.434e12, 0), Vec2(0, 9.7e3),
        sf::Color(210, 180, 140), 12.0f));

    // Uranus
    sim.addBody(std::make_shared<CelestialBody>(
        "Uranus", 8.6810e25, Vec2(2.871e12, 0), Vec2(0, 6.8e3),
        sf::Color(173, 216, 230), 9.0f));

    // Neptune
    sim.addBody(std::make_shared<CelestialBody>(
        "Neptune", 1.02413e26, Vec2(4.495e12, 0), Vec2(0, 5.4e3),
        sf::Color(63, 84, 186), 8.5f));

    sim.setTimeStep(1800.0); // 30 minutes — finer step keeps Moon orbit stable
}

void loadTwoBody(GravitySimulator& sim) {
    sim.clear();

    sim.addBody(std::make_shared<CelestialBody>(
        "Star", 1.0e30, Vec2(0, 0), Vec2(0, 0),
        sf::Color::Yellow, 15.0f));

    sim.addBody(std::make_shared<CelestialBody>(
        "Planet", 1.0e25, Vec2(100e9, 0), Vec2(0, 25.0e3),
        sf::Color(100, 149, 237), 6.0f));

    sim.setTimeStep(3600.0);
}

void loadBinaryStars(GravitySimulator& sim) {
    sim.clear();

    sim.addBody(std::make_shared<CelestialBody>(
        "Star1", 1.5e30, Vec2(-50e9, 0), Vec2(0, -10.0e3),
        sf::Color(255, 200, 100), 12.0f));

    sim.addBody(std::make_shared<CelestialBody>(
        "Star2", 1.2e30, Vec2(50e9, 0), Vec2(0, 12.5e3),
        sf::Color(100, 200, 255), 10.0f));

    sim.addBody(std::make_shared<CelestialBody>(
        "Planet", 2.0e24, Vec2(0, 200e9), Vec2(18.0e3, 0),
        sf::Color::Green, 5.0f));

    sim.setTimeStep(3600.0);
}

} // namespace Scenarios