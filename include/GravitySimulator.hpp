#pragma once

#include "CelestialBody.hpp"
#include <vector>
#include <memory>

enum class IntegrationMethod {
    EULER,
    VELOCITY_VERLET,
    RK4
};

class GravitySimulator {
public:
    void update();

    // Body management
    void addBody(std::shared_ptr<CelestialBody> body);
    bool removeBody(const std::string& name);
    std::shared_ptr<CelestialBody> getBody(const std::string& name) const;
    const std::vector<std::shared_ptr<CelestialBody>>& getBodies() const { return m_bodies; }
    void clear();

    // Configuration
    void setIntegrationMethod(IntegrationMethod method) { m_method = method; }
    IntegrationMethod getIntegrationMethod() const { return m_method; }
    void setTimeStep(double dt) { m_timeStep = dt; }
    double getTimeStep() const { return m_timeStep; }
    void setCollisionDetection(bool enabled) { m_detectCollisions = enabled; }
    bool getCollisionDetection() const { return m_detectCollisions; }
    void setSofteningLength(double eps) { m_softening = eps; }

private:
    std::vector<std::shared_ptr<CelestialBody>> m_bodies;
    IntegrationMethod m_method = IntegrationMethod::VELOCITY_VERLET;
    double m_timeStep = 3600.0; // seconds (1 hour default)
    bool m_detectCollisions = false;
    double m_softening = 1e6; // softening length in metres to prevent singularities

    // Gravitational constant (real SI value — all physics in real units)
    static constexpr double G = 6.67430e-11;

    // Compute acceleration on body i due to all other bodies
    Vec2 computeAcceleration(size_t i,
                             const std::vector<Vec2>& positions,
                             const std::vector<double>& masses) const;

    // Compute all accelerations at once
    std::vector<Vec2> computeAllAccelerations(
        const std::vector<Vec2>& positions,
        const std::vector<double>& masses) const;

    void integrateEuler(double dt);
    void integrateVelocityVerlet(double dt);
    void integrateRK4(double dt);
    void handleCollisions();
};