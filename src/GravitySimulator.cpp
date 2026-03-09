#include "GravitySimulator.hpp"
#include <algorithm>
#include <cmath>

// ---------- acceleration helpers ----------

Vec2 GravitySimulator::computeAcceleration(
        size_t i,
        const std::vector<Vec2>& positions,
        const std::vector<double>& masses) const {
    Vec2 acc{0.0, 0.0};
    const double eps2 = m_softening * m_softening;

    for (size_t j = 0; j < positions.size(); ++j) {
        if (j == i) continue;
        Vec2 r = positions[j] - positions[i];
        double dist2 = r.lengthSquared() + eps2;
        double dist  = std::sqrt(dist2);
        acc += r * (G * masses[j] / (dist2 * dist));
    }
    return acc;
}

std::vector<Vec2> GravitySimulator::computeAllAccelerations(
        const std::vector<Vec2>& positions,
        const std::vector<double>& masses) const {
    std::vector<Vec2> accs(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        accs[i] = computeAcceleration(i, positions, masses);
    }
    return accs;
}

// ---------- integrators ----------

void GravitySimulator::integrateEuler(double dt) {
    const size_t n = m_bodies.size();

    // Snapshot current state
    std::vector<Vec2> positions(n), velocities(n);
    std::vector<double> masses(n);
    for (size_t i = 0; i < n; ++i) {
        positions[i]  = m_bodies[i]->getPosition();
        velocities[i] = m_bodies[i]->getVelocity();
        masses[i]     = m_bodies[i]->getMass();
    }

    auto accs = computeAllAccelerations(positions, masses);

    for (size_t i = 0; i < n; ++i) {
        positions[i]  += velocities[i] * dt;
        velocities[i] += accs[i] * dt;

        m_bodies[i]->setPosition(positions[i]);
        m_bodies[i]->setVelocity(velocities[i]);
        m_bodies[i]->setAcceleration(accs[i]);
        m_bodies[i]->updateTrail();
    }
}

void GravitySimulator::integrateVelocityVerlet(double dt) {
    const size_t n = m_bodies.size();

    std::vector<Vec2> positions(n), velocities(n);
    std::vector<double> masses(n);
    for (size_t i = 0; i < n; ++i) {
        positions[i]  = m_bodies[i]->getPosition();
        velocities[i] = m_bodies[i]->getVelocity();
        masses[i]     = m_bodies[i]->getMass();
    }

    // a(t)
    auto accs0 = computeAllAccelerations(positions, masses);

    // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
    for (size_t i = 0; i < n; ++i) {
        positions[i] += velocities[i] * dt + accs0[i] * (0.5 * dt * dt);
    }

    // a(t+dt) using new positions
    auto accs1 = computeAllAccelerations(positions, masses);

    // v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    for (size_t i = 0; i < n; ++i) {
        velocities[i] += (accs0[i] + accs1[i]) * (0.5 * dt);

        m_bodies[i]->setPosition(positions[i]);
        m_bodies[i]->setVelocity(velocities[i]);
        m_bodies[i]->setAcceleration(accs1[i]);
        m_bodies[i]->updateTrail();
    }
}

void GravitySimulator::integrateRK4(double dt) {
    // Proper coupled RK4: advance ALL bodies together at each sub-step.
    const size_t n = m_bodies.size();

    std::vector<Vec2> pos0(n), vel0(n);
    std::vector<double> masses(n);
    for (size_t i = 0; i < n; ++i) {
        pos0[i]  = m_bodies[i]->getPosition();
        vel0[i]  = m_bodies[i]->getVelocity();
        masses[i] = m_bodies[i]->getMass();
    }

    // k1
    auto a1 = computeAllAccelerations(pos0, masses);
    // k1: dx = vel0, dv = a1

    // k2 positions/velocities at t + dt/2 using k1
    std::vector<Vec2> pos2(n), vel2(n);
    for (size_t i = 0; i < n; ++i) {
        pos2[i] = pos0[i] + vel0[i] * (dt * 0.5);
        vel2[i] = vel0[i] + a1[i] * (dt * 0.5);
    }
    auto a2 = computeAllAccelerations(pos2, masses);

    // k3 positions/velocities at t + dt/2 using k2
    std::vector<Vec2> pos3(n), vel3(n);
    for (size_t i = 0; i < n; ++i) {
        pos3[i] = pos0[i] + vel2[i] * (dt * 0.5);
        vel3[i] = vel0[i] + a2[i] * (dt * 0.5);
    }
    auto a3 = computeAllAccelerations(pos3, masses);

    // k4 positions/velocities at t + dt using k3
    std::vector<Vec2> pos4(n), vel4(n);
    for (size_t i = 0; i < n; ++i) {
        pos4[i] = pos0[i] + vel3[i] * dt;
        vel4[i] = vel0[i] + a3[i] * dt;
    }
    auto a4 = computeAllAccelerations(pos4, masses);

    // Weighted combination
    for (size_t i = 0; i < n; ++i) {
        Vec2 dxdt = (vel0[i] + vel2[i] * 2.0 + vel3[i] * 2.0 + vel4[i]) / 6.0;
        Vec2 dvdt = (a1[i]   + a2[i] * 2.0   + a3[i] * 2.0   + a4[i])   / 6.0;

        m_bodies[i]->setPosition(pos0[i] + dxdt * dt);
        m_bodies[i]->setVelocity(vel0[i] + dvdt * dt);
        m_bodies[i]->setAcceleration(a4[i]);
        m_bodies[i]->updateTrail();
    }
}

// ---------- collisions ----------

void GravitySimulator::handleCollisions() {
    if (!m_detectCollisions) return;

    std::vector<std::shared_ptr<CelestialBody>> result;
    std::vector<bool> absorbed(m_bodies.size(), false);

    for (size_t i = 0; i < m_bodies.size(); ++i) {
        if (absorbed[i]) continue;

        // Accumulate all bodies that collide with i (or with each other in the cluster)
        double totalMass = m_bodies[i]->getMass();
        Vec2 momentum = m_bodies[i]->getVelocity() * totalMass;
        Vec2 weightedPos = m_bodies[i]->getPosition() * totalMass;

        // Mass-weighted colour accumulation
        double colR = m_bodies[i]->getColor().r * totalMass;
        double colG = m_bodies[i]->getColor().g * totalMass;
        double colB = m_bodies[i]->getColor().b * totalMass;

        float maxRadius = m_bodies[i]->getDisplayRadius();
        std::string name = m_bodies[i]->getName();
        bool merged = false;

        for (size_t j = i + 1; j < m_bodies.size(); ++j) {
            if (absorbed[j]) continue;

            Vec2 diff = m_bodies[j]->getPosition() - m_bodies[i]->getPosition();
            double dist = diff.length();

            // Collision threshold: sum of display radii converted to metres
            // We use a heuristic: displayRadius * metersPerPixel isn't available here,
            // so we check if the physical distance is small relative to a configurable
            // softening length. A more robust approach would use physical radii.
            double collisionDist = m_softening * 10.0;
            if (dist < collisionDist) {
                absorbed[j] = true;
                merged = true;

                double mj = m_bodies[j]->getMass();
                momentum   += m_bodies[j]->getVelocity() * mj;
                weightedPos += m_bodies[j]->getPosition() * mj;
                totalMass  += mj;

                colR += m_bodies[j]->getColor().r * mj;
                colG += m_bodies[j]->getColor().g * mj;
                colB += m_bodies[j]->getColor().b * mj;

                maxRadius = std::max(maxRadius, m_bodies[j]->getDisplayRadius());
                name += "+" + m_bodies[j]->getName();
            }
        }

        if (merged) {
            Vec2 newPos = weightedPos / totalMass;
            Vec2 newVel = momentum / totalMass;
            sf::Color newColor(
                static_cast<sf::Uint8>(colR / totalMass),
                static_cast<sf::Uint8>(colG / totalMass),
                static_cast<sf::Uint8>(colB / totalMass)
            );
            // Radius grows with cube root of mass ratio
            float newRadius = maxRadius * static_cast<float>(
                std::cbrt(totalMass / m_bodies[i]->getMass()));

            result.push_back(std::make_shared<CelestialBody>(
                name, totalMass, newPos, newVel, newColor, newRadius));
        } else {
            result.push_back(m_bodies[i]);
        }
    }

    if (result.size() != m_bodies.size()) {
        m_bodies = std::move(result);
    }
}

// ---------- public ----------

void GravitySimulator::update() {
    switch (m_method) {
        case IntegrationMethod::EULER:           integrateEuler(m_timeStep); break;
        case IntegrationMethod::VELOCITY_VERLET: integrateVelocityVerlet(m_timeStep); break;
        case IntegrationMethod::RK4:             integrateRK4(m_timeStep); break;
    }
    handleCollisions();
}

void GravitySimulator::addBody(std::shared_ptr<CelestialBody> body) {
    m_bodies.push_back(std::move(body));
}

bool GravitySimulator::removeBody(const std::string& name) {
    auto it = std::find_if(m_bodies.begin(), m_bodies.end(),
        [&name](const auto& b) { return b->getName() == name; });
    if (it != m_bodies.end()) {
        m_bodies.erase(it);
        return true;
    }
    return false;
}

std::shared_ptr<CelestialBody> GravitySimulator::getBody(const std::string& name) const {
    for (const auto& b : m_bodies) {
        if (b->getName() == name) return b;
    }
    return nullptr;
}

void GravitySimulator::clear() {
    m_bodies.clear();
}