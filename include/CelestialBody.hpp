#pragma once

#include "Vec2.hpp"
#include <string>
#include <deque>
#include <SFML/Graphics.hpp>

class CelestialBody {
public:
    CelestialBody(const std::string& name, double mass,
                  const Vec2& position, const Vec2& velocity,
                  const sf::Color& color, float displayRadius);

    // Getters
    const std::string& getName() const { return m_name; }
    double getMass() const { return m_mass; }
    const Vec2& getPosition() const { return m_position; }
    const Vec2& getVelocity() const { return m_velocity; }
    const Vec2& getAcceleration() const { return m_acceleration; }
    const sf::Color& getColor() const { return m_color; }
    float getDisplayRadius() const { return m_displayRadius; }
    float getScaledRadius(double metersPerPixel) const;
    bool getDrawTrail() const { return m_drawTrail; }

    // Setters
    void setPosition(const Vec2& pos) { m_position = pos; }
    void setVelocity(const Vec2& vel) { m_velocity = vel; }
    void setAcceleration(const Vec2& acc) { m_acceleration = acc; }
    void setDrawTrail(bool draw) { m_drawTrail = draw; }
    void setMaxTrailLength(size_t length);

    // Trail management
    void updateTrail();
    void clearTrail();

    // Render the body and trail.
    // worldToScreen converts physics coordinates to screen pixels.
    void draw(sf::RenderWindow& window,
              const Vec2& viewCenter, double metersPerPixel) const;

private:
    std::string m_name;
    double m_mass;
    Vec2 m_position;      // metres (real physics units)
    Vec2 m_velocity;      // m/s
    Vec2 m_acceleration;  // m/s²
    sf::Color m_color;
    float m_displayRadius; // pixels (display only)

    std::deque<Vec2> m_trail;
    size_t m_maxTrailLength = 200;
    bool m_drawTrail = true;

    sf::Vector2f worldToScreen(const Vec2& worldPos,
                               const Vec2& viewCenter,
                               double metersPerPixel,
                               const sf::RenderWindow& window) const;
};