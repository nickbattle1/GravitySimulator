#include "CelestialBody.hpp"
#include <algorithm>

CelestialBody::CelestialBody(const std::string& name, double mass,
                             const Vec2& position, const Vec2& velocity,
                             const sf::Color& color, float displayRadius)
    : m_name(name), m_mass(mass), m_position(position), m_velocity(velocity),
      m_acceleration({0.0, 0.0}), m_color(color), m_displayRadius(displayRadius) {}

float CelestialBody::getScaledRadius(double metersPerPixel) const {
    static constexpr double kReferenceScale = 2.0e9;
    float scale = static_cast<float>(kReferenceScale / metersPerPixel);
    float scaled = m_displayRadius * scale;
    return std::clamp(scaled, 1.5f, m_displayRadius * 5.0f);
}

void CelestialBody::setMaxTrailLength(size_t length) {
    m_maxTrailLength = length;
    while (m_trail.size() > m_maxTrailLength) {
        m_trail.pop_front();
    }
}

void CelestialBody::updateTrail() {
    if (!m_drawTrail) return;
    m_trail.push_back(m_position);
    while (m_trail.size() > m_maxTrailLength) {
        m_trail.pop_front(); // O(1) with deque
    }
}

void CelestialBody::clearTrail() {
    m_trail.clear();
}

sf::Vector2f CelestialBody::worldToScreen(const Vec2& worldPos,
                                           const Vec2& viewCenter,
                                           double metersPerPixel,
                                           const sf::RenderWindow& window) const {
    Vec2 offset = (worldPos - viewCenter) / metersPerPixel;
    float screenX = static_cast<float>(offset.x) + window.getSize().x * 0.5f;
    float screenY = static_cast<float>(offset.y) + window.getSize().y * 0.5f;
    return {screenX, screenY};
}

void CelestialBody::draw(sf::RenderWindow& window,
                          const Vec2& viewCenter, double metersPerPixel) const {
    // Draw trail
    if (m_drawTrail && m_trail.size() > 1) {
        sf::VertexArray lines(sf::LineStrip, m_trail.size());
        for (size_t i = 0; i < m_trail.size(); ++i) {
            lines[i].position = worldToScreen(m_trail[i], viewCenter, metersPerPixel, window);
            float alpha = static_cast<float>(i) / static_cast<float>(m_trail.size()) * 255.0f;
            lines[i].color = sf::Color(m_color.r, m_color.g, m_color.b,
                                       static_cast<sf::Uint8>(alpha));
        }
        window.draw(lines);
    }

    // Draw body
    sf::Vector2f screenPos = worldToScreen(m_position, viewCenter, metersPerPixel, window);
    float radius = getScaledRadius(metersPerPixel);
    sf::CircleShape shape(radius);
    shape.setFillColor(m_color);
    shape.setOrigin(radius, radius);
    shape.setPosition(screenPos);
    window.draw(shape);

    // Draw name label
    // (Font rendering handled externally for simplicity)
}