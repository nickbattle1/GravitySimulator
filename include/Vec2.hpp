#pragma once

#include <cmath>
#include <SFML/System/Vector2.hpp>

// Double-precision 2D vector for physics calculations.
// Avoids float precision loss when dealing with astronomical scales.
struct Vec2 {
    double x = 0.0;
    double y = 0.0;

    Vec2() = default;
    Vec2(double x, double y) : x(x), y(y) {}

    Vec2 operator+(const Vec2& other) const { return {x + other.x, y + other.y}; }
    Vec2 operator-(const Vec2& other) const { return {x - other.x, y - other.y}; }
    Vec2 operator*(double scalar) const { return {x * scalar, y * scalar}; }
    Vec2 operator/(double scalar) const { return {x / scalar, y / scalar}; }

    Vec2& operator+=(const Vec2& other) { x += other.x; y += other.y; return *this; }
    Vec2& operator-=(const Vec2& other) { x -= other.x; y -= other.y; return *this; }
    Vec2& operator*=(double scalar) { x *= scalar; y *= scalar; return *this; }

    double lengthSquared() const { return x * x + y * y; }
    double length() const { return std::sqrt(lengthSquared()); }

    Vec2 normalized() const {
        double len = length();
        if (len < 1e-30) return {0.0, 0.0};
        return {x / len, y / len};
    }

    // Convert to SFML float vector for rendering
    sf::Vector2f toSFML() const {
        return sf::Vector2f(static_cast<float>(x), static_cast<float>(y));
    }

    static Vec2 fromSFML(const sf::Vector2f& v) {
        return {static_cast<double>(v.x), static_cast<double>(v.y)};
    }
};

inline Vec2 operator*(double scalar, const Vec2& v) { return v * scalar; }