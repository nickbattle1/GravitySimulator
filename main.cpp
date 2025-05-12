#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <SFML/Graphics.hpp>

// Universal gravitational constant (scaled for simulation)
const double G = 6.67430e-11;
const double SCALE_FACTOR = 1e9; // Scaling factor for visualization

// Available integration methods
enum class IntegrationMethod {
    EULER,
    VELOCITY_VERLET,
    RK4
};

// Class to represent a celestial body
class CelestialBody {
private:
    std::string name;
    double mass;
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f acceleration;
    sf::Color color;
    float radius;
    std::vector<sf::Vector2f> trail;
    size_t maxTrailLength = 100;
    bool drawTrail = true;

public:
    CelestialBody(const std::string& bodyName, double bodyMass, 
                 const sf::Vector2f& initialPos, const sf::Vector2f& initialVel,
                 const sf::Color& bodyColor, float bodyRadius)
        : name(bodyName), mass(bodyMass), position(initialPos), velocity(initialVel),
          acceleration(sf::Vector2f(0.0f, 0.0f)), color(bodyColor), radius(bodyRadius) {}

    // Getters and setters
    const std::string& getName() const { return name; }
    double getMass() const { return mass; }
    const sf::Vector2f& getPosition() const { return position; }
    const sf::Vector2f& getVelocity() const { return velocity; }
    const sf::Vector2f& getAcceleration() const { return acceleration; }
    const sf::Color& getColor() const { return color; }
    float getRadius() const { return radius; }
    bool getDrawTrail() const { return drawTrail; }
    void setDrawTrail(bool draw) { drawTrail = draw; }
    const std::vector<sf::Vector2f>& getTrail() const { return trail; }
    
    void setPosition(const sf::Vector2f& pos) { position = pos; }
    void setVelocity(const sf::Vector2f& vel) { velocity = vel; }
    void setAcceleration(const sf::Vector2f& acc) { acceleration = acc; }

    // Update trail with current position
    void updateTrail() {
        if (drawTrail) {
            trail.push_back(position);
            if (trail.size() > maxTrailLength) {
                trail.erase(trail.begin());
            }
        }
    }

    // Clear the trail
    void clearTrail() {
        trail.clear();
    }

    // Set maximum trail length
    void setMaxTrailLength(size_t length) {
        maxTrailLength = length;
        if (trail.size() > maxTrailLength) {
            trail.resize(maxTrailLength);
        }
    }

    // Draw the celestial body and its trail
    void draw(sf::RenderWindow& window, float zoomLevel = 1.0f) const {
        // Draw the trail
        if (drawTrail && trail.size() > 1) {
            sf::VertexArray lines(sf::LineStrip, trail.size());
            for (size_t i = 0; i < trail.size(); ++i) {
                lines[i].position = trail[i] * zoomLevel + sf::Vector2f(window.getSize().x / 2, window.getSize().y / 2);
                
                // Fade the trail from the body color to transparent
                float alpha = static_cast<float>(i) / trail.size() * 255;
                lines[i].color = sf::Color(color.r, color.g, color.b, static_cast<sf::Uint8>(alpha));
            }
            window.draw(lines);
        }

        // Draw the body
        sf::CircleShape shape(radius * zoomLevel);
        shape.setFillColor(color);
        shape.setOrigin(radius * zoomLevel, radius * zoomLevel);
        shape.setPosition(position * zoomLevel + sf::Vector2f(window.getSize().x / 2, window.getSize().y / 2));
        window.draw(shape);
    }
};

// Class to manage the physics simulation
class GravitySimulator {
private:
    std::vector<std::shared_ptr<CelestialBody>> bodies;
    IntegrationMethod integrationMethod = IntegrationMethod::VELOCITY_VERLET;
    double timeStep = 3600.0; // Default time step (in seconds): 1 hour
    bool detectCollisions = false;
    
    // Calculate gravitational force between two bodies
    sf::Vector2f calculateGravitationalForce(const CelestialBody& body1, const CelestialBody& body2) const {
        sf::Vector2f direction = body2.getPosition() - body1.getPosition();
        float distanceSquared = direction.x * direction.x + direction.y * direction.y;
        
        // Prevent division by zero or very small values
        if (distanceSquared < 1e-10) {
            return sf::Vector2f(0.0f, 0.0f);
        }
        
        float distance = std::sqrt(distanceSquared);
        
        // Normalize direction vector
        direction.x /= distance;
        direction.y /= distance;
        
        // Calculate force magnitude using Newton's law of universal gravitation
        float forceMagnitude = G * body1.getMass() * body2.getMass() / distanceSquared;
        
        // Force vector (direction * magnitude)
        return direction * forceMagnitude;
    }
    
    // Calculate the total acceleration for a body
    sf::Vector2f calculateAcceleration(const CelestialBody& body) const {
        sf::Vector2f totalForce(0.0f, 0.0f);
        
        // Sum up forces from all other bodies
        for (const auto& otherBody : bodies) {
            if (otherBody.get() != &body) {
                sf::Vector2f force = calculateGravitationalForce(body, *otherBody);
                totalForce += force;
            }
        }
        
        // a = F/m (Newton's second law)
        return totalForce / static_cast<float>(body.getMass());
    }
    
    // Handle collisions between bodies
    void handleCollisions() {
        if (!detectCollisions) return;
        
        std::vector<std::shared_ptr<CelestialBody>> newBodies;
        std::vector<bool> merged(bodies.size(), false);
        
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (merged[i]) continue;
            
            auto& body1 = bodies[i];
            sf::Vector2f totalMomentum = body1->getMass() * body1->getVelocity();
            double totalMass = body1->getMass();
            sf::Vector2f weightedPosition = body1->getMass() * body1->getPosition();
            std::string newName = body1->getName();
            float maxRadius = body1->getRadius();
            sf::Color avgColor = body1->getColor();
            
            bool collided = false;
            
            for (size_t j = i + 1; j < bodies.size(); ++j) {
                if (merged[j]) continue;
                
                auto& body2 = bodies[j];
                sf::Vector2f direction = body2->getPosition() - body1->getPosition();
                float distance = std::sqrt(direction.x * direction.x + direction.y * direction.y);
                
                if (distance < (body1->getRadius() + body2->getRadius())) {
                    // Collision detected
                    collided = true;
                    merged[j] = true;
                    
                    // Add to totals for the new merged body
                    totalMomentum += body2->getMass() * body2->getVelocity();
                    weightedPosition += body2->getMass() * body2->getPosition();
                    totalMass += body2->getMass();
                    newName += "+" + body2->getName();
                    maxRadius = std::max(maxRadius, body2->getRadius());
                    
                    // Mix colors
                    avgColor.r = static_cast<sf::Uint8>((avgColor.r + body2->getColor().r) / 2);
                    avgColor.g = static_cast<sf::Uint8>((avgColor.g + body2->getColor().g) / 2);
                    avgColor.b = static_cast<sf::Uint8>((avgColor.b + body2->getColor().b) / 2);
                }
            }
            
            if (collided) {
                // Create a new merged body
                sf::Vector2f newPosition = weightedPosition / static_cast<float>(totalMass);
                sf::Vector2f newVelocity = totalMomentum / static_cast<float>(totalMass);
                
                // Calculate new radius based on mass (assuming constant density)
                float newRadius = std::cbrt(totalMass / body1->getMass()) * maxRadius;
                
                auto newBody = std::make_shared<CelestialBody>(
                    newName, totalMass, newPosition, newVelocity, avgColor, newRadius
                );
                newBodies.push_back(newBody);
            } else {
                newBodies.push_back(body1);
            }
        }
        
        // Replace bodies with newBodies if any collisions occurred
        if (bodies.size() != newBodies.size()) {
            bodies = newBodies;
        }
    }

public:
    // Integration methods implementation
    void updateEuler(double dt) {
        // Calculate accelerations for all bodies
        for (auto& body : bodies) {
            body->setAcceleration(calculateAcceleration(*body));
        }
        
        // Update positions and velocities
        for (auto& body : bodies) {
            sf::Vector2f position = body->getPosition();
            sf::Vector2f velocity = body->getVelocity();
            sf::Vector2f acceleration = body->getAcceleration();
            
            // Update position: x = x0 + v0*dt
            position += velocity * static_cast<float>(dt);
            
            // Update velocity: v = v0 + a*dt
            velocity += acceleration * static_cast<float>(dt);
            
            body->setPosition(position);
            body->setVelocity(velocity);
            body->updateTrail();
        }
        
        handleCollisions();
    }
    
    void updateVelocityVerlet(double dt) {
        // Store initial accelerations
        std::vector<sf::Vector2f> initialAccelerations;
        for (auto& body : bodies) {
            sf::Vector2f acceleration = calculateAcceleration(*body);
            body->setAcceleration(acceleration);
            initialAccelerations.push_back(acceleration);
        }
        
        // Update positions and half-step velocities
        for (size_t i = 0; i < bodies.size(); ++i) {
            auto& body = bodies[i];
            sf::Vector2f position = body->getPosition();
            sf::Vector2f velocity = body->getVelocity();
            sf::Vector2f acceleration = initialAccelerations[i];
            
            // Update position: x = x0 + v0*dt + 0.5*a0*dt^2
            position += velocity * static_cast<float>(dt) + 
                        0.5f * acceleration * static_cast<float>(dt * dt);
            
            // Update velocity half-step: v_half = v0 + 0.5*a0*dt
            velocity += 0.5f * acceleration * static_cast<float>(dt);
            
            body->setPosition(position);
            body->setVelocity(velocity);
        }
        
        // Calculate new accelerations
        for (auto& body : bodies) {
            body->setAcceleration(calculateAcceleration(*body));
        }
        
        // Complete velocity updates
        for (size_t i = 0; i < bodies.size(); ++i) {
            auto& body = bodies[i];
            sf::Vector2f velocity = body->getVelocity();
            sf::Vector2f newAcceleration = body->getAcceleration();
            
            // Update velocity: v = v_half + 0.5*a_new*dt
            velocity += 0.5f * newAcceleration * static_cast<float>(dt);
            
            body->setVelocity(velocity);
            body->updateTrail();
        }
        
        handleCollisions();
    }
    
    void updateRK4(double dt) {
        struct Derivative {
            sf::Vector2f velocity;
            sf::Vector2f acceleration;
        };
        
        auto evaluate = [this](const CelestialBody& body, double t, const Derivative& d) -> Derivative {
            sf::Vector2f position = body.getPosition() + d.velocity * static_cast<float>(t);
            sf::Vector2f velocity = body.getVelocity() + d.acceleration * static_cast<float>(t);
            
            Derivative output;
            output.velocity = velocity;
            
            // Create a temporary body to calculate acceleration at the new position
            CelestialBody tempBody(body.getName(), body.getMass(), position, velocity, 
                                  body.getColor(), body.getRadius());
            output.acceleration = calculateAcceleration(tempBody);
            
            return output;
        };
        
        for (auto& body : bodies) {
            Derivative a, b, c, d;
            
            // Initial derivative
            a.velocity = body->getVelocity();
            a.acceleration = calculateAcceleration(*body);
            
            // Second derivative
            b = evaluate(*body, dt * 0.5, a);
            
            // Third derivative
            c = evaluate(*body, dt * 0.5, b);
            
            // Fourth derivative
            d = evaluate(*body, dt, c);
            
            // Weighted sum for position and velocity updates
            sf::Vector2f position = body->getPosition();
            sf::Vector2f velocity = body->getVelocity();
            
            sf::Vector2f dxdt = (a.velocity + 2.0f * (b.velocity + c.velocity) + d.velocity) / 6.0f;
            sf::Vector2f dvdt = (a.acceleration + 2.0f * (b.acceleration + c.acceleration) + d.acceleration) / 6.0f;
            
            position += dxdt * static_cast<float>(dt);
            velocity += dvdt * static_cast<float>(dt);
            
            body->setPosition(position);
            body->setVelocity(velocity);
            body->updateTrail();
        }
        
        handleCollisions();
    }
    
    // Main update method that calls the appropriate integration method
    void update() {
        switch (integrationMethod) {
            case IntegrationMethod::EULER:
                updateEuler(timeStep);
                break;
            case IntegrationMethod::VELOCITY_VERLET:
                updateVelocityVerlet(timeStep);
                break;
            case IntegrationMethod::RK4:
                updateRK4(timeStep);
                break;
        }
    }
    
    // Add a celestial body to the simulation
    void addBody(const std::shared_ptr<CelestialBody>& body) {
        bodies.push_back(body);
    }
    
    // Remove a celestial body by name
    bool removeBody(const std::string& name) {
        for (auto it = bodies.begin(); it != bodies.end(); ++it) {
            if ((*it)->getName() == name) {
                bodies.erase(it);
                return true;
            }
        }
        return false;
    }
    
    // Get a body by name
    std::shared_ptr<CelestialBody> getBody(const std::string& name) const {
        for (const auto& body : bodies) {
            if (body->getName() == name) {
                return body;
            }
        }
        return nullptr;
    }
    
    // Get all bodies
    const std::vector<std::shared_ptr<CelestialBody>>& getBodies() const {
        return bodies;
    }
    
    // Set integration method
    void setIntegrationMethod(IntegrationMethod method) {
        integrationMethod = method;
    }
    
    // Get current integration method
    IntegrationMethod getIntegrationMethod() const {
        return integrationMethod;
    }
    
    // Set time step
    void setTimeStep(double step) {
        timeStep = step;
    }
    
    // Get current time step
    double getTimeStep() const {
        return timeStep;
    }
    
    // Set collision detection
    void setCollisionDetection(bool detect) {
        detectCollisions = detect;
    }
    
    // Get collision detection status
    bool getCollisionDetection() const {
        return detectCollisions;
    }
    
    // Clear all bodies
    void clear() {
        bodies.clear();
    }
};

// Command-line interface for the simulator
class SimulatorCLI {
private:
    GravitySimulator simulator;
    bool running = false;
    float zoomLevel = 1.0f;
    sf::Vector2f viewCenter;
    sf::Clock clock;
    float simulationSpeed = 1.0f;

    // Initialize the solar system
    void initializeSolarSystem() {
        simulator.clear();
        
        // Sun (position at origin)
        auto sun = std::make_shared<CelestialBody>(
            "Sun", 1.989e30, sf::Vector2f(0, 0), sf::Vector2f(0, 0),
            sf::Color::Yellow, 20.0f
        );
        simulator.addBody(sun);
        
        // Mercury
        auto mercury = std::make_shared<CelestialBody>(
            "Mercury", 3.3011e23, sf::Vector2f(57.9e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 47.4e3), sf::Color(169, 169, 169), 3.5f
        );
        simulator.addBody(mercury);
        
        // Venus
        auto venus = std::make_shared<CelestialBody>(
            "Venus", 4.8675e24, sf::Vector2f(108.2e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 35.0e3), sf::Color(255, 198, 73), 8.0f
        );
        simulator.addBody(venus);
        
        // Earth
        auto earth = std::make_shared<CelestialBody>(
            "Earth", 5.972e24, sf::Vector2f(149.6e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 29.8e3), sf::Color(0, 0, 255), 9.0f
        );
        simulator.addBody(earth);
        
        // Moon (orbiting Earth)
        auto moon = std::make_shared<CelestialBody>(
            "Moon", 7.342e22, sf::Vector2f(149.6e9 + 0.384e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 29.8e3 + 1.022e3), sf::Color(200, 200, 200), 2.5f
        );
        simulator.addBody(moon);
        
        // Mars
        auto mars = std::make_shared<CelestialBody>(
            "Mars", 6.4171e23, sf::Vector2f(227.9e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 24.1e3), sf::Color(255, 0, 0), 4.8f
        );
        simulator.addBody(mars);
        
        // Jupiter
        auto jupiter = std::make_shared<CelestialBody>(
            "Jupiter", 1.8982e27, sf::Vector2f(778.5e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 13.1e3), sf::Color(255, 140, 0), 15.0f
        );
        simulator.addBody(jupiter);
    }

    // Initialize a simple two-body system
    void initializeTwoBodySystem() {
        simulator.clear();
        
        // Star
        auto star = std::make_shared<CelestialBody>(
            "Star", 1.0e30, sf::Vector2f(0, 0), sf::Vector2f(0, 0),
            sf::Color::Yellow, 15.0f
        );
        simulator.addBody(star);
        
        // Planet
        auto planet = std::make_shared<CelestialBody>(
            "Planet", 1.0e25, sf::Vector2f(100e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 25.0e3), sf::Color::Blue, 5.0f
        );
        simulator.addBody(planet);
    }

    // Create a binary star system
    void initializeBinaryStarSystem() {
        simulator.clear();
        
        // Star 1
        auto star1 = std::make_shared<CelestialBody>(
            "Star1", 1.5e30, sf::Vector2f(-50e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, -10.0e3), sf::Color(255, 200, 100), 12.0f
        );
        simulator.addBody(star1);
        
        // Star 2
        auto star2 = std::make_shared<CelestialBody>(
            "Star2", 1.2e30, sf::Vector2f(50e9, 0) / SCALE_FACTOR,
            sf::Vector2f(0, 12.5e3), sf::Color(100, 200, 255), 10.0f
        );
        simulator.addBody(star2);
        
        // Planet orbiting the binary system
        auto planet = std::make_shared<CelestialBody>(
            "Planet", 2.0e24, sf::Vector2f(0, 200e9) / SCALE_FACTOR,
            sf::Vector2f(18.0e3, 0), sf::Color::Green, 4.0f
        );
        simulator.addBody(planet);
    }

    // Process user input from command line
    void processCommand(const std::string& command) {
        std::string cmd = command;
        std::string arg;
        
        // Split command and arguments
        size_t spacePos = command.find(' ');
        if (spacePos != std::string::npos) {
            cmd = command.substr(0, spacePos);
            arg = command.substr(spacePos + 1);
        }
        
        if (cmd == "help") {
            std::cout << "Available commands:\n"
                      << "help - Show this help message\n"
                      << "add <name> <mass> <posX> <posY> <velX> <velY> <colorR> <colorG> <colorB> <radius> - Add new body\n"
                      << "remove <name> - Remove a body by name\n"
                      << "list - List all bodies\n"
                      << "clear - Remove all bodies\n"
                      << "solar - Initialize solar system\n"
                      << "binary - Initialize binary star system\n"
                      << "two-body - Initialize simple two-body system\n"
                      << "method <euler|verlet|rk4> - Set integration method\n"
                      << "timestep <seconds> - Set simulation time step\n"
                      << "speed <factor> - Set simulation speed multiplier\n"
                      << "collisions <on|off> - Enable/disable collision detection\n"
                      << "trails <on|off> - Enable/disable trails for all bodies\n"
                      << "traillength <number> - Set maximum trail length\n"
                      << "quit - Exit the simulation\n";
        } else if (cmd == "add") {
            try {
                std::istringstream iss(arg);
                std::string name;
                double mass, posX, posY, velX, velY;
                int colorR, colorG, colorB;
                float radius;
                
                iss >> name >> mass >> posX >> posY >> velX >> velY >> colorR >> colorG >> colorB >> radius;
                
                auto body = std::make_shared<CelestialBody>(
                    name, mass, sf::Vector2f(posX, posY), sf::Vector2f(velX, velY),
                    sf::Color(colorR, colorG, colorB), radius
                );
                simulator.addBody(body);
                std::cout << "Added body: " << name << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error adding body. Format: add <name> <mass> <posX> <posY> <velX> <velY> <colorR> <colorG> <colorB> <radius>" << std::endl;
            }
        } else if (cmd == "remove") {
            if (simulator.removeBody(arg)) {
                std::cout << "Removed body: " << arg << std::endl;
            } else {
                std::cout << "Body not found: " << arg << std::endl;
            }
        } else if (cmd == "list") {
            const auto& bodies = simulator.getBodies();
            std::cout << "Bodies in simulation (" << bodies.size() << "):" << std::endl;
            for (const auto& body : bodies) {
                const sf::Vector2f& pos = body->getPosition();
                const sf::Vector2f& vel = body->getVelocity();
                std::cout << "- " << body->getName() << ": "
                          << "Mass=" << body->getMass() << ", "
                          << "Pos=(" << pos.x << ", " << pos.y << "), "
                          << "Vel=(" << vel.x << ", " << vel.y << "), "
                          << "Radius=" << body->getRadius() << std::endl;
            }
        } else if (cmd == "clear") {
            simulator.clear();
            std::cout << "Cleared all bodies" << std::endl;
        } else if (cmd == "solar") {
            initializeSolarSystem();
            std::cout << "Initialized solar system" << std::endl;
        } else if (cmd == "binary") {
            initializeBinaryStarSystem();
            std::cout << "Initialized binary star system" << std::endl;
        } else if (cmd == "two-body") {
            initializeTwoBodySystem();
            std::cout << "Initialized two-body system" << std::endl;
        } else if (cmd == "method") {
            if (arg == "euler") {
                simulator.setIntegrationMethod(IntegrationMethod::EULER);
                std::cout << "Set integration method to Euler" << std::endl;
            } else if (arg == "verlet") {
                simulator.setIntegrationMethod(IntegrationMethod::VELOCITY_VERLET);
                std::cout << "Set integration method to Velocity Verlet" << std::endl;
            } else if (arg == "rk4") {
                simulator.setIntegrationMethod(IntegrationMethod::RK4);
                std::cout << "Set integration method to Runge-Kutta 4" << std::endl;
            } else {
                std::cout << "Unknown integration method. Use 'euler', 'verlet', or 'rk4'" << std::endl;
            }
        } else if (cmd == "timestep") {
            try {
                double step = std::stod(arg);
                simulator.setTimeStep(step);
                std::cout << "Set time step to " << step << " seconds" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Invalid time step. Format: timestep <seconds>" << std::endl;
            }
        } else if (cmd == "speed") {
            try {
                float speed = std::stof(arg);
                if (speed > 0.0f) {
                    simulationSpeed = speed;
                    std::cout << "Set simulation speed to " << speed << "x" << std::endl;
                } else {
                    std::cout << "Speed must be positive" << std::endl;
                }
            } catch (const std::exception& e) {
                std::cout << "Invalid speed. Format: speed <factor>" << std::endl;
            }
        } else if (cmd == "collisions") {
            if (arg == "on") {
                simulator.setCollisionDetection(true);
                std::cout << "Enabled collision detection" << std::endl;
            } else if (arg == "off") {
                simulator.setCollisionDetection(false);
                std::cout << "Disabled collision detection" << std::endl;
            } else {
                std::cout << "Invalid option. Use 'on' or 'off'" << std::endl;
            }
        } else if (cmd == "trails") {
            bool enableTrails = (arg == "on");
            for (const auto& body : simulator.getBodies()) {
                body->setDrawTrail(enableTrails);
            }
            std::cout << (enableTrails ? "Enabled" : "Disabled") << " trails for all bodies" << std::endl;
        } else if (cmd == "traillength") {
            try {
                size_t length = std::stoul(arg);
                for (const auto& body : simulator.getBodies()) {
                    body->setMaxTrailLength(length);
                }
                std::cout << "Set trail length to " << length << " points" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Invalid trail length. Format: traillength <number>" << std::endl;
            }
        } else if (cmd == "quit") {
            running = false;
        } else if (!cmd.empty()) {
            std::cout << "Unknown command: " << cmd << ". Type 'help' for available commands." << std::endl;
        }
    }

public:
    // Run the simulator with visualization
    void run() {
        running = true;
        
        // Initialize solar system as default
        initializeSolarSystem();
        
        // Create SFML window
        sf::RenderWindow window(sf::VideoMode(1200, 800), "Gravity Simulator");
        window.setFramerateLimit(60);
        
        // Font for text display
        sf::Font font;
        if (!font.loadFromFile("arial.ttf")) {
            // Try to load a system font if arial.ttf is not available
            std::cout << "Warning: Could not load arial.ttf. Using default font." << std::endl;
        }
        
        sf::Text infoText;
        infoText.setFont(font);
        infoText.setCharacterSize(14);
        infoText.setFillColor(sf::Color::White);
        infoText.setPosition(10, 10);
        
        sf::Text helpText;
        helpText.setFont(font);
        helpText.setCharacterSize(14);
        helpText.setFillColor(sf::Color::White);
        helpText.setPosition(10, window.getSize().y - 100);
        helpText.setString(
            "Controls:\n"
            "Arrows: Pan view | +/-: Zoom in/out | Space: Pause/Resume | Enter: Command mode"
        );
        
        bool paused = false;
        bool commandMode = false;
        std::string currentCommand;
        sf::Text commandText;
        commandText.setFont(font);
        commandText.setCharacterSize(16);
        commandText.setFillColor(sf::Color::Yellow);
        commandText.setPosition(10, window.getSize().y - 40);
        
        clock.restart();
        
        while (running && window.isOpen()) {
            float deltaTime = clock.restart().asSeconds();
            
            // Handle events
            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    window.close();
                    running = false;
                } else if (event.type == sf::Event::KeyPressed) {
                    if (commandMode) {
                        if (event.key.code == sf::Keyboard::Return) {
                            // Execute command and exit command mode
                            processCommand(currentCommand);
                            currentCommand.clear();
                            commandMode = false;
                        } else if (event.key.code == sf::Keyboard::Escape) {
                            // Cancel command mode
                            currentCommand.clear();
                            commandMode = false;
                        } else if (event.key.code == sf::Keyboard::BackSpace) {
                            // Delete last character
                            if (!currentCommand.empty()) {
                                currentCommand.pop_back();
                            }
                        }
                    } else {
                        // Regular key controls
                        switch (event.key.code) {
                            case sf::Keyboard::Space:
                                paused = !paused;
                                break;
                            case sf::Keyboard::Return:
                                commandMode = true;
                                break;
                            case sf::Keyboard::Add:
                            case sf::Keyboard::Equal:
                                zoomLevel *= 1.1f;
                                break;
                            case sf::Keyboard::Subtract:
                            case sf::Keyboard::Dash:
                                zoomLevel /= 1.1f;
                                break;
                            case sf::Keyboard::Left:
                                viewCenter.x -= 10.0f / zoomLevel;
                                break;
                            case sf::Keyboard::Right:
                                viewCenter.x += 10.0f / zoomLevel;
                                break;
                            case sf::Keyboard::Up:
                                viewCenter.y -= 10.0f / zoomLevel;
                                break;
                            case sf::Keyboard::Down:
                                viewCenter.y += 10.0f / zoomLevel;
                                break;
                            case sf::Keyboard::R:
                                // Reset view
                                viewCenter = sf::Vector2f(0, 0);
                                zoomLevel = 1.0f;
                                break;
                            case sf::Keyboard::Escape:
                                window.close();
                                running = false;
                                break;
                            default:
                                break;
                        }
                    }
                } else if (event.type == sf::Event::TextEntered && commandMode) {
                    // Handle text input for command mode
                    if (event.text.unicode >= 32 && event.text.unicode < 128) {
                        currentCommand += static_cast<char>(event.text.unicode);
                    }
                }
            }
            
            // Update simulation if not paused
            if (!paused) {
                // Run simulation at current speed
                float simTime = deltaTime * simulationSpeed;
                double timeStep = simulator.getTimeStep();
                
                // Run multiple steps if necessary
                double accumulatedTime = simTime;
                while (accumulatedTime >= timeStep / 1000.0) {
                    simulator.update();
                    accumulatedTime -= timeStep / 1000.0;
                }
            }
            
            // Clear the window
            window.clear(sf::Color(10, 10, 50));
            
            // Draw all bodies
            for (const auto& body : simulator.getBodies()) {
                body->draw(window, zoomLevel);
            }
            
            // Update and draw info text
            std::string methodName;
            switch (simulator.getIntegrationMethod()) {
                case IntegrationMethod::EULER:
                    methodName = "Euler";
                    break;
                case IntegrationMethod::VELOCITY_VERLET:
                    methodName = "Velocity Verlet";
                    break;
                case IntegrationMethod::RK4:
                    methodName = "Runge-Kutta 4";
                    break;
            }
            
            infoText.setString(
                "Bodies: " + std::to_string(simulator.getBodies().size()) + "\n" +
                "Integration: " + methodName + "\n" +
                "Time Step: " + std::to_string(simulator.getTimeStep()) + " s\n" +
                "Speed: " + std::to_string(simulationSpeed) + "x\n" +
                "Collisions: " + std::string(simulator.getCollisionDetection() ? "On" : "Off") + "\n" +
                "Zoom: " + std::to_string(zoomLevel) + "x\n" +
                "Status: " + (paused ? "Paused" : "Running")
            );
            window.draw(infoText);
            
            // Draw help text
            window.draw(helpText);
            
            // Draw command text if in command mode
            if (commandMode) {
                commandText.setString("> " + currentCommand + "_");
                window.draw(commandText);
            }
            
            // Display the window
            window.display();
        }
    }
};

int main() {
    SimulatorCLI simulator;
    simulator.run();
    return 0;
}