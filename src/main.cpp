#include "GravitySimulator.hpp"
#include "Scenarios.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

static std::string formatSI(double value) {
    const char* suffixes[] = {"", "k", "M", "G", "T", "P", "E"};
    int idx = 0;
    double v = std::abs(value);
    while (v >= 1000.0 && idx < 6) { v /= 1000.0; ++idx; }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << (value < 0 ? -v : v) << suffixes[idx];
    return oss.str();
}

static bool loadFont(sf::Font& font) {
    const char* paths[] = {
        "arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/noto/NotoSans-Regular.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/SFNSText.ttf",
        "C:\\Windows\\Fonts\\arial.ttf",
    };
    for (const auto& path : paths) {
        if (font.loadFromFile(path)) return true;
    }
    return false;
}

class Application {
public:
    void run() {
        sf::RenderWindow window(sf::VideoMode(1400, 900), "Gravity Simulator");
        window.setFramerateLimit(60);

        sf::Font font;
        bool fontLoaded = loadFont(font);
        if (!fontLoaded) {
            std::cerr << "Warning: could not load any font. Text will not display.\n";
        }

        Scenarios::loadSolarSystem(m_sim);

        double metersPerPixel = 2.0e9;
        Vec2 viewCenter{0.0, 0.0};

        sf::Clock clock;
        double accumulator = 0.0;

        bool paused = false;
        bool commandMode = false;
        std::string currentCommand;
        float simulationSpeed = 86400.0f;

        while (window.isOpen()) {
            float frameDelta = clock.restart().asSeconds();
            if (frameDelta > 0.1f) frameDelta = 0.1f;

            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    window.close();
                }
                else if (event.type == sf::Event::MouseWheelScrolled) {
                    double mouseX = static_cast<double>(event.mouseWheelScroll.x);
                    double mouseY = static_cast<double>(event.mouseWheelScroll.y);
                    double halfW = window.getSize().x * 0.5;
                    double halfH = window.getSize().y * 0.5;

                    double worldX = viewCenter.x + (mouseX - halfW) * metersPerPixel;
                    double worldY = viewCenter.y + (mouseY - halfH) * metersPerPixel;

                    if (event.mouseWheelScroll.delta > 0)
                        metersPerPixel *= 0.9;
                    else
                        metersPerPixel *= 1.1;

                    viewCenter.x = worldX - (mouseX - halfW) * metersPerPixel;
                    viewCenter.y = worldY - (mouseY - halfH) * metersPerPixel;
                }
                else if (event.type == sf::Event::KeyPressed) {
                    if (commandMode) {
                        if (event.key.code == sf::Keyboard::Return) {
                            processCommand(currentCommand, simulationSpeed,
                                           metersPerPixel, viewCenter);
                            currentCommand.clear();
                            commandMode = false;
                        } else if (event.key.code == sf::Keyboard::Escape) {
                            currentCommand.clear();
                            commandMode = false;
                        } else if (event.key.code == sf::Keyboard::BackSpace &&
                                   !currentCommand.empty()) {
                            currentCommand.pop_back();
                        }
                    } else {
                        handleKeyPress(event.key.code, paused, commandMode,
                                       metersPerPixel, viewCenter,
                                       simulationSpeed, window);
                    }
                }
                else if (event.type == sf::Event::TextEntered && commandMode) {
                    if (event.text.unicode >= 32 && event.text.unicode < 128) {
                        currentCommand += static_cast<char>(event.text.unicode);
                    }
                }
            }

            if (!paused) {
                accumulator += frameDelta * simulationSpeed;
                double dt = m_sim.getTimeStep();
                int steps = 0;
                const int maxSteps = 500;
                while (accumulator >= dt && steps < maxSteps) {
                    m_sim.update();
                    accumulator -= dt;
                    ++steps;
                }
                if (steps >= maxSteps) accumulator = 0.0;
            }

            window.clear(sf::Color(5, 5, 30));

            for (const auto& body : m_sim.getBodies()) {
                body->draw(window, viewCenter, metersPerPixel);
            }

            if (fontLoaded) {
                for (const auto& body : m_sim.getBodies()) {
                    Vec2 offset = (body->getPosition() - viewCenter) / metersPerPixel;
                    float sx = static_cast<float>(offset.x) + window.getSize().x * 0.5f;
                    float sy = static_cast<float>(offset.y) + window.getSize().y * 0.5f;

                    sf::Text label;
                    label.setFont(font);
                    label.setString(body->getName());
                    label.setCharacterSize(11);
                    label.setFillColor(sf::Color(200, 200, 200, 180));
                    label.setPosition(sx + body->getScaledRadius(metersPerPixel) + 3.0f, sy - 6.0f);
                    window.draw(label);
                }

                drawHUD(window, font, paused, simulationSpeed, metersPerPixel,
                        commandMode, currentCommand);
            }

            window.display();
        }
    }

private:
    GravitySimulator m_sim;
    std::string m_lastFeedback;
    sf::Clock m_feedbackClock;

    void handleKeyPress(sf::Keyboard::Key key, bool& paused, bool& commandMode,
                        double& metersPerPixel, Vec2& viewCenter,
                        float& simulationSpeed, const sf::RenderWindow& window) {
        double panStep = metersPerPixel * 40.0;
        switch (key) {
            case sf::Keyboard::Space:  paused = !paused; break;
            case sf::Keyboard::Return: commandMode = true; break;
            case sf::Keyboard::Add:
            case sf::Keyboard::Equal:  metersPerPixel *= 0.9; break;
            case sf::Keyboard::Subtract:
            case sf::Keyboard::Dash:   metersPerPixel *= 1.1; break;
            case sf::Keyboard::Left:   viewCenter.x -= panStep; break;
            case sf::Keyboard::Right:  viewCenter.x += panStep; break;
            case sf::Keyboard::Up:     viewCenter.y -= panStep; break;
            case sf::Keyboard::Down:   viewCenter.y += panStep; break;
            case sf::Keyboard::Period:
                simulationSpeed = std::min(simulationSpeed * 2.0f, 1e9f);
                m_lastFeedback = "Speed: " + formatSI(simulationSpeed) + " sim-s/s";
                m_feedbackClock.restart();
                break;
            case sf::Keyboard::Comma:
                simulationSpeed = std::max(simulationSpeed * 0.5f, 1.0f);
                m_lastFeedback = "Speed: " + formatSI(simulationSpeed) + " sim-s/s";
                m_feedbackClock.restart();
                break;
            case sf::Keyboard::R:
                viewCenter = {0.0, 0.0};
                metersPerPixel = 2.0e9;
                break;
            case sf::Keyboard::Escape:
                break;
            default: break;
        }
    }

    void drawHUD(sf::RenderWindow& window, const sf::Font& font,
                 bool paused, float speed, double metersPerPixel,
                 bool commandMode, const std::string& currentCommand) {
        std::string methodName;
        switch (m_sim.getIntegrationMethod()) {
            case IntegrationMethod::EULER:           methodName = "Euler"; break;
            case IntegrationMethod::VELOCITY_VERLET: methodName = "Velocity Verlet"; break;
            case IntegrationMethod::RK4:             methodName = "Runge-Kutta 4"; break;
        }

        std::ostringstream info;
        info << "Bodies: " << m_sim.getBodies().size() << "\n"
             << "Integration: " << methodName << "\n"
             << "Time step: " << formatSI(m_sim.getTimeStep()) << " s\n"
             << "Speed: " << formatSI(speed) << " sim-s/s\n"
             << "Collisions: " << (m_sim.getCollisionDetection() ? "On" : "Off") << "\n"
             << "Scale: " << formatSI(metersPerPixel) << " m/px\n"
             << "Status: " << (paused ? "PAUSED" : "Running");

        sf::Text infoText;
        infoText.setFont(font);
        infoText.setCharacterSize(13);
        infoText.setFillColor(sf::Color(220, 220, 220));
        infoText.setPosition(10.f, 10.f);
        infoText.setString(info.str());

        if (commandMode) {
            sf::RectangleShape overlay(sf::Vector2f(
                static_cast<float>(window.getSize().x),
                static_cast<float>(window.getSize().y)));
            overlay.setFillColor(sf::Color(0, 0, 0, 190));
            window.draw(overlay);

            float panelW = 800.f;
            float panelH = 460.f;
            float panelX = (window.getSize().x - panelW) / 2.f;
            float panelY = (window.getSize().y - panelH) / 2.f;

            sf::RectangleShape panel(sf::Vector2f(panelW, panelH));
            panel.setFillColor(sf::Color(20, 20, 45, 248));
            panel.setOutlineColor(sf::Color(70, 70, 130));
            panel.setOutlineThickness(2.f);
            panel.setPosition(panelX, panelY);
            window.draw(panel);

            sf::Text title;
            title.setFont(font);
            title.setCharacterSize(20);
            title.setFillColor(sf::Color(180, 180, 255));
            title.setStyle(sf::Text::Bold);
            title.setString("Command Mode");
            sf::FloatRect titleBounds = title.getLocalBounds();
            title.setPosition(
                panelX + (panelW - titleBounds.width) / 2.f,
                panelY + 10.f);
            window.draw(title);

            struct CmdEntry { const char* cmd; const char* desc; };
            static const CmdEntry rows[] = {
                {"help",                                                          "Show all commands"},
                {"add <name> <mass> <posX> <posY> <velX> <velY> <R> <G> <B> <r>", "Add body (SI units: kg, m, m/s)"},
                {"remove <name>",                                                 "Remove a body"},
                {"list",                                                          "List all bodies with positions"},
                {"clear",                                                         "Remove all bodies"},
                {"solar",                                                         "Load solar system"},
                {"binary",                                                        "Load binary star system"},
                {"two-body",                                                      "Load simple two-body system"},
                {"method <euler|verlet|rk4>",                                     "Set integration method"},
                {"timestep <seconds>",                                            "Set physics time step"},
                {"speed <factor>",                                                "Set simulation speed multiplier"},
                {"zoom <metres_per_pixel>",                                       "Set zoom level directly"},
                {"center <x> <y>",                                                "Centre view on coordinates"},
                {"collisions <on|off>",                                           "Toggle collision detection"},
                {"trails <on|off>",                                               "Toggle motion trails"},
                {"traillength <n>",                                               "Set trail point count"},
                {"quit",                                                          "Exit"},
            };
            constexpr int rowCount = sizeof(rows) / sizeof(rows[0]);

            float col1X = panelX + 18.f;
            float col2X = panelX + 468.f;
            float tableTop = panelY + 42.f;

            sf::Text hdrCmd;
            hdrCmd.setFont(font);
            hdrCmd.setCharacterSize(12);
            hdrCmd.setStyle(sf::Text::Bold);
            hdrCmd.setFillColor(sf::Color(160, 180, 220));
            hdrCmd.setPosition(col1X, tableTop);
            hdrCmd.setString("Command");
            window.draw(hdrCmd);

            sf::Text hdrDesc;
            hdrDesc.setFont(font);
            hdrDesc.setCharacterSize(12);
            hdrDesc.setStyle(sf::Text::Bold);
            hdrDesc.setFillColor(sf::Color(160, 180, 220));
            hdrDesc.setPosition(col2X, tableTop);
            hdrDesc.setString("Description");
            window.draw(hdrDesc);

            float sepY = tableTop + 18.f;
            sf::RectangleShape sep(sf::Vector2f(panelW - 36.f, 1.f));
            sep.setFillColor(sf::Color(70, 70, 130));
            sep.setPosition(col1X, sepY);
            window.draw(sep);

            float rowY = sepY + 4.f;
            float rowH = 19.f;
            for (int i = 0; i < rowCount; ++i) {
                if (i % 2 == 0) {
                    sf::RectangleShape rowBg(sf::Vector2f(panelW - 36.f, rowH));
                    rowBg.setFillColor(sf::Color(35, 35, 60, 80));
                    rowBg.setPosition(col1X, rowY);
                    window.draw(rowBg);
                }

                sf::Text cmdCol;
                cmdCol.setFont(font);
                cmdCol.setCharacterSize(11);
                cmdCol.setFillColor(sf::Color(90, 210, 255));
                cmdCol.setPosition(col1X + 4.f, rowY + 2.f);
                cmdCol.setString(rows[i].cmd);
                window.draw(cmdCol);

                if (rows[i].desc[0] != '\0') {
                    sf::Text descCol;
                    descCol.setFont(font);
                    descCol.setCharacterSize(11);
                    descCol.setFillColor(sf::Color(200, 200, 210));
                    descCol.setPosition(col2X, rowY + 2.f);
                    descCol.setString(rows[i].desc);
                    window.draw(descCol);
                }

                rowY += rowH;
            }

            float inputY = rowY + 10.f;
            sf::RectangleShape inputBg(sf::Vector2f(panelW - 40.f, 28.f));
            inputBg.setFillColor(sf::Color(12, 12, 30));
            inputBg.setOutlineColor(sf::Color(100, 100, 180));
            inputBg.setOutlineThickness(1.f);
            inputBg.setPosition(panelX + 20.f, inputY);
            window.draw(inputBg);

            sf::Text cmdText;
            cmdText.setFont(font);
            cmdText.setCharacterSize(14);
            cmdText.setFillColor(sf::Color::Yellow);
            cmdText.setPosition(panelX + 26.f, inputY + 4.f);
            cmdText.setString("> " + currentCommand + "_");
            window.draw(cmdText);

            sf::Text hint;
            hint.setFont(font);
            hint.setCharacterSize(11);
            hint.setFillColor(sf::Color(120, 120, 140));
            hint.setPosition(panelX + 20.f, inputY + 36.f);
            hint.setString("Enter: Execute  |  Escape: Cancel");
            window.draw(hint);
        } else {
            window.draw(infoText);

            sf::Text helpText;
            helpText.setFont(font);
            helpText.setCharacterSize(12);
            helpText.setFillColor(sf::Color(160, 160, 160));
            helpText.setPosition(10.f, window.getSize().y - 75.f);
            helpText.setString(
                "Arrows: Pan  |  +/-/Scroll: Zoom  |  Space: Pause  |  R: Reset view\n"
                ",  Slow down  |  .  Speed up  |  Enter: Command mode");
            window.draw(helpText);

            if (!m_lastFeedback.empty() &&
                m_feedbackClock.getElapsedTime().asSeconds() < 4.0f) {
                float elapsed = m_feedbackClock.getElapsedTime().asSeconds();
                float alpha = (elapsed > 3.0f)
                    ? 255.0f * (4.0f - elapsed)
                    : 255.0f;

                sf::Text fbText;
                fbText.setFont(font);
                fbText.setCharacterSize(14);
                fbText.setFillColor(sf::Color(
                    100, 255, 100, static_cast<sf::Uint8>(alpha)));
                fbText.setPosition(10.f, window.getSize().y - 30.f);
                fbText.setString(m_lastFeedback);
                window.draw(fbText);
            }
        }
    }

    void setFeedback(const std::string& msg) {
        m_lastFeedback = msg;
        m_feedbackClock.restart();
        std::cout << msg << "\n";
    }

    void processCommand(const std::string& command, float& speed,
                        double& metersPerPixel, Vec2& viewCenter) {
        if (command.empty()) return;

        std::istringstream iss(command);
        std::string cmd;
        iss >> cmd;

        if (cmd == "help") {
            std::cout << "\n=== Gravity Simulator Commands ===\n"
                << "help                                           - Show this help\n"
                << "add <name> <mass> <posX> <posY> <velX> <velY> <R> <G> <B> <radius>\n"
                << "                                               - Add body (SI units)\n"
                << "remove <name>                                  - Remove body\n"
                << "list                                           - List all bodies\n"
                << "clear                                          - Remove all bodies\n"
                << "solar                                          - Solar system scenario\n"
                << "binary                                         - Binary star scenario\n"
                << "two-body                                       - Simple two-body\n"
                << "method <euler|verlet|rk4>                      - Integration method\n"
                << "timestep <seconds>                             - Physics time step\n"
                << "speed <factor>                                 - Simulation speed\n"
                << "zoom <meters_per_pixel>                        - Set zoom level\n"
                << "center <x_metres> <y_metres>                   - Center view\n"
                << "collisions <on|off>                            - Toggle collisions\n"
                << "trails <on|off>                                - Toggle trails\n"
                << "traillength <number>                           - Trail point count\n"
                << "quit                                           - Exit\n\n";
            setFeedback("Help printed to console.");
        }
        else if (cmd == "add") {
            std::string name;
            double mass, px, py, vx, vy;
            int cr, cg, cb;
            float radius;
            if (iss >> name >> mass >> px >> py >> vx >> vy >> cr >> cg >> cb >> radius) {
                m_sim.addBody(std::make_shared<CelestialBody>(
                    name, mass, Vec2(px, py), Vec2(vx, vy),
                    sf::Color(cr, cg, cb), radius));
                setFeedback("Added: " + name);
            } else {
                setFeedback("Usage: add <name> <mass> <posX> <posY> <velX> <velY> <R> <G> <B> <radius>");
            }
        }
        else if (cmd == "remove") {
            std::string name;
            iss >> name;
            setFeedback(m_sim.removeBody(name)
                ? "Removed: " + name
                : "Not found: " + name);
        }
        else if (cmd == "list") {
            std::cout << "\nBodies (" << m_sim.getBodies().size() << "):\n";
            for (const auto& b : m_sim.getBodies()) {
                auto& p = b->getPosition();
                auto& v = b->getVelocity();
                std::cout << "  " << b->getName()
                          << "  mass=" << formatSI(b->getMass()) << " kg"
                          << "  pos=(" << formatSI(p.x) << ", " << formatSI(p.y) << ") m"
                          << "  vel=(" << formatSI(v.x) << ", " << formatSI(v.y) << ") m/s\n";
            }
            setFeedback("Listed " + std::to_string(m_sim.getBodies().size()) + " bodies (see console).");
        }
        else if (cmd == "clear") {
            m_sim.clear();
            setFeedback("Cleared all bodies.");
        }
        else if (cmd == "solar") {
            Scenarios::loadSolarSystem(m_sim);
            metersPerPixel = 2.0e9;
            viewCenter = {0, 0};
            setFeedback("Loaded solar system.");
        }
        else if (cmd == "binary") {
            Scenarios::loadBinaryStars(m_sim);
            metersPerPixel = 1.0e9;
            viewCenter = {0, 0};
            setFeedback("Loaded binary stars.");
        }
        else if (cmd == "two-body") {
            Scenarios::loadTwoBody(m_sim);
            metersPerPixel = 1.0e9;
            viewCenter = {0, 0};
            setFeedback("Loaded two-body.");
        }
        else if (cmd == "method") {
            std::string m;
            iss >> m;
            if (m == "euler") {
                m_sim.setIntegrationMethod(IntegrationMethod::EULER);
                setFeedback("Method: Euler");
            } else if (m == "verlet") {
                m_sim.setIntegrationMethod(IntegrationMethod::VELOCITY_VERLET);
                setFeedback("Method: Velocity Verlet");
            } else if (m == "rk4") {
                m_sim.setIntegrationMethod(IntegrationMethod::RK4);
                setFeedback("Method: RK4");
            } else {
                setFeedback("Use: euler, verlet, rk4");
            }
        }
        else if (cmd == "timestep") {
            double dt;
            if (iss >> dt && dt > 0) {
                m_sim.setTimeStep(dt);
                setFeedback("Timestep: " + formatSI(dt) + " s");
            } else {
                setFeedback("Usage: timestep <positive_seconds>");
            }
        }
        else if (cmd == "speed") {
            float s;
            if (iss >> s && s > 0) {
                speed = s;
                setFeedback("Speed: " + formatSI(s) + " sim-s/s");
            } else {
                setFeedback("Usage: speed <positive_factor>");
            }
        }
        else if (cmd == "zoom") {
            double z;
            if (iss >> z && z > 0) {
                metersPerPixel = z;
                setFeedback("Zoom: " + formatSI(z) + " m/px");
            } else {
                setFeedback("Usage: zoom <meters_per_pixel>");
            }
        }
        else if (cmd == "center") {
            double cx, cy;
            if (iss >> cx >> cy) {
                viewCenter = {cx, cy};
                setFeedback("Centered at (" + formatSI(cx) + ", " + formatSI(cy) + ")");
            } else {
                setFeedback("Usage: center <x> <y>");
            }
        }
        else if (cmd == "collisions") {
            std::string opt; iss >> opt;
            if (opt == "on") {
                m_sim.setCollisionDetection(true);
                setFeedback("Collisions on.");
            } else if (opt == "off") {
                m_sim.setCollisionDetection(false);
                setFeedback("Collisions off.");
            } else {
                setFeedback("Use: on, off");
            }
        }
        else if (cmd == "trails") {
            std::string opt; iss >> opt;
            bool on = (opt == "on");
            for (const auto& b : m_sim.getBodies()) b->setDrawTrail(on);
            setFeedback(std::string("Trails ") + (on ? "on." : "off."));
        }
        else if (cmd == "traillength") {
            size_t len;
            if (iss >> len) {
                for (const auto& b : m_sim.getBodies()) b->setMaxTrailLength(len);
                setFeedback("Trail length: " + std::to_string(len));
            } else {
                setFeedback("Usage: traillength <number>");
            }
        }
        else if (cmd == "quit") {
            std::exit(0);
        }
        else {
            setFeedback("Unknown command: " + cmd + " (see command list)");
        }
    }
};

int main() {
    Application app;
    app.run();
    return 0;
}
