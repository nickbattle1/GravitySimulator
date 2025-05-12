# Gravity Simulator

A C++ physics simulation program that visualizes gravitational interactions between celestial bodies using Newton's laws of universal gravitation.

## Features

- Real-time simulation of gravitational interactions
- Multiple numerical integration methods:
  - Euler
  - Velocity Verlet
  - Runge-Kutta 4 (RK4)
- Interactive visualization using SFML graphics library
- Pre-configured simulations:
  - Solar system
  - Two-body system
  - Binary star system
- Adjustable simulation parameters:
  - Time step
  - Simulation speed
  - Collision detection
- Motion trails for tracking object paths
- Command-line interface for real-time control

## Controls

- **Arrows**: Pan view
- **+/-**: Zoom in/out
- **R**: Reset view
- **Space**: Pause/Resume simulation
- **Enter**: Enter command mode
- **Escape**: Exit command mode or quit application

## Command Mode

Enter command mode by pressing Enter during simulation. Available commands include:

- `help`: Show all available commands
- `add <name> <mass> <posX> <posY> <velX> <velY> <colorR> <colorG> <colorB> <radius>`: Add new body
- `remove <name>`: Remove a body by name
- `list`: List all bodies in the simulation
- `clear`: Remove all bodies
- `solar`: Initialize solar system
- `binary`: Initialize binary star system
- `two-body`: Initialize simple two-body system
- `method <euler|verlet|rk4>`: Set integration method
- `timestep <seconds>`: Set simulation time step
- `speed <factor>`: Set simulation speed multiplier
- `collisions <on|off>`: Enable/disable collision detection
- `trails <on|off>`: Enable/disable trails for all bodies
- `traillength <number>`: Set maximum trail length
- `quit`: Exit the simulation

## Technical Details

- The simulation uses a scaled version of the universal gravitational constant (G) to make the visualization practical
- Bodies are represented as colored circles with size proportional to their radius
- Physics calculations include:
  - Newton's law of universal gravitation (F = G * m1 * m2 / r²)
  - Numerical integration of equations of motion
  - Optional handling of collisions with conservation of momentum
- Motion trails display the path of objects with a fading effect

## Requirements

- C++ compiler with C++11 support
- SFML graphics library
- Recommended: Arial font (arial.ttf) for text display

## Building

Compile with g++ or your preferred C++ compiler, linking against the SFML library:

```bash
g++ -o gravitysim main.cpp -lsfml-graphics -lsfml-window -lsfml-system
```

## Usage

Run the compiled executable:

```bash
./gravitysim
```

The simulation starts with the solar system by default. 