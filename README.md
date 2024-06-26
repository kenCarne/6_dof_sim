# 6_dof_sim
# Spacecraft Simulation

This project simulates the dynamics of a spacecraft in orbit around Earth, including its position, velocity, orientation, and angular velocity. The simulation uses Keplerian elements to initialize the spacecraft's state and then performs orbital maintenance maneuvers over time.

# Features

- **Vector3 Class**: Represents 3D vectors with basic vector operations such as addition, subtraction, scalar multiplication, dot product, cross product, magnitude, and normalization.
- **Quaternion Class**: Represents rotations using quaternions, with operations such as addition, multiplication, scalar multiplication, normalization, conjugation, and rotating vectors.
- **Keplerian to Cartesian Conversion**: Converts Keplerian elements (semi-major axis, eccentricity, inclination, etc.) to Cartesian coordinates (position and velocity).
- **Cartesian to Keplerian Conversion**: Converts Cartesian coordinates back to Keplerian elements.
- **Spacecraft Class**: Represents the spacecraft's state, including position, velocity, orientation, angular velocity, mass, moment of inertia, drag coefficient, reference area, atmospheric density, reflectivity, SRP area, distance from the sun, and delta-V. The spacecraft can apply forces and torques, and calculate gravitational, drag, and solar radiation pressure forces.
- **Simulation Loop**: Runs a simulation for a specified duration and time step, updating the spacecraft's state and performing maintenance maneuvers to correct its orbit.

# Dependencies

- Standard C++ Libraries: `iostream`, `fstream`, `cmath`

# Usage

1. **Compile the code**:
    ```bash
    g++ -o spacecraft_simulation spacecraft_simulation.cpp
    ```

2. **Run the simulation**:
    ```bash
    ./spacecraft_simulation
    ```

3. **Output**:
    - The simulation outputs the spacecraft's state at each time step to the console and saves it to a CSV file named `spacecraft_simulation.csv`. The CSV file includes the following columns:
        - Time
        - Position (X, Y, Z)
        - Orientation (W, X, Y, Z)
        - Velocity (X, Y, Z)
        - Angular Velocity (X, Y, Z)
        - Semi-Major Axis
        - Eccentricity
        - Inclination
        - RAAN (Right Ascension of Ascending Node)
        - Argument of Periapsis
        - Mean Anomaly
        - Delta-V

4. **Plot**:
  - plot the results with six_dof_sim_plotter.py
  - `spacecraft_simulation.csv`is referenced in file
