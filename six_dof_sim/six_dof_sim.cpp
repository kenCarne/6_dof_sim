#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>



// Vector3 class to represent 3D vectors
class Vector3 {
public:
    double x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vector3 operator+(const Vector3& v) const {
        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    Vector3 operator-(const Vector3& v) const {
        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    Vector3 operator*(double scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    Vector3 cross(const Vector3& v) const {
        return Vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    double dot(const Vector3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vector3 normalize() const {
        double mag = magnitude();
        return Vector3(x / mag, y / mag, z / mag);
    }
};

// Quaternion class to represent rotations
class Quaternion {
public:
    double w, x, y, z;

    Quaternion() : w(1), x(0), y(0), z(0) {}
    Quaternion(double w, double x, double y, double z) : w(w), x(x), y(y), z(z) {}

    Quaternion operator+(const Quaternion& q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }

    Quaternion operator*(const Quaternion& q) const {
        return Quaternion(
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        );
    }

    Quaternion operator*(double scalar) const {
        return Quaternion(w * scalar, x * scalar, y * scalar, z * scalar);
    }

    Quaternion normalize() const {
        double mag = std::sqrt(w * w + x * x + y * y + z * z);
        return Quaternion(w / mag, x / mag, y / mag, z / mag);
    }

    Vector3 rotate(const Vector3& v) const {
        Quaternion qv(0, v.x, v.y, v.z);
        Quaternion qr = (*this) * qv * conjugate();
        return Vector3(qr.x, qr.y, qr.z);
    }

    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }
};

// Function to convert Keplerian elements to position and velocity vectors
void keplerianToCartesian(double a, double e, double i, double omega, double w, double M, double mu,
    Vector3& position, Vector3& velocity) {
    // Convert degrees to radians
    i *= M_PI / 180.0;
    omega *= M_PI / 180.0;
    w *= M_PI / 180.0;
    M *= M_PI / 180.0;

    // Solve Kepler's equation for E (eccentric anomaly) using Newton-Raphson
    double E = M;
    double E_prev;
    do {
        E_prev = E;
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
    } while (fabs(E - E_prev) > 1e-6);

    // Calculate true anomaly
    double v = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    // Calculate distance
    double r = a * (1 - e * cos(E));

    // Position in orbital plane
    double x_orb = r * cos(v);
    double y_orb = r * sin(v);

    // Velocity in orbital plane
    double vx_orb = sqrt(mu / a) * -sin(E);
    double vy_orb = sqrt(mu / a) * sqrt(1 - e * e) * cos(E);

    // Rotate to inertial frame
    double cos_w = cos(w);
    double sin_w = sin(w);
    double cos_i = cos(i);
    double sin_i = sin(i);
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);

    position.x = (cos_omega * cos_w - sin_omega * sin_w * cos_i) * x_orb +
        (-cos_omega * sin_w - sin_omega * cos_w * cos_i) * y_orb;
    position.y = (sin_omega * cos_w + cos_omega * sin_w * cos_i) * x_orb +
        (-sin_omega * sin_w + cos_omega * cos_w * cos_i) * y_orb;
    position.z = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb;

    velocity.x = (cos_omega * cos_w - sin_omega * sin_w * cos_i) * vx_orb +
        (-cos_omega * sin_w - sin_omega * cos_w * cos_i) * vy_orb;
    velocity.y = (sin_omega * cos_w + cos_omega * sin_w * cos_i) * vx_orb +
        (-sin_omega * sin_w + cos_omega * cos_w * cos_i) * vy_orb;
    velocity.z = (sin_w * sin_i) * vx_orb + (cos_w * sin_i) * vy_orb;
}

// Spacecraft class to represent the spacecraft state
class Spacecraft {
public:
    Vector3 position;
    Vector3 velocity;
    Quaternion orientation;
    Vector3 angularVelocity;
    double mass;
    Vector3 momentOfInertia;

    Spacecraft()
        : position(), velocity(), orientation(), angularVelocity(),
        mass(1.0), momentOfInertia(1.0, 1.0, 1.0) {}

    void applyForce(const Vector3& force, double dt) {
        Vector3 acceleration = force * (1.0 / mass);
        velocity = velocity + acceleration * dt;
        position = position + velocity * dt;
    }

    void applyTorque(const Vector3& torque, double dt) {
        Vector3 angularAcceleration = Vector3(
            torque.x / momentOfInertia.x,
            torque.y / momentOfInertia.y,
            torque.z / momentOfInertia.z
        );
        angularVelocity = angularVelocity + angularAcceleration * dt;

        Quaternion deltaOrientation = Quaternion(
            0,
            angularVelocity.x * dt * 0.5,
            angularVelocity.y * dt * 0.5,
            angularVelocity.z * dt * 0.5
        ) * orientation;

        orientation = (orientation + deltaOrientation).normalize();
    }

    void update(double dt) {
        position = position + velocity * dt;
        Quaternion deltaOrientation = Quaternion(
            0,
            angularVelocity.x * dt * 0.5,
            angularVelocity.y * dt * 0.5,
            angularVelocity.z * dt * 0.5
        ) * orientation;
        orientation = (orientation + deltaOrientation).normalize();
    }
};

int main() {
    Spacecraft spacecraft;
    std::ofstream outFile("spacecraft_simulation.csv");
    outFile << "Time,Position_X,Position_Y,Position_Z,Orientation_W,Orientation_X,Orientation_Y,Orientation_Z,Velocity_X,Velocity_Y,Velocity_Z,AngularVelocity_X,AngularVelocity_Y,AngularVelocity_Z\n";

    // Keplerian elements (example values)
    double a = 7000;         // Semi-major axis (km)
    double e = 0.001;        // Eccentricity
    double i = 98.7;         // Inclination (degrees)
    double omega = 257.7;    // Right ascension of the ascending node (degrees)
    double w = 0.0;          // Argument of periapsis (degrees)
    double M = 0.0;          // Mean anomaly (degrees)
    double mu = 398600.4418; // Standard gravitational parameter for Earth (km^3/s^2)

    keplerianToCartesian(a, e, i, omega, w, M, mu, spacecraft.position, spacecraft.velocity);

    // Simulate for 10 seconds with a time step of 0.1 seconds
    double simulationTime = 10.0;
    double timeStep = 0.1;

    for (double t = 0; t < simulationTime; t += timeStep) {
        // Apply some arbitrary force and torque
        Vector3 force(0, 0, 1);
        Vector3 torque(0.01, 0.01, 0);

        spacecraft.applyForce(force, timeStep);
        spacecraft.applyTorque(torque, timeStep);
        spacecraft.update(timeStep);

        // Print the spacecraft state
        std::cout << "Time: " << t << " s\n";
        std::cout << "Position: (" << spacecraft.position.x << ", "
            << spacecraft.position.y << ", " << spacecraft.position.z << ")\n";
        std::cout << "Orientation: (" << spacecraft.orientation.w << ", "
            << spacecraft.orientation.x << ", " << spacecraft.orientation.y << ", "
            << spacecraft.orientation.z << ")\n";
        std::cout << "Velocity: (" << spacecraft.velocity.x << ", "
            << spacecraft.velocity.y << ", " << spacecraft.velocity.z << ")\n";
        std::cout << "Angular Velocity: (" << spacecraft.angularVelocity.x << ", "
            << spacecraft.angularVelocity.y << ", " << spacecraft.angularVelocity.z << ")\n";
        std::cout << "-------------------------\n";

        // Save the spacecraft state to CSV file
        outFile << t << ","
            << spacecraft.position.x << ","
            << spacecraft.position.y << ","
            << spacecraft.position.z << ","
            << spacecraft.orientation.w << ","
            << spacecraft.orientation.x << ","
            << spacecraft.orientation.y << ","
            << spacecraft.orientation.z << ","
            << spacecraft.velocity.x << ","
            << spacecraft.velocity.y << ","
            << spacecraft.velocity.z << ","
            << spacecraft.angularVelocity.x << ","
            << spacecraft.angularVelocity.y << ","
            << spacecraft.angularVelocity.z << "\n";
    }

    outFile.close();
    return 0;
}

//TODO

//add Sensors -> create random matrix generator size of tolerance for sensors
//add forces -> Drag, SRP, Gravity 
//incorporate obital maintence burns
//add EKF for sensors and perabations
//add different refernce frames
//gui to select data, enter satellite starting point ect? -> becomes python script, subprocess run C++ six_dof_sim.cpp get output -> plot, gui is plot interface
