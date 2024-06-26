import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("C:\\Users\\extre\\source\\repos\\six_dof_sim\\spacecraft_simulation.csv")

# Plot position data
plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df["Position_X"], label="Position X")
plt.plot(df["Time"], df["Position_Y"], label="Position Y")
plt.plot(df["Time"], df["Position_Z"], label="Position Z")
plt.title("Spacecraft Position Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.legend()
plt.grid(True)
plt.show()

# Plot orientation data
plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df["Orientation_W"], label="Orientation W")
plt.plot(df["Time"], df["Orientation_X"], label="Orientation X")
plt.plot(df["Time"], df["Orientation_Y"], label="Orientation Y")
plt.plot(df["Time"], df["Orientation_Z"], label="Orientation Z")
plt.title("Spacecraft Orientation Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Orientation (quaternion components)")
plt.legend()
plt.grid(True)
plt.show()

# Plot velocity data
plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df["Velocity_X"], label="Velocity X")
plt.plot(df["Time"], df["Velocity_Y"], label="Velocity Y")
plt.plot(df["Time"], df["Velocity_Z"], label="Velocity Z")
plt.title("Spacecraft Velocity Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.legend()
plt.grid(True)
plt.show()

# Create subplots for Keplerian Elements
fig, axs = plt.subplots(6, 1, figsize=(10, 36))

# Plot SemiMajorAxis
axs[0].plot(df["Time"], df["SemiMajorAxis"], label="SemiMajorAxis")
axs[0].set_title("SemiMajorAxis Over Time")
axs[0].set_xlabel("Time (s)")
axs[0].set_ylabel("SemiMajorAxis (km)")
axs[0].legend()
axs[0].grid(True)

# Plot Eccentricity
axs[1].plot(df["Time"], df["Eccentricity"], label="Eccentricity")
axs[1].set_title("Eccentricity Over Time")
axs[1].set_xlabel("Time (s)")
axs[1].set_ylabel("Eccentricity")
axs[1].legend()
axs[1].grid(True)

# Plot Inclination
axs[2].plot(df["Time"], df["Inclination"], label="Inclination")
axs[2].set_title("Inclination Over Time")
axs[2].set_xlabel("Time (s)")
axs[2].set_ylabel("Inclination (degrees)")
axs[2].legend()
axs[2].grid(True)

# Plot Right Ascension of the Ascending Node (RAAN)
axs[3].plot(df["Time"], df["RAAN"], label="RAAN")
axs[3].set_title("RAAN Over Time")
axs[3].set_xlabel("Time (s)")
axs[3].set_ylabel("RAAN (degrees)")
axs[3].legend()
axs[3].grid(True)

# Plot Argument of Periapsis
axs[4].plot(df["Time"], df["ArgumentOfPeriapsis"], label="Argument of Periapsis")
axs[4].set_title("Argument of Periapsis Over Time")
axs[4].set_xlabel("Time (s)")
axs[4].set_ylabel("Argument of Periapsis (degrees)")
axs[4].legend()
axs[4].grid(True)

# Plot Mean Anomaly
axs[5].plot(df["Time"], df["MeanAnomaly"], label="Mean Anomaly")
axs[5].set_title("Mean Anomaly Over Time")
axs[5].set_xlabel("Time (s)")
axs[5].set_ylabel("Mean Anomaly (degrees)")
axs[5].legend()
axs[5].grid(True)

# Adjust layout
plt.tight_layout()
plt.show()

# Plot angular velocity data
plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df["AngularVelocity_X"], label="Angular Velocity X")
plt.plot(df["Time"], df["AngularVelocity_Y"], label="Angular Velocity Y")
plt.plot(df["Time"], df["AngularVelocity_Z"], label="Angular Velocity Z")
plt.title("Spacecraft Angular Velocity Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Angular Velocity (rad/s)")
plt.legend()
plt.grid(True)
plt.show()

# Create subplots for deltaV and its components
fig, axs = plt.subplots(4, 1, figsize=(10, 24))

# Plot deltaV
axs[0].plot(df["Time"], df["DeltaV"], label="DeltaV")
axs[0].set_title("DeltaV Over Time")
axs[0].set_xlabel("Time (s)")
axs[0].set_ylabel("DeltaV (m/s)")
axs[0].legend()
axs[0].grid(True)

# Plot deltaV components
axs[1].plot(df["Time"], df["DeltaV_X"], label="DeltaV X")
axs[1].set_title("DeltaV X Over Time")
axs[1].set_xlabel("Time (s)")
axs[1].set_ylabel("DeltaV X (m/s)")
axs[1].legend()
axs[1].grid(True)

axs[2].plot(df["Time"], df["DeltaV_Y"], label="DeltaV Y")
axs[2].set_title("DeltaV Y Over Time")
axs[2].set_xlabel("Time (s)")
axs[2].set_ylabel("DeltaV Y (m/s)")
axs[2].legend()
axs[2].grid(True)

axs[3].plot(df["Time"], df["DeltaV_Z"], label="DeltaV Z")
axs[3].set_title("DeltaV Z Over Time")
axs[3].set_xlabel("Time (s)")
axs[3].set_ylabel("DeltaV Z (m/s)")
axs[3].legend()
axs[3].grid(True)

# Adjust layout
plt.tight_layout()
plt.show()
