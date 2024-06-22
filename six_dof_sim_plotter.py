import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("C:\\Users\\extre\\Documents\\GitHub\\6_dof_sim\\six_dof_sim\\spacecraft_simulation.csv")

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

# Plot angular velocity data
plt.figure(figsize=(10, 6))
plt.plot(df["Time"], df)
