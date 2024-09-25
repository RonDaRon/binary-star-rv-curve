import numpy as np
import matplotlib.pyplot as plt

# Constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M_sun = 1.989e30  # Solar mass (kg)
R_sun = 6.957e8
year_in_seconds = 3.154e7  # One year in seconds

# Orbital parameters
M1 = 10.0 * M_sun  # Mass of star 1 (kg)
M2 = 5.0 * M_sun  # Mass of star 2 (kg)
a = 250.0 * R_sun  # Semi-major axis (m)
e = 0.5  # Eccentricity
inclination = np.radians(90)  # Inclination in radians
omega = 0.8  # Argument of periapsis in radians
P = 2 * np.pi * np.sqrt(a**3 / (G * (M1 + M2)))  # Orbital period (seconds)
print('P = {0} [days]'.format(np.round(P/(60*60*24), 2)))
T = P / year_in_seconds  # Orbital period in years

# Function to compute the radial velocity amplitude (K)
def velocity_amplitude(M, P, e, i, M_total):
    return (2 * np.pi * G / P) ** (1/3) * M / (M_total ** (2/3)) * np.sin(i) / np.sqrt(1 - e**2)

# Compute K1 and K2 (velocity amplitudes)
K1 = velocity_amplitude(M2, P, e, inclination, M1 + M2)
K2 = velocity_amplitude(M1, P, e, inclination, M1 + M2)

# Function to compute the true anomaly
def true_anomaly(E, e):
    return 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

# Function to compute the radial velocity for each star.
def radial_velocity(K, e, omega, time, P):
    M_mean = 2 * np.pi * time / P  # Mean anomaly
    E = M_mean  # Eccentric anomaly (simplified, without solving Kepler's equation)
    theta = true_anomaly(E, e)  # True anomaly
    return K * (np.cos(omega + theta) + e * np.cos(omega))

# Time array over two complete orbits.
time = np.linspace(0, 2 * T, 1000)  # in years

# Radial velocities for both stars.
v1 = radial_velocity(K1, e, omega, time * year_in_seconds, P)
v2 = radial_velocity(K2, e, omega + np.pi, time * year_in_seconds, P)

# Convert to km/s.
v1_kms = v1 / 1000
v2_kms = v2 / 1000

# Plot the radial velocity curves for both stars.
plt.figure(figsize=(10, 6))
plt.plot(time, v1_kms, label="Star 1", color='red')
plt.plot(time, v2_kms, label="Star 2", color='blue')
plt.axhline(0, color='black', linewidth=0.5)
plt.title("Radial Velocity Curves for a Binary Star System")
plt.xlabel("Time (years)")
plt.ylabel("Radial Velocity (km/s)")
plt.yticks(np.linspace(-200, 100, 7))
plt.legend()
plt.grid(False)
plt.savefig('binary.png')
plt.close()