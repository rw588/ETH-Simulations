import numpy as np
import matplotlib.pyplot as plt

# Define constants
GRATING_PERIOD = 10
WAVELENGTH = 1
TALBOT_LENGTH = 2*GRATING_PERIOD**2 / WAVELENGTH  # This is an arbitrary value, you may replace it with your actual value
NUMBER_OF_PERIODS = 25  # Number of periods in simulation
STEPS = 1000  # Number of steps in simulation
Z0 = TALBOT_LENGTH * NUMBER_OF_PERIODS / 2  # We'll center the Gaussian along z at Z0
L = TALBOT_LENGTH * 3  # Width of the Gaussian drop off, adjust as needed

# Create 2D grid for x and z
x, z = np.linspace(0, NUMBER_OF_PERIODS * TALBOT_LENGTH, STEPS), np.linspace(0, NUMBER_OF_PERIODS * TALBOT_LENGTH, STEPS)
X, Z = np.meshgrid(x, z)

# Generate pattern
pattern = np.cos(2 * np.pi / GRATING_PERIOD * X) * np.cos(2 * np.pi / TALBOT_LENGTH * (Z - Z0)) * np.exp(-(Z - Z0)**2 / L**2)

# Define the third grating plane, starts upright
third_grating = np.ones_like(X)

# Tip and rotation angles
tip_angle = np.pi / 3  # Tip angle in radians, adjust as needed
rotation_angle = np.pi / 3  # Rotation angle in radians, adjust as needed

# Compute the rotation and tipping of the third grating using rotation matrices
R_z = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                [np.sin(rotation_angle), np.cos(rotation_angle), 0],
                [0, 0, 1]])

R_x = np.array([[1, 0, 0],
                [0, np.cos(tip_angle), -np.sin(tip_angle)],
                [0, np.sin(tip_angle), np.cos(tip_angle)]])

rotation_matrix = R_z @ R_x  # Matrix multiplication

# Apply rotation and tipping to the third grating
third_grating_rotated_and_tipped = np.zeros_like(X)
for i in range(STEPS):
    for j in range(STEPS):
        coordinates = np.array([X[i, j], 0, Z[i, j]])  # Original coordinates
        new_coordinates = rotation_matrix @ coordinates  # Apply rotation and tipping
        third_grating_rotated_and_tipped[i, j] = new_coordinates[2]  # We're interested in the z-coordinate

# Define variable transmittance across the third grating
transmittance = np.sin(2 * np.pi / GRATING_PERIOD * X)  # For example, sinusoidal variation

# Compute sampled_pattern
sampled_pattern = pattern * third_grating_rotated_and_tipped * transmittance

# Compute total signal received by the third grating at each z-coordinate
total_signal = np.sum(sampled_pattern, axis=1)
total_signal /= np.max(total_signal)
print(total_signal)

# Compute maximum and minimum values for each z-coordinate
max_val = np.max(sampled_pattern, axis=1)
min_val = np.min(sampled_pattern, axis=1)

# Compute SampledContrast
SampledContrast = (max_val - min_val) / (max_val + min_val)
print(SampledContrast)


# Compute total signal received by the third grating at each z-coordinate
total_signal = np.sum(sampled_pattern, axis=1)

# Plot the pattern, the sampled pattern, and the total signal
fig, ax = plt.subplots(3, 1, figsize=(10, 18))

pattern_plot = ax[0].contourf(X, Z, pattern, cmap='viridis')
ax[0].set_title('Pattern')
ax[0].set_xlabel('x')
ax[0].set_ylabel('z')
fig.colorbar(pattern_plot, ax=ax[0], orientation='vertical')

sampled_plot = ax[1].contourf(X, Z, sampled_pattern, cmap='viridis')
ax[1].set_title('Sampled Pattern')
ax[1].set_xlabel('x')
ax[1].set_ylabel('z')
fig.colorbar(sampled_plot, ax=ax[1], orientation='vertical')

ax[2].plot(z, total_signal)
ax[2].set_title('Total Signal')
ax[2].set_xlabel('z')
ax[2].set_ylabel('Signal')

plt.tight_layout()
plt.show()

