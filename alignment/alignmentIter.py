import numpy as np
import matplotlib.pyplot as plt

# Define constants
GRATING_PERIOD = 100
WAVELENGTH = 5
TALBOT_LENGTH = GRATING_PERIOD**2 / WAVELENGTH
NUMBER_OF_PERIODS = 20
STEPS = 1000
Z0 = TALBOT_LENGTH*NUMBER_OF_PERIODS/2
L = TALBOT_LENGTH * 5

# Create 2D grid for x and z
x, y, z = np.linspace(0, NUMBER_OF_PERIODS * GRATING_PERIOD, STEPS), np.linspace(0, NUMBER_OF_PERIODS * GRATING_PERIOD, STEPS), np.linspace(0, NUMBER_OF_PERIODS * TALBOT_LENGTH, STEPS)
X, Y, Z = np.meshgrid(x, y, z)

# Generate pattern
def Pattern(x, y, z):
    return np.cos(2 * np.pi / GRATING_PERIOD * x)**2 * np.cos(2 * np.pi / TALBOT_LENGTH * (z - Z0))**2 * np.exp(-(z - Z0)**2 / L**2)

# Define the third grating plane, starts upright
third_grating = np.ones_like(X)

# Tip and rotation angles
tip_angle = 0*np.pi / 6
rotation_angle_x = 0*np.pi / 12
rotation_angle_y = np.pi / 3

# Compute the rotation and tipping of the third grating using rotation matrices
R_z = np.array([[np.cos(rotation_angle_x), -np.sin(rotation_angle_x), 0],
                [np.sin(rotation_angle_x), np.cos(rotation_angle_x), 0],
                [0, 0, 1]])

R_x = np.array([[1, 0, 0],
                [0, np.cos(tip_angle), -np.sin(tip_angle)],
                [0, np.sin(tip_angle), np.cos(tip_angle)]])

R_y = np.array([[np.cos(rotation_angle_y), 0, np.sin(rotation_angle_y)],
                [0, 1, 0],
                [-np.sin(rotation_angle_y), 0, np.cos(rotation_angle_y)]])

rotation_matrix = R_y @ R_x @ R_z

# Start with a perfectly aligned plane: z = 0 #normal along the z axis
plane = [0, 0, 1]
rotated_plane = rotation_matrix @ plane

def plane_to_grid(coeff, xlim, ylim, steps):
    A, B, C = coeff
    x = np.linspace(-xlim, xlim, steps)
    y = np.linspace(-ylim, ylim, steps)
    X, Y = np.meshgrid(x, y)
    Z = -(A*X + B*Y) / C
    return Z, X, Y

xlim = NUMBER_OF_PERIODS * GRATING_PERIOD
ylim = NUMBER_OF_PERIODS * GRATING_PERIOD
steps = 10000

# Create plane array
Z, X, Y = plane_to_grid(rotated_plane, xlim, ylim, steps)

z_values = np.arange(0, NUMBER_OF_PERIODS * TALBOT_LENGTH, TALBOT_LENGTH / 25)
mean_pattern_values = []

for z_val in z_values:
    pattern = Pattern(X, Y, z_val)
    mean_pattern_values.append(np.mean(pattern))

# Plot mean pattern values vs z values
plt.plot(z_values, mean_pattern_values)
plt.xlabel('z value')
plt.ylabel('Mean pattern value')
plt.show()
