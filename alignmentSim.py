import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure



# Define constants
GRATING_PERIOD = 100
WAVELENGTH = 5
TALBOT_LENGTH = GRATING_PERIOD**2 / WAVELENGTH  # This is an arbitrary value, you may replace it with your actual value
print('talbot: ', TALBOT_LENGTH)
print('wavelength: ', WAVELENGTH)
print('grating: ', GRATING_PERIOD)
NUMBER_OF_PERIODS = 20  # Number of periods in simulation
STEPS = 1000  # Number of steps in simulation
Z0 = TALBOT_LENGTH*NUMBER_OF_PERIODS/2
L = TALBOT_LENGTH * 2

# Create 2D grid for x and z
x, y, z = np.linspace(0, NUMBER_OF_PERIODS * GRATING_PERIOD, STEPS), np.linspace(0, NUMBER_OF_PERIODS * GRATING_PERIOD, STEPS), np.linspace(0, NUMBER_OF_PERIODS * TALBOT_LENGTH, STEPS)
X, Y, Z = np.meshgrid(x, y, z)

# Generate pattern
def Pattern(x, y, z):
    #independant of y
    return np.cos(2 * np.pi / GRATING_PERIOD * x)**2 * np.cos(2 * np.pi / TALBOT_LENGTH * (z - Z0))**2 * np.exp(-(z- Z0)**2 / L**2)

# Define the third grating plane, starts upright
third_grating = np.ones_like(X)
#print('third grating: ', third_grating)

# Tip and rotation angles
tip_angle = np.pi / 6  # Tip angle in radians, adjust as needed
rotation_angle_x = np.pi / 12  # Rotation around x-axis in radians, adjust as needed
rotation_angle_y = np.pi / 24  # Rotation around y-axis in radians, adjust as needed

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

rotation_matrix = R_y @ R_x @ R_z  # Matrix multiplication

print('rotation matrix: ', rotation_matrix)

###############################################################################################

#start with a perfectly aligned plane: z = 0 #normal along the z axis

plane = [0, 0, 1]

rotated_plane = rotation_matrix @ plane

print('rotated plane', rotated_plane)

def plane_to_grid(coeff, xlim, ylim, steps):
    A, B, C = coeff
    x = np.linspace(-xlim, xlim, steps)
    y = np.linspace(-ylim, ylim, steps)
    X, Y = np.meshgrid(x, y)
    Z = -(A*X + B*Y) / C
    return Z

xlim = NUMBER_OF_PERIODS * GRATING_PERIOD
ylim = NUMBER_OF_PERIODS * GRATING_PERIOD
steps = 10000

#create plane array
Z = plane_to_grid(rotated_plane, xlim, ylim, steps)

fig = Figure()
canvas = FigureCanvas(fig)
plt.imshow(Z)

plt.show()

###############################################################################################
'''
# Apply rotation and tipping to the third grating
third_grating_rotated_and_tipped = np.zeros_like(X)
for i in range(STEPS):
    for j in range(STEPS):
        coordinates = np.array([X[i, j], 0, Z[i, j]])  # Original coordinates
        new_coordinates = rotation_matrix @ coordinates  # Apply rotation and tipping
        third_grating_rotated_and_tipped[i, j] = new_coordinates[2]  # We're interested in the z-coordinate

# Define variable transmittance across the third grating
# Assuming X and GRATING_PERIOD are defined
transmittance = np.sin(2 * np.pi / GRATING_PERIOD * Z)  # For example, sinusoidal variation

# If transmittance is a 2D array, apply np.where operation along axis=0 (vertically)
transmittance_shifted = np.roll(transmittance, shift=-1, axis=0) # shift values up
transmittance_shifted[-1] = transmittance[-1] # last row should be same as original

transmittance = np.where(transmittance_shifted >= 0, 1, 0)

plt.imshow(transmittance)
print(transmittance)

#fig2 = plt.figure()
#plt.plot(transmittance[1:100,1])
print('mean: ', np.mean(transmittance[:]))
#plt.show()


# Now we simulate the "sampling" by the third grating plane
sampled_pattern = pattern * third_grating_rotated_and_tipped * transmittance

# Compute total signal received by the third grating at each z-coordinate
total_signal = np.sum(sampled_pattern, axis=1)
total_signal /= max(total_signal)
# Compute maximum and minimum values for each z-coordinate
max_val = np.max(sampled_pattern, axis=1)
min_val = np.min(sampled_pattern, axis=1)

# Compute SampledContrast
SampledContrast = (max_val - min_val) / (max_val + min_val)

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

'''