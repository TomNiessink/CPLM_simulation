import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import rotate


# Definition of optical elements
H = np.array([[1, 0], [0, 0]])  # Horizontal Polarizer [1]
V = np.array([[0, 0], [0, 1]])  # Vertical Polarizer [1]
RHCP = 0.5 * np.array([[1, 1j], [-1j, 1]])  # Right Hand Circular Polarizer [1], not used
LHCP = 0.5 * np.array([[1, -1j], [1j, 1]])  # Left Hand Circular Polarizer [1], not used
QWPH = np.exp(1j * np.pi / 4) * np.array([[1, 0], [0, -1j]])  # Quarter wave plate horizontal [1], not used
QWPV = np.exp(-1j * np.pi / 4) * np.array([[1, 0], [0, 1j]])  # Quarter wave plate vertical [1], not used

# Linear polarizer (angular)
def LP(theta):
    return np.array([[np.cos(theta)**2, np.cos(theta) * np.sin(theta)],
                     [np.cos(theta) * np.sin(theta), np.sin(theta)**2]])

# Arbitrary birefringent material as phase retarder (Linear) [2]
def BFO(theta, eta):
    return np.exp(-1j * eta / 2) * np.array([[np.cos(theta)**2 + np.exp(1j * eta) * np.sin(theta)**2,
                                               (1 - np.exp(1j * eta)) * np.cos(theta) * np.sin(theta)],
                                              [(1 - np.exp(1j * eta)) * np.cos(theta) * np.sin(theta),
                                               np.sin(theta)**2 + np.exp(1j * eta) * np.cos(theta)**2]])

def rod(r, l, theta, Gridsize_m, Gridsize_px):
    # Function which generates a rod-shaped crystal with
    # Length l (m), radius r (m), angle theta to the optical axis (rad),
    # and a defined grid (Gridsize_m, m and Gridsize_px, px)

    H = np.zeros((Gridsize_px, Gridsize_px))
    res = Gridsize_m / Gridsize_px

    crys = np.zeros((round(2 * r / res) + 1, round(l / res)))
    angles = np.zeros((round(2 * r / res) + 1, round(l / res)))
    angles[:, :] = theta

    x = np.arange(-r, r + res, res)
    h_x = 2 * r * np.sin(np.arccos(x / r))
    for i in range(round(l / res)):
        crys[:, i] = h_x

    crys = rotate(crys, np.degrees(theta), reshape=True)
    angles = rotate(angles, np.degrees(theta), reshape=True)

    szcrys = crys.shape
    szH = H.shape
    angH = np.copy(H)

    H[(szH[0] // 2 - szcrys[0] // 2):(szH[0] // 2 - szcrys[0] // 2 + szcrys[0]),
      (szH[1] // 2 - szcrys[1] // 2):(szH[1] // 2 - szcrys[1] // 2 + szcrys[1])] = crys

    angH[(szH[0] // 2 - szcrys[0] // 2):(szH[0] // 2 - szcrys[0] // 2 + szcrys[0]),
         (szH[1] // 2 - szcrys[1] // 2):(szH[1] // 2 - szcrys[1] // 2 + szcrys[1])] = angles

    return H, angH


def MCross(r, ne, no, Gridsize_m, Gridsize_px):
    # Function which generates a Maltese cross birefringent object
    # With radius r, and a birefingence dn which is dependent on a given ne,
    # no. Output is an angle table and a retardance table H, which is equal
    # to the thickness of the sphere * local birefringence.
    
    res = Gridsize_m / Gridsize_px
    lb = -1 * Gridsize_m / 2
    rb = Gridsize_m / 2

    x = np.arange(lb, rb + res, res)  # m
    y = np.arange(lb, rb + res, res)  # m
    z = np.arange(lb, rb + res, res)  # m

    r_ij = np.zeros((len(x), len(y)))
    d_ij = np.zeros((len(x), len(y)))
    ang_ij = np.zeros((len(x), len(y)))
    ne_p_ijk = np.zeros((len(x), len(y), len(z)))
    dne_p_ijk = np.zeros((len(x), len(y), len(z)))
    r_ijk = np.zeros((len(x), len(y), len(z)))
    theta_ijk = np.zeros((len(x), len(y), len(z)))

    for i in range(len(x)):
        for j in range(len(y)):
            ang_ij[i, j] = np.arctan2(y[j], x[i])
            for k in range(len(z)):
                r_ij[i, j] = np.sqrt(x[i]**2 + y[j]**2)
                d_ij[i, j] = 2 * np.sqrt(abs(r**2 - r_ij[i, j]**2))
                r_ijk[i, j, k] = np.sqrt(z[k]**2 + r_ij[i, j]**2)

                if r_ijk[i, j, k] > r:
                    ne_p_ijk[i, j, k] = 0
                else:
                    theta_ijk[i, j, k] = np.pi / 2 - np.arctan((d_ij[i, j] / 2 - z[k]) / r_ij[i, j])
                    ne_p_ijk[i, j, k] = 1 / np.sqrt((np.cos(theta_ijk[i, j, k])**2 / no**2) +
                                                     (np.sin(theta_ijk[i, j, k])**2 / ne**2))
                    dne_p_ijk[i, j, k] = no - ne_p_ijk[i, j, k]

    H = np.sum(dne_p_ijk, axis=2) * res
    H[np.isnan(H)] = 0
    angH = ang_ij
    angH[np.isnan(angH)] = 0


    return H, angH


# Load LED spectrum
LED_spectrum = pd.read_csv(r'PATH_HERE', delimiter="\t")
LED_spectrum = np.array(LED_spectrum)
lambda_vals = LED_spectrum[:, 0] / 1E9  # Wavelengths array
intensities = LED_spectrum[:, 1]  # Intensities per wavelength
intensities /= max(intensities)  # Normalize intensities

# Downsample LED spectrum
lambda_vals = lambda_vals[::4]
intensities = intensities[::4]

# Define crystal
Gridsize_px = 100  # Gridsize in pixels
Gridsize_m = 10E-6  # Size of pixels in grid

# For Maltese Cross
#BF = 1
#r = 2.5E-6
#ne = 1.46
#no = 1.40
#d, ang = MCross(r, ne, no, Gridsize_m, Gridsize_px)

# For Rod-like crystal
r = 1E-6
l = 5E-6
BF = 0.015
theta = np.pi/4
d, ang = rod(r,l,theta,Gridsize_m,Gridsize_px)


# Calculate image
I_out_Hg = np.zeros((len(lambda_vals), 101, 101))  # Preallocate data to save time

Sz_Hg = d.shape  # The size of the image

for i in range(len(lambda_vals)):  # Loop over all wavelengths in the wavelength table
    lam = lambda_vals[i]
    I_in = np.ones(2) * intensities[i]  # Retrieve intensity of certain wavelength
    H_dphi = (2 * np.pi / lam) * (BF) * d  # Calculate the retardance for the full image
    dphiqwl = (2 * np.pi / lam) * 550E-9  # Calculate the retardance for the additional quarter lambda plate

    for j in range(Sz_Hg[0]):  # Loop over one axis
        for k in range(Sz_Hg[1]):  # Loop over other axis
            dphi = H_dphi[j, k]  # Take the retardance of the current pixel
            # Calculate output
            #print("Shapes:", I_in.shape, H.shape, BFO(ang[j, k], dphi).shape, BFO(np.pi / 4, dphiqwl).shape, LP(np.pi / 2).shape)
            I_out_Hg[i, j, k] = np.linalg.norm(I_in @ H @ BFO(ang[j, k], dphi) @ BFO(np.pi / 4, dphiqwl) @ LP(np.pi / 2))


# Color analysis
image = np.zeros((Sz_Hg[0], (Sz_Hg[1]), 3))  # Preallocate data space


for j in range(Sz_Hg[0]):
    for k in range(Sz_Hg[1]):
        Spectrum = I_out_Hg[:, j, k]
        lb=np.where(lambda_vals == 613E-9)
        lb=lb[0][0]
        rb=np.where(lambda_vals == 801E-9)
        rb=rb[0][0]
        image[j, k, 0] = np.sum(Spectrum[lb:rb])  # Red bin
        lb=np.where(lambda_vals == 513E-9)
        lb=lb[0][0]
        rb=np.where(lambda_vals == 601E-9)
        rb=rb[0][0]
        image[j, k, 1] = np.sum(Spectrum[lb:rb])#np.sum(Spectrum[np.where(lambda_vals == 513):np.where(lambda_vals == 601)])  # Green bin
        lb=np.where(lambda_vals == 381E-9)
        lb=lb[0][0]
        rb=np.where(lambda_vals == 501E-9)
        rb=rb[0][0]
        image[j, k, 2] = np.sum(Spectrum[lb:rb])#np.sum(Spectrum[np.where(lambda_vals == 381):np.where(lambda_vals == 501)])  # Blue bin

# Adjust gain
image[:, :, 0] = image[:, :, 0]
image[:, :, 1] = image[:, :, 1] * 1.1
image[:, :, 2] = image[:, :, 2] * 4.2
image /= 15

plt.figure()
plt.imshow(image)
plt.show()



