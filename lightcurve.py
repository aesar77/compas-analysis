import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import astropy.constants as ap
import astropy.units as un 

# t = np.linspace(5,35,1000) * un.hour
# t = t.to(un.s)
# M = (8.5 * ap.M_sun) 
# m = (0.35 * ap.M_sun)
# phi = 0.33 * un.arcsec
# a = 17 * un.au #(au)
# # a = a.to(un.m)
# tE = 2 * (a/ap.c) * np.sqrt(M/((M + m)*un.arcsec)) 
# print((tE).unit)
# a = a.to(un.m)
# Rs = 2 * ap.G * M / ap.c**2 
# print(ap.c.unit)
# tzero = 20 * un.hour
# tzero = tzero.to(un.s)
# print(Rs.unit)
# uzero = np.sqrt(a/(2*Rs) * phi)
# print((uzero**2).unit)
# print((((t - tzero)/tE)**2).unit)
# u_all= np.sqrt((uzero**2 +  ((t - tzero)/tE)**2)/un.arcsec) /3600
# print(u_all.unit)
# print(u_all[0])

t = np.linspace(5,35,1000) * un.hour
M = (8.5 * ap.M_sun) 
m = (0.35 * ap.M_sun)
phi = 0.33
a = 17 * un.au #(au)
tE = 2 * (a/ap.c) * np.sqrt(M/((M + m))) 
print((tE).unit)
Rs = 2 * ap.G * M / ap.c**2 
print(ap.c.unit)
tzero = 20 * un.hour
print(Rs.unit)
uzero = np.sqrt(a/(2*Rs) * phi)
print((uzero**2).unit)
print((((t - tzero)/tE)**2).unit)
u_all= np.sqrt((uzero**2 +  ((t - tzero)/tE)**2))
print(u_all.unit)
print(u_all[0])

# Define the parameters
rho_star = 0.81  # Example value for angular size of the source star (relative to Einstein radius)
import numpy as np
from scipy.special import ellipe, ellipk, ellipkinc, ellipeinc
A = []
# Step 1: Compute k and n
for u in u_all:
    print(u)
    u = u.value
    k = np.sqrt(4 / (4 + (u - rho_star)**2))  # Elliptic modulus
    n = 4 * rho_star / (u + rho_star)**2      # Parameter for Pi(n; k)

    # Step 2: Compute the elliptic integrals
    E_k = ellipe(k)                          # E(k): Elliptic integral of the second kind
    K_k = ellipk(k)                          # K(k): Elliptic integral of the first kind

    # Step 3: Compute each term in A(u, rho_star)
    term_10 = (1 / (2 * np.pi)) * ((u + rho_star) / rho_star**2) * np.sqrt(4 + (u - rho_star)**2) * E_k
    term_11 = (- (u - rho_star) / rho_star**2) * (8 + u**2 - rho_star**2) / np.sqrt(4 + (u - rho_star)**2) * K_k

    # For the third term, scipy doesn't directly have elliptic Pi(n; k).
    # So, we approximate it (or ignore it if we need to simplify).
    term_12 = 0  # Placeholder if Pi(n; k) is not included or approximated.

    # Total magnification A(u, rho_star)
    A_total = term_10 + term_11 + term_12

    # Output the result
    print(f"Magnification A(u, rho_star): {A_total}")
    A.append(A_total)

# Output the result
print(f"Magnification A(u, rho_star): {A_total}")

# A = 2/np.pi * ( (1 + 1/rho**2) * np.arcsin(1/np.sqrt(1 + 1/rho**2)) + 1/rho)

plt.plot(t, A)
plt.show()