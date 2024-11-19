#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:21:03 2024

@author: nils
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import trapezoid  # Instead of trapz
from scipy.interpolate import griddata
import argparse
import math

# Column indices for data parsing
mT_index = 0
pT_index = 1
phi_index = 2
y_index = 3
dNd3p_index = 5


def infer_rounding_precision(values):
    """
    Infers the appropriate number of decimal places based on the smallest difference
    between consecutive unique values in a list.
    
    Parameters:
    - values: List or array of numerical values.
    
    Returns:
    - The number of decimal places to use for rounding.
    """
    unique_values = sorted(set(values))
    
    if len(unique_values) < 2:
        return 0  # No precision needed if there's only one unique value
    
    # Find the smallest difference between consecutive unique values
    differences = np.diff(unique_values)
    min_difference = np.min(differences)
    
    # Calculate how many decimal places are needed to represent this difference
    if min_difference == 0:
        return 0  # No need to round if the values are already identical
    
    # Infer precision from the smallest difference
    precision = max(0, int(np.ceil(-np.log10(min_difference))))
    
    return precision

def read_data(file_path):
    """
    Reads the data from the file and returns all rows as a list of lists.
    Each row contains the full data line split into its respective columns.
    Automatically infers the rounding precision for y_p, pT, and phi_p columns.
    
    Parameters:
    - file_path: Path to the input file.
    """
    data = []
    y_values = []
    pT_values = []
    phi_p_values = []

    with open(file_path, 'r') as file:
        # Skip header and comment lines
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue
            
            # Split the line into columns and convert to floats
            columns = [float(value) for value in line.split()]
            
            if len(columns) == 6:  # Ensure 6 columns in each line
                data.append(columns)
                # Store y, pT, and phi_p values for precision inference
                pT_values.append(columns[pT_index])
                y_values.append(columns[y_index])
                phi_p_values.append(columns[phi_index])
            else:
                print(f"Skipping invalid line: {line.strip()}")

    # Infer rounding precision for y_p, pT, and phi_p columns
    y_round_digits = infer_rounding_precision(y_values)
    pT_round_digits = infer_rounding_precision(pT_values)
    phi_round_digits = infer_rounding_precision(phi_p_values)

    print(f"Inferred rounding precision: y_p: {y_round_digits}, pT: {pT_round_digits}, phi_p: {phi_round_digits}")
    
    # Round y_p, pT, and phi_p columns based on inferred precision
    for row in data:
        row[pT_index] = round(row[pT_index], pT_round_digits)
        row[y_index] = round(row[y_index], y_round_digits)
        row[phi_index] = round(row[phi_index], phi_round_digits)

    # Sort the data by pT, y_p, and phi_p (columns 1, 3, and 2 respectively)
    data.sort(key=lambda x: (x[pT_index], x[y_index], x[phi_index]))
    
    y_values = [row[y_index] for row in data]
    
    # Convert y_values to a set to get unique values
    unique_y_values = set(y_values)
    
    # Print the unique y values directly
    # print("Unique y values after rounding:", sorted(unique_y_values))

    return data
# def read_data(file_path):
#     """
#     Reads the data from the file and returns all rows as a list of lists.
#     Each row contains the full data line split into its respective columns.
#     Automatically infers the rounding precision for y_p, pT, and phi_p columns.
    
#     Parameters:
#     - file_path: Path to the input file.
#     """
#     data = []
#     y_values = []
#     pT_values = []
#     phi_p_values = []
#     phi_p = 0.0
#     pT = 0.0
#     y_p = 0.0
#     counter = 0
#     with open(file_path, 'r') as file:
#         # Skip header and comment lines
#         for line in file:
#             if not line.strip() or line.startswith("#"):
#                 continue
            
#             # Split the line into columns and convert to floats
#             columns = [float(value) for value in line.split()]
            
#             if len(columns) == 6:  # Ensure 6 columns in each line
#                 data.append(columns)
#                 # Store y, pT, and phi_p values for precision inference
#                 if counter > 0:
#                     if columns[pT_index] != pT:
#                         pT = columns[pT_index]
#                         pT_values.append(columns[pT_index])
#                     if columns[y_index] != y_p:
#                         y_p = columns[y_index]
#                         y_values.append(columns[y_index])
#                     if columns[phi_index] != phi_p:
#                         phi_p = columns[phi_index]
#                         phi_p_values.append(columns[y_index])
                    
#                 # pT_values.append(columns[pT_index])
#                 # y_values.append(columns[y_index])
#                 # phi_p_values.append(columns[phi_index])
#                 counter=counter+1
#             else:
#                 print(f"Skipping invalid line: {line.strip()}")

#     print("y values from the file:", y_values)
#     # Infer rounding precision for y_p, pT, and phi_p columns
#     y_round_digits = infer_rounding_precision(y_values)
#     pT_round_digits = infer_rounding_precision(pT_values)
#     phi_round_digits = infer_rounding_precision(phi_p_values)

#     print(f"Inferred rounding precision: y_p: {y_round_digits}, pT: {pT_round_digits}, phi_p: {phi_round_digits}")
    
#     # Round y_p, pT, and phi_p columns based on inferred precision
#     for row in data:
#         row[pT_index] = round(row[pT_index], pT_round_digits)
#         row[y_index] = round(row[y_index], y_round_digits)
#         row[phi_index] = round(row[phi_index], phi_round_digits)

#     # Sort the data by pT, y_p, and phi_p (columns 1, 3, and 2 respectively)
#     data.sort(key=lambda x: (x[pT_index], x[y_index], x[phi_index]))

#     y_values = [row[y_index] for row in data]
    
#     # Convert y_values to a set to get unique values
#     unique_y_values = set(y_values)
    
#     # Print the unique y values directly
#     print("Unique y values after rounding:", sorted(unique_y_values))

#     return data
# def check_unique_y_values(data):
#     """
#     Extract and print the unique y values as they are, to check if the data
#     is being processed or read correctly.
#     """
    # y_values = [row[y_index] for row in data]
    
    # # Convert y_values to a set to get unique values
    # unique_y_values = set(y_values)
    
    # # Print the unique y values directly
    # print("Unique y values found in the file:", unique_y_values)

#     # If you want, return the sorted list of unique values
#     return sorted(unique_y_values)


def interpolate_dN3p(data):
    """
    Interpolates dN/d^3p in terms of p_T, y, and phi_p using the actual unique values.
    """

    y_values_2 = [row[y_index] for row in data]
    
    # Convert y_values to a set to get unique values
    
    # Print the unique y values directly
    # print("Unique y values again:", sorted(unique_y_values_2))

    # Extract relevant columns
    pT = np.array([row[pT_index] for row in data])
    y_p = np.array([row[y_index] for row in data])
    phi_p = np.array([row[phi_index] for row in data])
    dNd3p = np.array([row[dNd3p_index] for row in data])


    print("dN/d^3p values before interpolation:\n", dNd3p[:10])  # Print the first 10 values

    # Find the unique points for each axis
    unique_pT = np.unique(pT)
    unique_y = np.unique(y_p)
    unique_phi_p = np.unique(phi_p)

    # print(f"Unique pT values: {unique_pT}")
    # print(f"Unique y values: {unique_y}")
    # print(f"Unique phi_p values: {unique_phi_p}")

    # Use the actual unique points for grid creation
    pT_grid, y_grid, phi_p_grid = np.meshgrid(unique_pT, unique_y, unique_phi_p)

    # Perform the interpolation using griddata
    interpolated_dN3p = griddata(
        points=(pT, y_p, phi_p),
        values=dNd3p,
        xi=(pT_grid, y_grid, phi_p_grid),
        method='linear',  # 'linear', 'nearest', or 'cubic' can be chosen based on your needs
        fill_value=0  # Optional: to handle extrapolation if needed
    )

    # print("Interpolated dN/d^3p values:\n", interpolated_dN3p)

    print("Shape of interpolated pT_grid:", pT_grid.shape)
    print("Shape of interpolated y_grid:", y_grid.shape)
    print("Shape of interpolated phi_p_grid:", phi_p_grid.shape)
    print("Shape of interpolated dN/d^3p:", interpolated_dN3p.shape)
    for idx, label in [(0, 'y'), (1, 'phi'), (2, 'pT')]:
        plt.figure()
    
    if idx == 0:
        plt.imshow(interpolated_dN3p[:, :, 0], extent=(pT_grid.min(), pT_grid.max(), phi_p_grid.min(), phi_p_grid.max()))
        plt.colorbar(label=f"dN/d^3p at {label} = {y_grid[0, 0, 0]}")
        plt.title(f"Interpolated dN/d^3p slice for fixed {label}")
        plt.xlabel(f"{label}")
        plt.ylabel("pT")
    elif idx == 1:
        plt.imshow(interpolated_dN3p[:, 0, :], extent=(pT_grid.min(), pT_grid.max(), y_grid.min(), y_grid.max()))
        plt.colorbar(label=f"dN/d^3p at {label} = {phi_p_grid[0, 0, 0]}")
        plt.title(f"Interpolated dN/d^3p slice for fixed {label}")
        plt.xlabel(f"{label}")
        plt.ylabel("pT")
    elif idx == 2:
        plt.imshow(interpolated_dN3p[0, :, :], extent=(y_grid.min(), y_grid.max(), phi_p_grid.min(), phi_p_grid.max()))
        plt.colorbar(label=f"dN/d^3p at {label} = {pT_grid[0, 0, 0]}")
        plt.title(f"Interpolated dN/d^3p slice for fixed {label}")
        plt.xlabel(f"{label}")
        plt.ylabel("phi")
        plt.show()


    return pT_grid, y_grid, phi_p_grid, interpolated_dN3p

def dNdy(data):
    """
    Calculates dN/dy using interpolation and numerical integration.
    dN/dy = int pT dpT dphi dN/d^3p
    """
    # First, perform interpolation over the full dataset
    pT_grid, _,phi_p_grid, interpolated_dN3p = interpolate_dN3p(data)
    
    # Get unique y values from the data (this should be 1 in your case)
    y_values = sorted(set(row[y_index] for row in data))
    
    dNdy_hist = []

    for y_value in y_values:
        # Filter the interpolated data for this specific y_value
        # Since y is fixed, we will take the relevant slice of the interpolated data
        # However, for a fixed y, we just append the interpolated integral for all y
        integral_pT_phi = trapezoid(trapezoid(interpolated_dN3p, x=phi_p_grid, axis=0), x=pT_grid[:, 0])
        
        # Append the result for this y value
        dNdy_hist.append([y_value, integral_pT_phi])

    return dNdy_hist

def dNpTdpT(data):
    """
    Calculates dN/pTdpT using interpolation and numerical integration.
    dN/pTdpT = int dy dphi dN/d^3p
    """
    # Perform interpolation over the full dataset
    pT_grid, y_grid, phi_p_grid, interpolated_dN3p = interpolate_dN3p(data)
    
    # Get unique pT values from the data
    pT_values = sorted(set(row[pT_index] for row in data))
    
    dNpTdpT_hist = []

    for pT_value in pT_values:
        # Filter the interpolated data for this specific pT_value
        # Find the index for the specific pT value in the grid
        pT_idx = np.where(np.isclose(pT_grid[:, 0, 0], pT_value))[0][0]  # Find the index for the specific pT value
        
        # Numerically integrate over phi_p (azimuthal angle) and y (rapidity)
        # First integrate over phi_p
        integral_phi = trapezoid(interpolated_dN3p[pT_idx, :, :], x=phi_p_grid[:, 0])

        # Now integrate over y
        integral_y_phi = trapezoid(integral_phi, x=y_grid[:, 0, 0])

        # Append the result for this pT value (no division by pT)
        dNpTdpT_hist.append([pT_value, integral_y_phi])

    return dNpTdpT_hist

def N(data):
    """
    Calculates N using interpolation and numerical integration.
    N = int pT dpT dy dphi dN/d^3p
    """

    pT_grid, y_grid, phi_p_grid, interpolated_dN3p = interpolate_dN3p(data)
    integral_phi = trapezoid(interpolated_dN3p, x=phi_p_grid, axis=2)
    integral_phi_y = trapezoid(integral_phi, x=y_grid[:, 0, 0], axis=1)
    integral_pT_y_phi = trapezoid(integral_phi_y, x=pT_grid[:, 0, 0])
    return integral_pT_y_phi

def flow(data, particle_number, flow_order):
    """
    Calculates the flow harmonics v1 or v2 using interpolation and numerical integration.
    v1 = 1/N * int dy dphi cos(phi) dN/d^3p
    v2 = 1/N * int dy dphi cos(2*phi) dN/d^3p
    """
    if flow_order not in (1, 2):
        raise ValueError("Flow order must be 1 or 2!")

    pT_grid, y_grid, phi_p_grid, interpolated_dN3p = interpolate_dN3p(data)

    if flow_order == 1:
        flow_factor = np.cos(phi_p_grid)
    elif flow_order == 2:
        flow_factor = np.cos(2 * phi_p_grid)

    dN3p_flow = interpolated_dN3p * flow_factor
    integral_phi = trapezoid(dN3p_flow, x=phi_p_grid, axis=2)
    integral_phi_y = trapezoid(integral_phi, x=y_grid[:, 0, 0], axis=1)
    integral_pT_y_phi = trapezoid(integral_phi_y, x=pT_grid[:, 0, 0])
    
    flow_harmonic = integral_pT_y_phi / particle_number
    return flow_harmonic

def v1(data, particle_number):
    return flow(data, particle_number, 1)

def v2(data, particle_number):
    return flow(data, particle_number, 2)

def plot_original_dNd3p(data, fixed_y=None, fixed_pT=None, fixed_phi=None):
    """
    Plots the original dN/d^3p values for a fixed y, pT, or phi.
    
    Parameters:
    - data: List of rows containing the parsed data.
    - fixed_y: (Optional) Fixed y value to slice the data by.
    - fixed_pT: (Optional) Fixed pT value to slice the data by.
    - fixed_phi: (Optional) Fixed phi value to slice the data by.
    """
    pT_values = np.array([row[pT_index] for row in data])
    y_values = np.array([row[y_index] for row in data])
    phi_values = np.array([row[phi_index] for row in data])
    dNd3p_values = np.array([row[dNd3p_index] for row in data])

    # Plot for fixed y
    if fixed_y is not None:
        mask = np.isclose(y_values, fixed_y)
        pT_filtered = pT_values[mask]
        phi_filtered = phi_values[mask]
        dNd3p_filtered = dNd3p_values[mask]

        plt.figure()
        plt.scatter(pT_filtered, dNd3p_filtered, c=phi_filtered, cmap='viridis', label=f"Fixed y = {fixed_y}")
        plt.colorbar(label="phi")
        plt.xlabel(r"$p_T$")
        plt.ylabel(r"$dN/d^3p$")
        plt.title(f"Original dN/d^3p at y = {fixed_y}")
        plt.show()

    # Plot for fixed pT
    if fixed_pT is not None:
        mask = np.isclose(pT_values, fixed_pT)
        y_filtered = y_values[mask]
        phi_filtered = phi_values[mask]
        dNd3p_filtered = dNd3p_values[mask]

        plt.figure()
        plt.scatter(y_filtered, dNd3p_filtered, c=phi_filtered, cmap='plasma', label=f"Fixed pT = {fixed_pT}")
        plt.colorbar(label="phi")
        plt.xlabel("y")
        plt.ylabel(r"$dN/d^3p$")
        plt.title(f"Original dN/d^3p at pT = {fixed_pT}")
        plt.show()

    # Plot for fixed phi
    if fixed_phi is not None:
        mask = np.isclose(phi_values, fixed_phi)
        pT_filtered = pT_values[mask]
        y_filtered = y_values[mask]
        dNd3p_filtered = dNd3p_values[mask]

        plt.figure()
        plt.scatter(pT_filtered, dNd3p_filtered, c=y_filtered, cmap='cividis', label=f"Fixed phi = {fixed_phi}")
        plt.colorbar(label="y")
        plt.xlabel(r"$p_T$")
        plt.ylabel(r"$dN/d^3p$")
        plt.title(f"Original dN/d^3p at phi = {fixed_phi}")
        plt.show()

#######################
# USAGE
#######################

parser = argparse.ArgumentParser(description='Reads output from Mitica and plots yield and flow')
parser.add_argument('--file', type=str, required=True, help='Path to the Mitica output')
args = parser.parse_args()

# Step 2: Read the file and process data
file_path = args.file
data = read_data(file_path)


# plot_original_dNd3p(data, fixed_y=-0.0)  # You can change this to any specific y value
# plot_original_dNd3p(data, fixed_pT=0.0)  # You can change this to any specific pT value
# plot_original_dNd3p(data, fixed_phi=0.0)  # You can change this to any specific phi value

# Compute dN/dy, dN/pTdpT, and N
n = N(data)
print(f"Total N: {n}")

# dNdy_histogram = dNdy(data)
# dNpTdpT_histogram = dNpTdpT(data)

# # Compute flow harmonics
# dv1dpT_histogram = v1(data, n)
# dv2dpT_histogram = v2(data, n)

# # Plot results
# x_values1, y_values1 = zip(*dNdy_histogram)
# x_values2, y_values2 = zip(*dNpTdpT_histogram)
# x_values3, y_values3 = zip(*dv1dpT_histogram)
# x_values4, y_values4 = zip(*dv2dpT_histogram)

# # Plot the dN/dy spectrum
# plt.figure(1)
# plt.plot(x_values1, y_values1, marker='o', linestyle='-', color='steelblue')
# plt.title('dN/dy Spectrum')
# plt.xlabel('y')
# plt.ylabel('dN/dy')

# # Plot the dN/pTdpT spectrum
# plt.figure(2)
# plt.plot(x_values2, y_values2, marker='o', linestyle='-', color='seagreen')
# plt.yscale('log')
# plt.title(r'$dN/p_Tdp_T$ Spectrum (log)')
# plt.xlabel(r'$p_T$')
# plt.ylabel(r'$dN/p_Tdp_T$')

# # Plot the dv1/dpT spectrum
# plt.figure(3)
# plt.plot(x_values3, y_values3, marker='o', linestyle='-', color='palevioletred')
# plt.title('dv1/dpT Spectrum')
# plt.xlabel('pT')
# plt.ylabel('dv1/dpT')

# # Plot the dv2/dpT spectrum
# plt.figure(4)
# plt.plot(x_values4, y_values4, marker='o', linestyle='-', color='mediumvioletred')
# plt.title('dv2/dpT Spectrum')
# plt.xlabel('pT')
# plt.ylabel('dv2/dpT')

# plt.show()