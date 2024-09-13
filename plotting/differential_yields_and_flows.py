#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:21:03 2024

@author: nils
"""
import numpy as np
import math
import matplotlib.pyplot as plt


tolerance=1e-6

pT_index = 1
phi_index = 2
y_index = 3
dNd3p_index = 5

def read_data(file_path):
    """
    Reads the data from the file and returns all rows as a list of lists.
    Each row contains the full data line split into its respective columns.
    """
    data = []
    
    with open(file_path, 'r') as file:
        # Skip header and comment lines
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue
            
            # Split the line into columns and convert to floats
            columns = [float(value) for value in line.split()]
            data.append(columns)
            
    # Sort the data by pT, y_p, and phi_p (columns 1, 3, and 2 respectively)
    data.sort(key=lambda x: (x[pT_index], x[y_index], x[phi_index]))
    
    return data


def find_bins(column, tolerance=1e-6):
    """
    Finds unique bins in the y_p column, considering a numerical precision tolerance.
    """
    bins = []
    
    for value in column:
        # Check if the value is close to any existing bin
        if not any(math.isclose(value, bin_value, abs_tol=tolerance) for bin_value in bins):
            bins.append(value)
    
    return bins


def group_rows_by_bins(data, bins, column, tolerance=1e-6):
    """
    Groups the rows of data by the bins in the corresponding column, 
    considering a numerical precision tolerance.
    """
    bin_to_rows = {bin_value: [] for bin_value in bins}
    counter = 0
    for row in data:
        counter += 1
        value_in_row = row[column]
        for bin_value in bins:
            if math.isclose(value_in_row, bin_value, abs_tol=tolerance):
                bin_to_rows[bin_value].append(row)
                break  # Once added to a bin, no need to check further bins
    return bin_to_rows


def dNdy(file_path):
    """
    We are using:
    dN/dy = int pTdpT dphi dN/d3p
    """
    data = read_data(file_path)
    
    # Extract the y_p column from the data (3rd column, index 2)
    y_column = [row[y_index] for row in data]
    y_bins = find_bins(y_column)
    
    # Group rows by y_p bins
    grouped_data = group_rows_by_bins(data, y_bins, y_index)
    
    
    dNdy_hist = []
    
    
    for y_bin in y_bins:
        data_in_bin = grouped_data[y_bin]
        
        if not data_in_bin:  # Skip if no data in this bin
            continue
        
        phi_start = data_in_bin[0][phi_index]
        integral = 0.0
        dy = y_bins[1] - y_bins[0]
        
        
        # Stores (pT, phi_integral) for the current y
        data_phi_integrated = []
        
        # Integration over phi_p
        for i in range(len(data_in_bin)-1):
            current_phi = data_in_bin[i][phi_index]
            current_dNd3p = data_in_bin[i][dNd3p_index]
            next_phi = data_in_bin[i+1][phi_index]
            next_dNd3p = data_in_bin[i+1][dNd3p_index]
            
            if not math.isclose(phi_start, next_phi, abs_tol = tolerance):
                # Apply trapezoidal rule for integration
                integral += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
            
            else:
                data_phi_integrated.append([data_in_bin[i][pT_index], integral])
                integral = 0
                
            if i == len(data_in_bin)-1:
                data_phi_integrated.append([data_in_bin[i][pT_index], integral])
            
        integral = 0
        
        # Integration over pT
        for i in range(len(data_phi_integrated)-1):
            current_pT = data_phi_integrated[i][0]
            current_integrand = data_phi_integrated[i][1]
            next_pT = data_phi_integrated[i+1][0]
            next_integrand = data_phi_integrated[i+1][1]
            
            # Apply trapezoidal rule for integration
            integral += current_pT * 0.5 * (next_pT - current_pT) * (next_integrand + current_integrand)
             
        dNdy_hist.append([y_bin, integral/dy])

    return dNdy_hist


def dNpTdpT(file_path): 
    """
    We are using:
    dN/pTdpT = int dy dphi dN/d3p
    """
    data = read_data(file_path)
    
    # Extract the y_p column from the data (3rd column, index 2)
    pT_column = [row[pT_index] for row in data]
    pT_bins = find_bins(pT_column)
    dpT = pT_bins[1] - pT_bins[0]
    
    # Group rows by p_T bins
    grouped_data = group_rows_by_bins(data, pT_bins, pT_index)
    
    
    dNpTdpT_hist = []
    
    
    for pT_bin in pT_bins:
        data_in_bin = grouped_data[pT_bin]
        
        if not data_in_bin:  # Skip if no data in this bin
            continue
        
        phi_start = data_in_bin[0][phi_index]
        integral = 0.0
        
        # Stores (y, phi_integral) for the current pT
        data_phi_integrated = []
        
        # for pT=0 set the result to nan as it diverges
        if math.isclose(pT_bin, 0.0, abs_tol=tolerance):
            dNpTdpT_hist.append([0.0, np.nan])
        else:
            # Integration over phi_p
            for i in range(len(data_in_bin)-1):
                current_phi = data_in_bin[i][phi_index]
                current_dNd3p = data_in_bin[i][dNd3p_index]
                next_phi = data_in_bin[i+1][phi_index]
                next_dNd3p = data_in_bin[i+1][dNd3p_index]
                
                if not math.isclose(phi_start, next_phi, abs_tol = tolerance):
                    # Apply trapezoidal rule for integration
                    integral += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
                
                else:
                    data_phi_integrated.append([data_in_bin[i][y_index], integral])
                    integral = 0
                    
                if i == len(data_in_bin)-1:
                    data_phi_integrated.append([data_in_bin[i][y_index], integral])
                
            integral = 0
            
            # Integration over y
            for i in range(len(data_phi_integrated)-1):
                current_y = data_phi_integrated[i][0]
                current_integrand = data_phi_integrated[i][1]
                next_y = data_phi_integrated[i+1][0]
                next_integrand = data_phi_integrated[i+1][1]
                
                # Apply trapezoidal rule for integration
                integral += 0.5 * (next_y - current_y) * (next_integrand + current_integrand)
                 
            dNpTdpT_hist.append([pT_bin, integral/(dpT * pT_bin)])

    return dNpTdpT_hist


def N(file_path):
    """
    We are using:
    dN/pTdpT = int pTdpT dy dphi dN/d3p
    """
    data = read_data(file_path)
    
    phi_start = data[0][phi_index]
    integral = 0.0
    
    
    # Stores (pT, y, phi_integral)
    data_phi_integrated = []
    
    # Integration over phi_p
    for i in range(len(data)-1):
        current_phi = data[i][phi_index]
        current_dNd3p = data[i][dNd3p_index]
        next_phi = data[i+1][phi_index]
        next_dNd3p = data[i+1][dNd3p_index]
        
        if not math.isclose(phi_start, next_phi, abs_tol = tolerance):
            # Apply trapezoidal rule for integration
            integral += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
        
        else:
            data_phi_integrated.append([data[i][pT_index], data[i][y_index], integral])
            integral = 0
            
        if i == len(data)-1:
            data_phi_integrated.append([data[i][pT_index], data[i][y_index], integral])
        
    integral = 0
    
    
    # Stores (pT, phi_y_integral)
    data_phi_y_integrated = []
    
    y_start = data_phi_integrated[0][1]
    # Integration over y
    for i in range(len(data_phi_integrated)-1):
        current_y = data_phi_integrated[i][1]
        current_integrand = data_phi_integrated[i][2]
        next_y = data_phi_integrated[i+1][1]
        next_integrand = data_phi_integrated[i+1][2]
        
        if not math.isclose(y_start, next_y, abs_tol = tolerance):
            # Apply trapezoidal rule for integration
            integral += 0.5 * (next_y - current_y) * (next_integrand + current_integrand)
        
        else:
            data_phi_y_integrated.append([data_phi_integrated[i][0], integral])
            integral = 0
            
        if i == len(data_phi_integrated)-1:
            data_phi_y_integrated.append([data_phi_integrated[i][0], integral])
        
    integral = 0
    del data_phi_integrated

    
    # Integration over pT
    for i in range(len(data_phi_y_integrated)-1):
        current_pT = data_phi_y_integrated[i][0]
        current_integrand = data_phi_y_integrated[i][1]
        next_pT = data_phi_y_integrated[i+1][0]
        next_integrand = data_phi_y_integrated[i+1][1]
        
        # Apply trapezoidal rule for integration
        integral += current_pT * 0.5 * (next_pT - current_pT) * (next_integrand + current_integrand)      

    return integral
        

def flow(file_path, particle_number, flow_order):
    """
    We are using:
    v1 = 1/N * int dy dphi cos(phi) dN/d3p
    v2 = 1/N * int dy dphi cos(2*phi) dN/d3p
    """
    if not (flow_order == 1 or flow_order == 2):
        raise ValueError("Floworder must be 1 or 2!")
    
    data = read_data(file_path)
    
    # Extract the p_T column from the data (1st column, index 0)
    pT_column = [row[pT_index] for row in data]
    pT_bins = find_bins(pT_column)
    
    # Group rows by p_T bins
    grouped_data = group_rows_by_bins(data, pT_bins, pT_index)
    
    dv1dpT_hist = []
    
    
    for pT_bin in pT_bins:
        data_in_bin = grouped_data[pT_bin]
        
        if not data_in_bin:  # Skip if no data in this bin
            continue
        
        phi_start = data_in_bin[0][phi_index]
        integral = 0.0
        
        # Stores (pT, y, phi_integral)
        data_phi_integrated = []
        
        # Integration over phi_p
        for i in range(len(data_in_bin)-1):
            current_phi = data_in_bin[i][phi_index]
            current_dNd3p = data_in_bin[i][dNd3p_index]
            next_phi = data_in_bin[i+1][phi_index]
            next_dNd3p = data_in_bin[i+1][dNd3p_index]
            
            if not math.isclose(phi_start, next_phi, abs_tol = tolerance):
                # Apply trapezoidal rule for integration
                if flow_order == 1:
                    integral += 0.5 * (next_phi - current_phi) * (next_dNd3p * math.cos(next_phi) + current_dNd3p * math.cos(current_phi))
                elif flow_order == 2:
                    integral += 0.5 * (next_phi - current_phi) * (next_dNd3p * math.cos(2 * next_phi) + current_dNd3p * math.cos(2 * current_phi))
            
            else:
                data_phi_integrated.append([data_in_bin[i][pT_index], data_in_bin[i][y_index], integral])
                integral = 0
                
            if i == len(data_in_bin)-1:
                data_phi_integrated.append([data_in_bin[i][pT_index], data_in_bin[i][y_index], integral])
            
        integral = 0
        
        data_phi_y_integrated = []
        
        y_start = data_phi_integrated[0][1]
        
        # Integration over y
        for i in range(len(data_phi_integrated)-1):
            current_y = data_phi_integrated[i][1]
            current_integrand = data_phi_integrated[i][2]
            next_y = data_phi_integrated[i+1][1]
            next_integrand = data_phi_integrated[i+1][2]
            
            if not math.isclose(y_start, next_y, abs_tol = tolerance):
                # Apply trapezoidal rule for integration
                integral += 0.5 * (next_y - current_y) * (next_integrand + current_integrand)
                
            else:
                data_phi_y_integrated.append([data_phi_integrated[i][0], integral])
                integral = 0
             
            if i == len(data_phi_integrated)-1:
                data_phi_y_integrated.append([data_phi_integrated[i][0], integral])
                
        integral = 0
        del data_phi_integrated
        
                
        dv1dpT_hist.append([pT_bin, integral/particle_number])

    return dv1dpT_hist


def v1(file_path, particle_number):
    return flow(file_path, particle_number, 1)


def v2(file_path, particle_number):
    return flow(file_path, particle_number, 2)
    


#######################
# USAGE
#######################

file_path = '/Users/nils/Desktop/Projects/Polarization/Masoud_Integration_script/y_test_for_nils.dat'

dNdy_histogram = dNdy(file_path)
dNpTdpT_histogram = dNpTdpT(file_path)
n = N(file_path)
dv1dpT_histogram = v1(file_path, n)
dv2dpT_histogram = v2(file_path, n)



# Plot
x_values1, y_values1 = zip(*dNdy_histogram)
x_values2, y_values2 = zip(*dNpTdpT_histogram)
x_values3, y_values3 = zip(*dv1dpT_histogram)
x_values4, y_values4 = zip(*dv2dpT_histogram)

# Plot the dN/dy spectrum
plt.figure(1)
plt.plot(x_values1, y_values1, marker='o', linestyle='-', color='steelblue')
plt.title('dN/dy Spectrum')
plt.xlabel('y')
plt.ylabel('dN/dy')

# Plot the dN/pTdpT spectrum
plt.figure(2)
plt.plot(x_values2, y_values2, marker='o', linestyle='-', color='seagreen')
plt.yscale('log')
plt.title('dN/pTdpT Spectrum')
plt.xlabel('ptdpT')
plt.ylabel('dN/pTdpT')

# Plot the dv1/dpT spectrum
plt.figure(3)
plt.plot(x_values3, y_values3, marker='o', linestyle='-', color='palevioletred')
plt.title('dv1/dpT Spectrum')
plt.xlabel('pT')
plt.ylabel('dv1/dpT')

# Plot the dv2/dpT spectrum
plt.figure(4)
plt.plot(x_values4, y_values4, marker='o', linestyle='-', color='mediumvioletred')
plt.title('dv2/dpT Spectrum')
plt.xlabel('pT')
plt.ylabel('dv2/dpT')

plt.show()
