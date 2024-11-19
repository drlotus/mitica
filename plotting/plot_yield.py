import numpy as np
import math
import matplotlib.pyplot as plt
import argparse


tolerance = 1e-6

pT_index = 1
phi_index = 2
y_index = 3
dNd3p_index = 5

def calculate_column_tolerance(column):
    """
    Calculate a dynamic tolerance for a given column based on the range of values in the column.
    Uses relative tolerance, but ensures a minimum tolerance.
    
    Parameters:
    - column: list of values for which to calculate tolerance
    - rel_tol: relative tolerance (fraction of the value range or typical scale)
    - min_tol: minimum tolerance to avoid too strict tolerance for very small values
    """
    if len(column) < 2:
        return 0.0  # If there's not enough data, no tolerance can be calculated

    column_range = max(column) - min(column)
    n = len(column)
    
    # Return tolerance as the average difference between values in the range
    return column_range / n

def get_tolerances(data):
    tolerances = {}
    columns = [y_index, pT_index, phi_index]  # Globally defined column indices
    for col in columns:
        col_values = [row[col] for row in data]
        tolerances[col] = calculate_column_tolerance(col_values)
    return tolerances

# Function to read the data from a file
def read_data(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.strip() or line.startswith("#"):
                continue
            columns = [float(value) for value in line.split()]
            data.append(columns)
    # Sort the data for easier integration
    data.sort(key=lambda x: (x[pT_index], x[y_index], x[phi_index]))
    return data

# Improved find_bins function with variable tolerance per column
def find_bins(column, tolerance):
    bins = []
    for value in column:
        if not any(math.isclose(value, bin_value, abs_tol=tolerance) for bin_value in bins):
            bins.append(value)
    return bins

# Group data by the bins in the given column
def group_rows_by_bins(data, bins, column, tolerance=1e-6):
    bin_to_rows = {bin_value: [] for bin_value in bins}
    for row in data:
        value_in_row = row[column]
        for bin_value in bins:
            if math.isclose(value_in_row, bin_value, abs_tol=tolerance):
                bin_to_rows[bin_value].append(row)
                break  # Once added to a bin, no need to check further bins
    return bin_to_rows

# Calculation of dN/dy with three tolerances
# def dNdy(data, tol_y=None, tol_pT=None, tol_phi=None):
def dNdy(data):
    y_column = [row[y_index] for row in data]
    y_bins = sorted(set(y_column))  # Use unique y values as bins
    
    # Group rows by exact y values
    grouped_data = {y_bin: [row for row in data if row[y_index] == y_bin] for y_bin in y_bins}

    dNdy_hist = []
    for y_bin in y_bins:
        data_in_y_bin = grouped_data[y_bin]
        if not data_in_y_bin:
            print("No data in y bin")
            continue
        
        # Extract unique pT values for this y bin and sort them
        pT_column = [row[pT_index] for row in data_in_y_bin]
        pT_bins = sorted(set(pT_column))  # Group by distinct pT values
        
        integral_over_pT = 0.0  # This will accumulate over pT
        for i in range(len(pT_bins) - 1):  # Loop over pT bins in pairs
            current_pT = pT_bins[i]
            next_pT = pT_bins[i + 1]
            
            # Extract data rows corresponding to the current pT bin
            data_in_current_pT_bin = [row for row in data_in_y_bin if math.isclose(row[pT_index], current_pT)]
            data_in_next_pT_bin = [row for row in data_in_y_bin if math.isclose(row[pT_index], next_pT)]
            
            if not data_in_current_pT_bin or not data_in_next_pT_bin:
                print(f"No data in pT bins: {current_pT}, {next_pT}")
                continue

            # Integrate over phi using trapezoidal rule
            integral_over_phi = 0.0
            for j in range(len(data_in_current_pT_bin) - 1):
                current_phi = data_in_current_pT_bin[j][phi_index]
                current_dNd3p = data_in_current_pT_bin[j][dNd3p_index]
                next_phi = data_in_current_pT_bin[j + 1][phi_index]
                next_dNd3p = data_in_current_pT_bin[j + 1][dNd3p_index]

                # Apply trapezoidal rule for phi integration
                integral_over_phi += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
            
            # Now integrate over pT using the result from phi integration
            pT_mid = 0.5 * (current_pT + next_pT)
            delta_pT = next_pT - current_pT
            integral_over_pT += pT_mid * delta_pT * integral_over_phi

        # Append the result for this y_bin (integral over pT)
        dNdy_hist.append([y_bin, integral_over_pT])

    return dNdy_hist

# Calculation of dN/pTdpT
def dNpTdpT(data):
    """
    Calculate dN/pTdpT by integrating over y and phi_p for each distinct pT bin.
    
    We are using:
    dN/pTdpT = ∫ dy ∫ dphi dN/d^3p
    """
    # Step 1: Extract distinct pT values and group data accordingly
    pT_column = [row[pT_index] for row in data]
    pT_bins = sorted(set(pT_column))  # Use unique pT values

    # Group data by pT values
    grouped_data = {pT_bin: [row for row in data if row[pT_index] == pT_bin] for pT_bin in pT_bins}
    
    dNpTdpT_hist = []
    
    # Step 2: Loop over each pT bin
    for pT_bin in pT_bins:
        data_in_bin = grouped_data[pT_bin]
        
        if not data_in_bin:  # Skip if no data in this pT bin
            continue
        
        integral = 0.0
        
        # Step 3: Integrate over y and phi
        y_column = [row[y_index] for row in data_in_bin]
        y_bins = sorted(set(y_column))  # Group by distinct y values
        
        for y_bin in y_bins:  # Loop over each y bin
            # Filter the data for this y_bin
            data_in_y_bin = [row for row in data_in_bin if row[y_index] == y_bin]
            
            if len(data_in_y_bin) < 2:  # Need at least two points for trapezoidal integration
                continue
            
            # Integrate over phi using the trapezoidal rule
            integral_over_phi = 0.0
            for i in range(len(data_in_y_bin) - 1):
                current_phi = data_in_y_bin[i][phi_index]
                current_dNd3p = data_in_y_bin[i][dNd3p_index]
                next_phi = data_in_y_bin[i + 1][phi_index]
                next_dNd3p = data_in_y_bin[i + 1][dNd3p_index]
                
                # Apply trapezoidal rule for phi integration
                integral_over_phi += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
            
            # Accumulate the result for this y_bin
            integral += integral_over_phi
        
        # Step 4: Store the result (integral for this pT bin)
        dNpTdpT_hist.append([pT_bin, integral])

    return dNpTdpT_hist

def N(data):
    """
    Integrate over all three variables: y, pT, and phi.
    """
    # Step 1: Group data by distinct pT values
    pT_column = [row[pT_index] for row in data]
    pT_bins = sorted(set(pT_column))  # Group by unique pT values
    
    # Step 2: Group data by distinct y values
    y_column = [row[y_index] for row in data]
    y_bins = sorted(set(y_column))  # Group by unique y values
    
    total_integral = 0.0
    
    # Step 3: Loop over pT and y bins
    for pT_bin in pT_bins:
        data_in_pT_bin = [row for row in data if row[pT_index] == pT_bin]
        
        for y_bin in y_bins:
            # Filter data for the current y bin
            data_in_y_bin = [row for row in data_in_pT_bin if row[y_index] == y_bin]
            
            if len(data_in_y_bin) < 2:
                continue  # Need at least two points for integration
            
            # Step 4: Perform integration over phi
            phi_integral = 0.0
            for i in range(len(data_in_y_bin) - 1):
                current_phi = data_in_y_bin[i][phi_index]
                current_dNd3p = data_in_y_bin[i][dNd3p_index]
                next_phi = data_in_y_bin[i + 1][phi_index]
                next_dNd3p = data_in_y_bin[i + 1][dNd3p_index]
                
                # Apply trapezoidal rule for phi integration
                phi_integral += 0.5 * (next_phi - current_phi) * (next_dNd3p + current_dNd3p)
            
            # Step 5: Accumulate results over pT bins using the result from phi integration
            total_integral += phi_integral * pT_bin  # Weight by pT for total integration
    
    return total_integral

# Flow harmonics v1 and v2
def flow(data, particle_number, harmonic):
    # Step 1: Extract distinct pT values and group data accordingly
    pT_column = [row[pT_index] for row in data]
    pT_bins = sorted(set(pT_column))  # Use unique pT values

    # Group data by pT values
    grouped_data = {pT_bin: [row for row in data if row[pT_index] == pT_bin] for pT_bin in pT_bins}
    
    weighted_flow_hist = []
    
    # Step 2: Loop over each pT bin
    for pT_bin in pT_bins:
        data_in_bin = grouped_data[pT_bin]
        
        if not data_in_bin:  # Skip if no data in this pT bin
            continue
        
        integral = 0.0
        
        # Step 3: Integrate over y and phi
        y_column = [row[y_index] for row in data_in_bin]
        y_bins = sorted(set(y_column))  # Group by distinct y values
        
        for y_bin in y_bins:  # Loop over each y bin
            # Filter the data for this y_bin
            data_in_y_bin = [row for row in data_in_bin if row[y_index] == y_bin]
            
            if len(data_in_y_bin) < 2:  # Need at least two points for trapezoidal integration
                continue
            
            # Integrate over phi using the trapezoidal rule with weighting
            integral_over_phi = 0.0
            for i in range(len(data_in_y_bin) - 1):
                current_phi = data_in_y_bin[i][phi_index]
                current_dNd3p = data_in_y_bin[i][dNd3p_index]
                next_phi = data_in_y_bin[i + 1][phi_index]
                next_dNd3p = data_in_y_bin[i + 1][dNd3p_index]
                
                # Weighting by cos(harmonic * φ)
                current_weight = math.cos(harmonic * current_phi)
                next_weight = math.cos(harmonic * next_phi)
                
                # Apply trapezoidal rule for phi integration with weighting
                integral_over_phi += 0.5 * (next_phi - current_phi) * (
                    current_weight * current_dNd3p + next_weight * next_dNd3p
                )
            
            # Accumulate the result for this y_bin
            integral += integral_over_phi / particle_number
        # Step 4: Store the result (integral for this pT bin)
        weighted_flow_hist.append([pT_bin, integral])

    return weighted_flow_hist

# Wrappers for flow harmonics v1 and v2
def v1(data, particle_number):
    return flow(data, particle_number, 1)

def v2(data, particle_number):
    return flow(data, particle_number, 2)


def nils_flow(data, particle_number, flow_order):
    """
    We are using:
    v1 = 1/N * int dy dphi cos(phi) dN/d3p
    v2 = 1/N * int dy dphi cos(2*phi) dN/d3p
    """
    if not (flow_order == 1 or flow_order == 2):
        raise ValueError("Floworder must be 1 or 2!")

    
    # Extract the p_T column from the data (1st column, index 0)
    pT_column = [row[pT_index] for row in data]
    pT_bins = find_bins(pT_column, 0)
    
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
        print(f"data_phi_integrated: {data_phi_integrated}")    
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



# Main program logic

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads output from Mitica and plots yield and flow')
    parser.add_argument('--file', type=str, required=True, help='Path to the Mitica output')
    args = parser.parse_args()

    # Read the data once
    file_path = args.file
    data = read_data(file_path)

    

    tolerances = get_tolerances(data)
    # Calculate dN/dy, dN/pTdpT, and N
    # dNdy_histogram = dNdy(data, tol_y=tolerances[y_index], tol_pT=tolerances[pT_index], tol_phi=tolerances[phi_index])
    dNdy_histogram = dNdy(data)
    dNpTdpT_histogram = dNpTdpT(data)
    n = N(data)
    print(f"N={n}")
    v1_pT_histogram = v1(data, n)

    # nils_flow(data, n, 2)

    v2_pT_histogram = v2(data, n)

    pT_min = 0.2
    pT_max = 3.0

    # Filter the dN/pTdpT_histogram to keep only pT values within the desired range
    dNpTdpT_histogram_cut = [entry for entry in dNpTdpT_histogram if pT_min <= entry[0] <= pT_max]
    
    pT_max = 2.0

    v1_pT_histogram_cut = [entry for entry in v1_pT_histogram if pT_min <= entry[0] <= pT_max]



    v2_pT_histogram_cut = [entry for entry in v2_pT_histogram if pT_min <= entry[0] <= pT_max]


    # Plot results
    x_values1, y_values1 = zip(*dNdy_histogram)
    x_values2, y_values2 = zip(*dNpTdpT_histogram)
    x_values3, y_values3 = zip(*v1_pT_histogram)
    x_values4, y_values4 = zip(*v2_pT_histogram)

    # Plot dN/dy spectrum
    plt.figure(1)
    plt.plot(x_values1, y_values1, marker='o', linestyle='-', color='steelblue')
    plt.title('dN/dy Spectrum')
    plt.xlabel('y')
    plt.ylabel('dN/dy')
    plt.savefig('dndy.png')

    # Plot dN/pTdpT spectrum
    plt.figure(2)
    plt.plot(x_values2, y_values2, marker='o', linestyle='-', color='seagreen')
    plt.yscale('log')
    plt.title('dN/pTdpT Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel(r'$dN/p_Tdp_T$')
    plt.savefig('pt.png')

    # Plot v1/pT spectrum
    plt.figure(3)
    plt.plot(x_values3, y_values3, marker='o', linestyle='-', color='palevioletred')
    plt.title(r'$v_1-p_T$ Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel('v1')
    plt.savefig('v1.png')

    # Plot v2/pT spectrum
    plt.figure(4)
    plt.plot(x_values4, y_values4, marker='o', linestyle='-', color='mediumvioletred')
    plt.title(r'$v_2-p_T$pT Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel(r'$v_2$')
    plt.savefig('v2.png')

    x_values5, y_values5 = zip(*dNpTdpT_histogram_cut)

    plt.figure(5)
    plt.plot(x_values5, y_values5, marker='o', linestyle='-', color='seagreen')
    plt.yscale('log')
    plt.title('dN/pTdpT Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel(r'$dN/p_Tdp_T$')
    plt.savefig('pt.png')

    x_values6, y_values6 = zip(*v1_pT_histogram_cut)

    # Plot v1/pT spectrum
    plt.figure(6)
    plt.plot(x_values6, y_values6, marker='o', linestyle='-', color='palevioletred')
    plt.title(r'$v_1-p_T$ Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel('v1')
    plt.savefig('v1.png')

    x_values7, y_values7 = zip(*v2_pT_histogram_cut)
    # Plot v2/pT spectrum
    plt.figure(7)
    plt.plot(x_values7, y_values7, marker='o', linestyle='-', color='mediumvioletred')
    plt.title(r'$v_2-p_T$pT Spectrum')
    plt.xlabel(r'$p_T$')
    plt.ylabel(r'$v_2$')
    plt.savefig('v2.png')



    plt.show()