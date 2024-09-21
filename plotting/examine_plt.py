import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata


# Load the file
file_path = './output/u_dot_n.dat'  # replace with your file path
columns = ['tau', 'x', 'y', 'eta', 'u.n', '|theta|', '|sigma|', '|tvort|', '|tshear|', '|acc|']
data = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)

t_res = 100

# Step 1: Define bins for tau
tau_bins = np.linspace(data['tau'].min(), data['tau'].max(), t_res)  # 10 bins for tau
data['tau_bin'] = pd.cut(data['tau'], bins=tau_bins, labels=False)  # Assign tau to bins

# Step 2: Group the data by tau bins and calculate the mean for each field
grouped_data = data.groupby('tau_bin').mean()

# Step 3: Calculate the center of each tau bin for plotting
tau_bin_centers = 0.5 * (tau_bins[:-1] + tau_bins[1:])

# Step 4: Interpolate the mean values of the fields as a function of tau_bin_centers

f_un  = interp1d(tau_bin_centers, grouped_data['u.n'], kind='cubic', fill_value="extrapolate")

# Step 5: Generate a finer set of tau values for interpolation
tau_fine = np.linspace(tau_bin_centers.min(), tau_bin_centers.max(), 100)

# Step 6: Plot the interpolated data
plt.figure(figsize=(10, 6))

# # Plot the original binned data
# plt.scatter(tau_bin_centers, grouped_data['|theta|'], color='red', label='Binned data')

# Plot the interpolated data
plt.plot(tau_fine, f_un(tau_fine), label='Interpolated u.n', color='blue')

plt.xlabel(r'$\tau$ ', fontsize=14)
plt.ylabel(r'$\langle u.n \rangle$ ', fontsize=14)
plt.legend()
plt.show()


f_theta = interp1d(tau_bin_centers, grouped_data['|theta|'], kind='cubic', fill_value="extrapolate")
f_sigma = interp1d(tau_bin_centers, grouped_data['|sigma|'], kind='cubic', fill_value="extrapolate")
f_acc= interp1d(tau_bin_centers, grouped_data['|acc|'], kind='cubic', fill_value="extrapolate")


# Step 6: Plot the interpolated data
plt.figure(figsize=(10, 6))


f_theta = interp1d(tau_bin_centers, grouped_data['|theta|'], kind='cubic', fill_value="extrapolate")
f_sigma = interp1d(tau_bin_centers, grouped_data['|sigma|'], kind='cubic', fill_value="extrapolate")
f_acc= interp1d(tau_bin_centers, grouped_data['|acc|'], kind='cubic', fill_value="extrapolate")


# Plot the interpolated data
plt.plot(tau_fine, f_theta(tau_fine), label=r'Interpolated $\theta$', color='blue')
plt.plot(tau_fine, f_sigma(tau_fine), label=r'Interpolated |$\sigma$|', color='red')
plt.plot(tau_fine, f_acc(tau_fine), label=r'Interpolated |a|', color='green')


plt.xlabel(r'$\tau$ (fm) ', fontsize=14)
plt.ylabel(r'$\partial u$ (GeV) ', fontsize=14)
plt.title(r'Interpolated u derivatives as functions of $\tau$', fontsize=16)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))


f_tvort = interp1d(tau_bin_centers, grouped_data['|tvort|'], kind='cubic', fill_value="extrapolate")
f_tshear = interp1d(tau_bin_centers, grouped_data['|tshear|'], kind='cubic', fill_value="extrapolate")


# Plot the interpolated data
plt.plot(tau_fine, f_tvort(tau_fine), label=r'Interpolated $\varpi$', color='blue')
plt.plot(tau_fine, f_tshear(tau_fine), label=r'Interpolated $\xi', color='red')


plt.xlabel(r'$\tau$ (fm) ', fontsize=14)
plt.ylabel(r'$\partial \beta$  ', fontsize=14)
plt.title(r'Interpolated $\beta$ derivatives as functions of $\tau$', fontsize=16)
plt.legend()
plt.show()


data = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)

# Step 1: Define bins for (x, y)
x_bins = np.linspace(data['x'].min(), data['x'].max(), 10)  # 10 bins for x
y_bins = np.linspace(data['y'].min(), data['y'].max(), 10)  # 10 bins for y



# Bin the data in x and y
data['x_bin'] = pd.cut(data['x'], bins=x_bins, labels=False)
data['y_bin'] = pd.cut(data['y'], bins=y_bins, labels=False)

# Step 2: Group by (tau, eta) and find the mean of each field
grouped_data = data.groupby(['tau', 'eta']).mean()

# Step 3: Interpolate the mean values as a function of (tau, eta)
# Use griddata to interpolate |theta| over the grid of (tau, eta)
tau_unique = np.linspace(data['tau'].min(), data['tau'].max(), 100)  # Fine tau grid
eta_unique = np.linspace(data['eta'].min(), data['eta'].max(), 100)  # Fine eta grid
tau_grid, eta_grid = np.meshgrid(tau_unique, eta_unique)



# Interpolate |theta| on the (tau, eta) grid
un_interp = griddata((grouped_data.index.get_level_values('tau'), 
                         grouped_data.index.get_level_values('eta')),
                        grouped_data['u.n'],
                        (tau_grid, eta_grid),
                        method='cubic')


# Step 4: Plot |theta| as a function of (tau, eta)
plt.figure(figsize=(10, 8))
plt.contourf(eta_grid, tau_grid, un_interp, levels=100, cmap='viridis')
plt.colorbar(label=r'$u\cdot n$')

# Use LaTeX in labels
plt.ylabel(r'$\tau$', fontsize=14)
plt.xlabel(r'$\eta$', fontsize=14)
# plt.title(r'Interpolated Mean $\langle |\theta| \rangle$ as a function of $\tau$ and $\eta$', fontsize=16)
plt.show()


# Interpolate |theta| on the (tau, eta) grid
sigma_interp = griddata((grouped_data.index.get_level_values('tau'), 
                         grouped_data.index.get_level_values('eta')),
                        grouped_data['|sigma|'],
                        (tau_grid, eta_grid),
                        method='cubic')

# Step 4: Plot |theta| as a function of (tau, eta)
plt.figure(figsize=(10, 8))
plt.contourf(eta_grid, tau_grid, sigma_interp, levels=100, cmap='viridis')
plt.colorbar(label=r'$\langle |\sigma| \rangle$')

# Use LaTeX in labels
plt.ylabel(r'$\tau$ (Proper Time)', fontsize=14)
plt.xlabel(r'$\eta$ (Spatial Rapidity)', fontsize=14)
plt.title(r'Interpolated Mean $\langle |\sigma| \rangle$ as a function of $\tau$ and $\eta$', fontsize=16)
plt.show()

data = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)

t_bins = np.linspace(data['tau'].min(), data['tau'].max(), 10)  # 10 bins for x
e_bins = np.linspace(data['eta'].min(), data['eta'].max(), 10)  # 10 bins for y

x_unique = np.linspace(data['x'].min(), data['x'].max(), 100)  # Fine tau grid
y_unique = np.linspace(data['y'].min(), data['y'].max(), 100)  # Fine eta grid
x_grid, y_grid = np.meshgrid(x_unique, y_unique)


un_interp = griddata((grouped_data.index.get_level_values('x'), 
                         grouped_data.index.get_level_values('y')),
                        grouped_data['u.n'],
                        (x_grid, y_grid),
                        method='cubic')


# Step 4: Plot |theta| as a function of (tau, eta)
plt.figure(figsize=(10, 8))
plt.contourf(x_grid, y_grid, un_interp, levels=100, cmap='viridis')
plt.colorbar(label=r'$u\cdot n$')


# Use LaTeX in labels
plt.ylabel(r'$x$', fontsize=14)
plt.xlabel(r'$y$', fontsize=14)
# plt.title(r'Interpolated Mean $\langle |\theta| \rangle$ as a function of $\tau$ and $\eta$', fontsize=16)
plt.show()

# Interpolate |theta| on the (tau, eta) grid
