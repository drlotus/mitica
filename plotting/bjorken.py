import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
filename = './output/bjorken_midrap.dat'
data = np.loadtxt(filename, skiprows=1)

# Extract the columns
mT = data[:, 0]
pT = data[:, 1]
dNd3p = data[:, 5]

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(mT, dNd3p, 'o-', label='$dNd3p$')

# Set y-axis to logarithmic scale
# plt.yscale('log')

# Add labels and title
plt.xlabel('$m_T$ (GeV)')
plt.ylabel('$dNd3p$ (GeV$^{-3}$)')
plt.title('$dNd3p$ vs $m_T$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()


 
plt.figure(figsize=(10, 6))
plt.plot(mT, dNd3p, 'o-', label='$dNd3p$')

# Set y-axis to logarithmic scale
plt.yscale('log')

# Add labels and title
plt.xlabel('m_T$ (GeV)')
plt.ylabel('$log(dNd3p)$ (GeV$^{-3}$)')
plt.title('$dNd3p$ vs $p_T$')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()
