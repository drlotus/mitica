import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
filename = './output/yield-proton-rn.dat'
data = np.loadtxt(filename, skiprows=1)

# Extract the columns
mT = data[:, 0]
pT = data[:, 1]
phi_p = data[:, 2]
y_p = data[:, 3]

dphi = 2*np.pi / phi_p.size

min_y = min(y_p, key=abs)
mass = min(mT)

# Filter the data for y_p = 0 and phi_p = pi/2
filtered_data = data[(y_p == min_y)]

filtered_data2 = data[(phi_p == 0)]


unique_mT = np.unique(filtered_data[:, 0])

unique_y_p = np.unique(filtered_data2[:, 3])


dN_over_mTdmTdy = []
for mt in unique_mT:
    mask = filtered_data[:, 0] == mt
    sum_dNd3p = np.sum(filtered_data[mask, 5]) 
    dN_over_mTdmTdy.append(sum_dNd3p)

# Convert to numpy array for plotting
dN_over_mTdmTdy = np.array(dN_over_mTdmTdy)


dN_over_dy = []
for y in unique_y_p:
    mask = filtered_data2[:, 3] == y
    sum_dNd3p = np.sum(filtered_data2[mask, 5]*filtered_data2[mask, 0]) 
    dN_over_dy.append(sum_dNd3p)

# Convert to numpy array for plotting
dN_over_mTdmTdy = np.array(dN_over_mTdmTdy)

dN_over_dy = np.array(dN_over_dy)

mt_min_m = [x - mass for x in unique_mT]



# Plotting the filtered data
plt.figure(figsize=(10, 6))
plt.plot(mt_min_m, dN_over_mTdmTdy, 'o-', label='$dN/(m_T  dm_T  dy)$')


# Add labels and title
plt.xlabel('$m_T-m$ (GeV)')
plt.ylabel('$dN /m_T d m_T d y$')
plt.title('(Summed over $\\phi_p$, $y_p = 0$)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()

# Plotting the filtered data
plt.figure(figsize=(10, 6))

plt.plot(mt_min_m, dN_over_mTdmTdy, 'o-', label='$dN/(m_T  dm_T  dy)$')

# Set y-axis to logarithmic scale
plt.yscale('log')

# Add labels and title


# Add labels and title
plt.xlabel('$m_T-m$ (GeV)')
plt.ylabel('$dN /m_T d m_T d y$')
plt.title('(Summed over $\\phi_p$, $y_p = 0$)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(unique_y_p, dN_over_dy, 'o-', label='$dN/dy$')


# Add labels and title
plt.xlabel('$y_p$')
plt.ylabel('$dN /d y$')
plt.title('')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()
