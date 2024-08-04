import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt


#User defined functions
####################
def Rap(eta, pt, m):
    return np.log((np.sqrt(m**2 + pt**2*np.cosh(eta)**2) + pt*np.sinh(eta)) / np.sqrt(m**2 + pt**2))
#####################


# Open the file and read the first line
with open('particleInformation.dat', 'r') as file:
    first_line = file.readline().strip()

# Split the line based on whitespace
values = first_line.split()

# Extract the values
m = float(values[0])
iphimax = int(values[6])
iptmax = int(values[5])
ietamax = int(values[2])
ptmin = float(values[3])
ptmax = float(values[4])
etamax = float(values[7])

################
npt = iptmax
nphi = iphimax
neta = ietamax
################
minrap = -0.5
maxrap = 0.5
###############

# Print the extracted values
print("iphimax:", iphimax)
print("iptmax:", iptmax)
print("ietamax:", ietamax)

# Read the data from the file
data = np.loadtxt('particleInformation.dat', skiprows=1)

nrows = len(data)

# Reshape the data to match the dimensions of dNdydptdphi
dNdydptdphi = data.reshape(ietamax,iptmax,iphimax)


# Calculate pt values
ipt_values = np.arange(iptmax)
pt_array = ptmin + (ptmax - ptmin) * (ipt_values / (iptmax - 1))**2

print(pt_array)


deltaeta = 2.0 * etamax / float(ietamax - 1)

# #Pseudo-rapidity array
y = np.arange(-etamax, etamax+deltaeta, deltaeta)

# Initialize rapidity array
ylist = np.empty(neta) #Rapidity

#print(y)
#print(ylist)


# Define the number of harmonics
nharmonics = 3

# Define the arrays and constants
intvn = np.zeros((iptmax, nharmonics, 2))
vn = np.zeros((nharmonics, 2, iptmax))

# Loop over pt
for ipt in range(npt):
    pt = pt_array[ipt]

    # Loop over phi
    for iphi in range(nphi):
        ylist = np.zeros(neta)
        etalist = np.zeros(neta)
        dndpt = np.zeros(neta)

        # Integrate over pseudorapidity
        for ieta in range(neta):
            rap_local = y[ieta]
            eta_local = rap_local
            etalist[ieta] = eta_local #Pseudorapidity
            y_local = Rap(eta_local, pt, m)
            ylist[ieta] = y_local #Rapidity

            #Check if you want dN/(dypTdpTdphi) or dN/(dydpTdphi)
            #and uncomment accordingly
            # Calculate dN/(dypTdpTdphi)
            dndpt[ieta] = dNdydptdphi[ieta][ipt][iphi]
            

            # Calculate dN/(dydpTdphi)
            #dndpt[ieta] = pt*dNdydptdphi[ieta][ipt][iphi]

        # Perform interpolation
        spline = interp1d(ylist, dndpt, kind='cubic')

        # Evaluate integral
        if minrap < maxrap:
            dNdp, _ = quad(spline, minrap, maxrap)
        else:
            dNdp = spline(maxrap)

        #print("dNdp=",dNdp)
        # Calculate phi and update intvn
        phi = iphi * 2 * np.pi / nphi
        for i in range(nharmonics):
            intvn[ipt][i][0] += np.cos(i * phi) * dNdp * 2 * np.pi / nphi
            intvn[ipt][i][1] += np.sin(i * phi) * dNdp * 2 * np.pi / nphi


    vn[0][0][ipt] = intvn[ipt][0][0]
    vn[0][1][ipt] = intvn[ipt][0][1]

    for i in range(1, nharmonics):
        vn[i][0][ipt] = intvn[ipt][i][0] / intvn[ipt][0][0]
        vn[i][1][ipt] = intvn[ipt][i][1] / intvn[ipt][0][0]

#################################

# Define the file name
file_name = 'vn_data.txt'

# Define the column headers
headers = ['pT', 'dN/pTdpT', 'v1', 'v2']  

# Write data to the file
with open(file_name, 'w') as file:
    # Write comment line with column headers
    file.write('# ' + ' '.join(headers) + '\n')

    # Write data rows
    for ipt, pt_value in enumerate(pt_array):
        row_data = ['{:.3e}'.format(pt_value)] + ['{:.3e}'.format(vn[i][0][ipt]) for i in range(nharmonics)]
        file.write(' '.join(map(str, row_data)) + '\n')

print(f'Data has been written to {file_name}.')

# Plot intvn[ipt][i][0] as a function of pT
plt.plot(pt_array, intvn[:,2, 0]/intvn[:, 0, 0], label=f'Harmonic {i+1}')

plt.xlabel('pT')
plt.ylabel('intvn[ipt][i][0]')
plt.title('intvn[ipt][i][0] vs. pT')
plt.legend()
plt.grid(True)
plt.show()

# Plot intvn[ipt][i][0] as a function of pT
plt.plot(pt_array, intvn[:, 0, 0], label=f'Harmonic {i+1}')

plt.xlabel('pT')
plt.ylabel('intvn[ipt][i][0]')
plt.title('intvn[ipt][i][0] vs. pT')
plt.legend()
plt.grid(True)
plt.show()
