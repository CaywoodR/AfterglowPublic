import numpy as np
import matplotlib.pyplot as plt

# Data
data = np.loadtxt('gamma.txt')
col1 = data[:, 0]
col2 = data[:, 1]
col3 = data[:,2]

# Plot
plt.plot(col1, col2/(5.0E+50), 'b-', label='L($\Theta$)')
plt.plot(col1,col3, label = 'G($\Theta$)')
plt.xlabel('$\Theta$')
plt.yscale('log')
plt.axhline(y=1, color='black', linestyle ='--')
plt.axvline(x= 5.45, color = 'red', linestyle = '--', label = 'R_dec')
plt.legend(loc='best')
plt.xlim(0,50)
plt.title('Input Assumption: Energy and Lorentz Factor')

plt.show()
