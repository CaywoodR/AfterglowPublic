#Fit a function to extraction, test many functions
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Function Selection 
def func(x, a, b, c):
    return b*((x**(-(3*a)))) 

# Read and preprocess
filename = "extract.txt"  # Update 
with open(filename, 'r') as file:
    lines = file.readlines()

# Clean lines
cleaned_lines = [re.sub(r'\[|\]', '', line.strip()) for line in lines]

# Convert 
cleaned_data = np.genfromtxt(cleaned_lines, skip_header=1)

# Extract 
x = cleaned_data[:, 0]
y1 = cleaned_data[:, 1]
y2_vectors = cleaned_data[:, 2:]


y2 = y2_vectors[:, 0]  # Change 0 to the index of the element you want to extract. Come back later after presentation

# Convert data to float
x = x.astype(float)
y1 = y1.astype(float)
y2 = y2.astype(float)

# Filter It's to late for plotting errors...
valid_indices = ~np.isnan(x) & ~np.isnan(y1) & ~np.isnan(y2)
x = x[valid_indices]
y1 = y1[valid_indices]
y2 = y2[valid_indices]

# Perform curve fit and pop coefficients at same time?
popt1, pcov1 = curve_fit(func, x, y1)
popt2, pcov2 = curve_fit(func, x, y2)

print("For y1:")
print("a =", popt1[0])
print("b =", popt1[1])
print("c =", popt1[2])

print("For y2:")
print("a =", popt2[0])
print("b =", popt2[1])
print("c =", popt2[2])

# Generate y values for the fitted curves and write functions to file
y1_fit = func(x, *popt1)
y2_fit = func(x, *popt2)
# Writing y1_fit and y2_fit in COLUMNS
with open("function_in.txt", 'w') as file:
    file.write("y1_fit\t\t\t y2_fit\n")  
    for val1, val2 in zip(y1_fit, y2_fit):
        file.write("{:<10}\t {:<10}\n".format(val1, val2))

#Plot
plt.figure(figsize=(8, 6))
plt.scatter(x, y1, label='E_tot')
plt.scatter(x, y2, label='Gamma')
plt.plot(x, y1_fit, 'r-', label='E_tot fit')
plt.plot(x, y2_fit, 'g-', label='Gamma fit')
plt.xlabel("h [$^{\circ}$]")
plt.ylabel('')
plt.xlim(0,90)
plt.yscale('log')
plt.title('(Gamma , E_tot) x h')
plt.legend()
plt.grid(False)
plt.show()
