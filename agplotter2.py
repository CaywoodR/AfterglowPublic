import math
import matplotlib.pyplot as plt
tvals = []
tvals2 = []
data_sets = [[] for _ in range(9)]
sum_fx_data=[]

with open('roft2.txt', 'r') as f:
    for line in f:
        data = line.split()
        tvals.append(float(data[0]))
        for i in range(8):
            data_sets[i].append(float(data[i+1]))

with open('sum_fx2.txt', 'r') as f:
    for line in f:
        data = line.split()
        tvals2.append(float(data[0]))  
        sum_fx_data.append(float(data[1]))
            
labels = ['R(t)', r'$\beta$(t)', r'$\Gamma$(t)', 'Doft', 'Bp', 'nup', 'Ppnpm', 'F(t)']

plt.ion()
'''for i in range(8):
    plt.loglog(tvals, data_sets[i], label=f'Data Set {i+1}')

plt.xlabel('Time (t)')
plt.ylabel('Data Sets')
plt.title('Plot of tvals vs 6 Data Sets')
plt.ylim(bottom=0)
plt.legend()
plt.show()'''

'''for i in range(0,1):
    plt.scatter(tvals, data_sets[i])
    plt.xscale('log')
    plt.yscale('log')'''
    
plt.loglog(tvals2, sum_fx_data, label = "Sum_Fx")

for i in range(7,8):
    plt.loglog(tvals, data_sets[i], label = f'{labels[i]}')
    plt.draw()
    plt.pause(0.01)
plt.pause(1.0)
plt.ioff()
plt.xlim(1,1000)
plt.ylabel("Flux [mJy]")
plt.xlabel("Time [Days]")
plt.title('Flux vs. Time')
plt.legend()
plt.show()
#print(len(tvals))
#print(len(sum_fx_data))