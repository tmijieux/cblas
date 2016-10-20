#!/usr/bin/python2

import sys
import numpy as np
import matplotlib.pyplot as plt


table = {}
lines = sys.stdin.read().splitlines()
########################################


X = []
Y = []
Z = []
name = None
for line in lines:
    try:
        float(line.split()[0])
        int(line.split()[0])
    except:
        if name is not None:
            table.update({name: (X, Y, Z)})
        Z = []
        Y = []
        X = []
        name = line
        continue
    
    x, y, z = line.strip().split()
    X.append(x)
    Y.append(y)
    Z.append(z)
table.update({name: (X, Y, Z)})
    

print table.keys()

n,MF,T=table["dgemm_block"]
plt.plot(n,MF, color='y', label="dgemm_block")
n,MF,T=table["dgemm_j"]
plt.plot(n,MF, color='r',label="dgemm_j")
plt.legend()
plt.xlabel('Matrix size')
plt.ylabel('Mflops/s')

plt.show()

########################################

n,MF,T=table["dgemm_block"]
plt.plot(n,T, color='y',label="dgemm_block")
n,MF,T=table["dgemm_j"]
plt.plot(n,T, color='r',label="dgemm_j")
plt.legend()
plt.xlabel('Matrix size')
plt.ylabel('Time (Microseconds)')
plt.show()


########################################

# ax1.set_title("Performance Comparison")    
# ax1.set_xlabel('Matrix size')
# ax1.set_ylabel('Time (microseconds)')


# leg = ax1.legend()
# ax2.set_xlabel('Matrix size')
# ax2.set_ylabel('MFlops/s')
