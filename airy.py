# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 21:44:56 2025

@author: arunb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import airy
import os
# -------------------------------
# Initial Parameters (Choose your initial values and method)
# -------------------------------
ss = 2
xx = 4
zz = 1
initial=-50
final=3
n = 1000 

#Choose Method -> "lie_trotter", "strang", "third", "sixth"
method="sixth"


# -------------------------------
# Do not change anything below this
# -------------------------------

total1=abs(initial)+ss
total2=final-ss
n1=int(np.floor((total1 / (total1+total2) * n)))
n2=n-n1
dt=-total1/(n1)
dt2=total2/(n2)


# -------------------------------
# Sigma and Gamma definitions for Airy
# -------------------------------

def sigmaLor(t, x, y, z):
    return t+x

def gammaLor(t, x, y, z):
    return t*z+y

def tauLor(t, x, y, z):
    return t*x*y +z

# Initialize arrays
x = np.zeros(n1+1)
y = np.zeros(n1+1)
z = np.zeros(n1+1)

x2 = np.zeros(n2)
y2 = np.zeros(n2)
z2 = np.zeros(n2)

# Initial conditions
x[0] = ss
y[0] = xx
z[0] = zz

x2[0] = ss
y2[0] = xx
z2[0] = zz

#complexcoefficient for sixth method
def compute_a_k(k):
    # Compute the exponent for the complex number
    return np.exp((1j * np.pi) / (2 * k + 1)) /(2**(1 / (2 * k + 1)) + 2 * np.exp((1j * np.pi )/ (2 * k + 1)))

aa = compute_a_k(1) * compute_a_k(2)
bb = compute_a_k(2) * (1 - 2 * compute_a_k(1))
cc = compute_a_k(1) * (1 - 2 * compute_a_k(2))
dd = (1 - 2 * compute_a_k(1)) * (1 - 2 * compute_a_k(2))

if method == "lie_trotter":
    for k in range(n1):
        x[k+1] = sigmaLor(dt, x[k], y[k], z[k] )
        y[k+1] = gammaLor(dt, x[k+1], y[k], z[k])
        z[k+1] = tauLor(dt, x[k+1], y[k+1], z[k])
    
    for k in range(n2-1):
        x2[k+1] = sigmaLor(dt2, x2[k], y2[k], z2[k] )
        y2[k+1] = gammaLor(dt2, x2[k+1], y2[k], z2[k])
        z2[k+1] = tauLor(dt2, x2[k+1], y2[k+1], z2[k])
        
elif method == "strang":
    for k in range(n1):
        z_1 = tauLor(dt/2, x[k], y[k], z[k])
        y_2 = gammaLor(dt/2, x[k], y[k], z_1)
        x[k+1] = sigmaLor(dt, x[k], y_2, z_1)
        y[k+1] = gammaLor(dt/2, x[k+1], y_2, z_1)
        z[k+1] = tauLor(dt/2, x[k+1], y[k+1], z_1) 
    
    
    for k in range(n2-1):
        z_1 = tauLor(dt2/2, x2[k], y2[k], z2[k])
        y_2 = gammaLor(dt2/2, x2[k], y2[k], z_1)
        x2[k+1] = sigmaLor(dt2, x2[k], y_2, z_1)
        y2[k+1] = gammaLor(dt2/2, x2[k+1], y_2, z_1)
        z2[k+1] = tauLor(dt2/2, x2[k+1], y2[k+1], z_1)
    
elif method == 'third':
    aT = 1/2 + 1j * np.sqrt(3)/6
    aTb=1/2 - 1j * np.sqrt(3)/6
        
    for k in range(n1):
        y_1 = gammaLor(aT*dt/2, x[k], y[k], z[k])
        z_2 = tauLor(aT*aT*dt/2, x[k], y_1, z[k])
        x_3 = sigmaLor(aT*aT*dt, x[k], y_1, z_2)
        z_4 = tauLor(aT*dt/2, x_3, y_1, z_2)
        x_5 = sigmaLor(aT*aTb*dt, x_3, y_1, z_4 )
        z_6 = tauLor(aT*aTb*dt/2, x_5, y_1, z_4)
        y_7 = gammaLor(dt/2, x_5, y_1, z_6,)
        z_8 = tauLor(aTb*aT*dt/2, x_5, y_7, z_6,)
        x_9 = sigmaLor(aT*aTb*dt, x_5, y_7, z_8)
        z_10 = tauLor(aTb*dt/2, x_9, y_7, z_8)
        x[k+1] = float(sigmaLor(aTb * aTb * dt, x_9, y_7, z_10).real)
        #x[k+1] = sigmaLor(aTb*aTb*dt, x_9, y_7, z_10)
        z[k+1] = float(tauLor(aTb*aTb*dt/2, x[k+1], y_7, z_10).real)
        y[k+1] = float(gammaLor(aTb*dt/2, x[k+1], y_7, z[k+1]).real) 
    
    
    for k in range(n2-1):
        y_1 = gammaLor(aT*dt2/2, x2[k], y2[k], z2[k])
        z_2 = tauLor(aT*aT*dt2/2, x2[k], y_1, z2[k])
        x_3 = sigmaLor(aT*aT*dt2, x2[k], y_1, z_2)
        z_4 = tauLor(aT*dt2/2, x_3, y_1, z_2)
        x_5 = sigmaLor(aT*aTb*dt2, x_3, y_1, z_4 )
        z_6 = tauLor(aT*aTb*dt2/2, x_5, y_1, z_4)
        y_7 = gammaLor(dt2/2, x_5, y_1, z_6,)
        z_8 = tauLor(aTb*aT*dt2/2, x_5, y_7, z_6,)
        x_9 = sigmaLor(aT*aTb*dt2, x_5, y_7, z_8)
        z_10 = tauLor(aTb*dt2/2, x_9, y_7, z_8)
        x2[k+1] = float(sigmaLor(aTb*aTb*dt2, x_9, y_7, z_10).real)
        z2[k+1] = float(tauLor(aTb*aTb*dt2/2, x2[k+1], y_7, z_10).real)
        y2[k+1] = float(gammaLor(aTb*dt2/2, x2[k+1], y_7, z2[k+1]).real)


elif method == 'sixth':
    def z1loop(inpt,a,b,c):
        z_1a = tauLor(inpt *aa * dt / 2, a, b, c)
        x_2a = sigmaLor(inpt *aa* dt, a, b, z_1a,)
        z_3a = tauLor(inpt  *(aa + bb) * dt / 2, x_2a,b, z_1a)
        x_4a = sigmaLor(inpt *bb * dt, x_2a, b, z_3a)
        z_5a = tauLor(inpt  *(aa + bb) * dt/2 , x_4a, b, z_3a)
        x_6a = sigmaLor(inpt  *aa * dt, x_4a,b, z_5a )
        z_7a = tauLor(inpt  *(aa + cc) * dt / 2, x_6a, b, z_5a)
        x_8a = sigmaLor(inpt  *cc * dt, x_6a, b, z_7a)
        z_9a = tauLor(inpt  *(cc + dd) * dt / 2, x_8a, b, z_7a)
        x_10a = sigmaLor(inpt  *dd * dt, x_8a,b, z_9a)
        z_11a = tauLor(inpt *(cc + dd) * dt / 2, x_10a,b, z_9a)
        x_12a = sigmaLor(inpt  *cc * dt, x_10a, b, z_11a)
        z_13a = tauLor(inpt  *(aa + cc) * dt / 2, x_12a,b, z_11a)
        x_14a = sigmaLor(inpt  *aa * dt, x_12a,b, z_13a )
        z_15a = tauLor(inpt  *(aa + bb) * dt/2, x_14a, b, z_13a)
        x_16a = sigmaLor(inpt  * bb * dt, x_14a, b, z_15a )
        z_17a = tauLor(inpt *(aa + bb) * dt / 2, x_16a, b, z_15a)
        xx = sigmaLor(inpt *aa * dt, x_16a, b, z_17a )
        zz = tauLor(inpt *aa * dt / 2, xx, b, z_17a)
        return xx , b, zz
    
    
    def z2loop(inpt,a,b,c):
        yy= gammaLor(inpt*dt/2, a, b, c)
        return a, yy, c
    
    def zaloop3( xp, yq, zq):
        x1, y1, z1=z2loop(aa, xp, yq, zq)
        x2, y2, z2=z1loop(aa, x1, y1, z1)
        x3, y3, z3=z2loop((aa + bb), x2, y2, z2)
        x4, y4, z4=z1loop(bb,x3, y3, z3)
        x5, y5, z5=z2loop((aa + bb), x4, y4, z4)
        x6, y6, z6=z1loop(aa,x5, y5, z5)
        x7, y7, z7=z2loop((aa +cc), x6, y6, z6)
        x8, y8, z8=z1loop(cc,x7, y7, z7)
        x9, y9, z9=z2loop((cc + dd), x8, y8, z8)
        x10, y10, z10=z1loop(dd,x9, y9, z9)
        x11, y11, z11=z2loop((cc + dd),x10, y10, z10)
        x12, y12, z12=z1loop(cc,x11, y11, z11)
        x13, y13, z13=z2loop((aa + cc),x12, y12, z12)
        x14, y14, z14=z1loop(aa, x13, y13, z13)
        x15, y15, z15=z2loop((aa + bb),x14, y14, z14)
        x16, y16, z16=z1loop(bb,x15, y15, z15)
        x17, y17, z17=z2loop((aa + bb),x16, y16, z16)
        x18, y18, z18=z1loop(aa, x17, y17, z17)
        x19, y19, z19=z2loop(aa,x18, y18, z18)
    
        return x19, y19, z19
    
    for k in range(n1):
        x[k+1], y[k+1],z[k+1]= zaloop3(x[k], y[k],z[k])
    
    
    ##left part
    def z1loopb(inpt,a,b,c):
        z_1a = tauLor(inpt *aa * dt2 / 2, a, b, c)
        x_2a = sigmaLor(inpt *aa* dt2, a, b, z_1a,)
        z_3a = tauLor(inpt  *(aa + bb) * dt2 / 2, x_2a,b, z_1a)
        x_4a = sigmaLor(inpt *bb * dt2, x_2a, b, z_3a)
        z_5a = tauLor(inpt  *(aa + bb) * dt2/2 , x_4a, b, z_3a)
        x_6a = sigmaLor(inpt  *aa * dt2, x_4a,b, z_5a )
        z_7a = tauLor(inpt  *(aa + cc) * dt2 / 2, x_6a, b, z_5a)
        x_8a = sigmaLor(inpt  *cc * dt2, x_6a, b, z_7a)
        z_9a = tauLor(inpt  *(cc + dd) * dt2 / 2, x_8a, b, z_7a)
        x_10a = sigmaLor(inpt  *dd * dt2, x_8a,b, z_9a)
        z_11a = tauLor(inpt *(cc + dd) * dt2 / 2, x_10a,b, z_9a)
        x_12a = sigmaLor(inpt  *cc * dt2, x_10a, b, z_11a)
        z_13a = tauLor(inpt  *(aa + cc) * dt2 / 2, x_12a,b, z_11a)
        x_14a = sigmaLor(inpt  *aa * dt2, x_12a,b, z_13a )
        z_15a = tauLor(inpt  *(aa + bb) * dt2/2, x_14a, b, z_13a)
        x_16a = sigmaLor(inpt  * bb * dt2, x_14a, b, z_15a )
        z_17a = tauLor(inpt *(aa + bb) * dt2 / 2, x_16a, b, z_15a)
        xx = sigmaLor(inpt *aa * dt2, x_16a, b, z_17a )
        zz = tauLor(inpt *aa * dt2 / 2, xx, b, z_17a)
        return xx , b, zz
    
    
    def z2loopb(inpt,a,b,c):
        yy= gammaLor(inpt*dt2/2, a, b, c)
        return a, yy, c
    
    def zaloop3b( xp, yq, zq):
        x1, y1, z1=z2loopb(aa, xp, yq, zq)
        x2, y2, z2=z1loopb(aa, x1, y1, z1)
        x3, y3, z3=z2loopb((aa + bb), x2, y2, z2)
        x4, y4, z4=z1loopb(bb,x3, y3, z3)
        x5, y5, z5=z2loopb((aa + bb), x4, y4, z4)
        x6, y6, z6=z1loopb(aa,x5, y5, z5)
        x7, y7, z7=z2loopb((aa +cc), x6, y6, z6)
        x8, y8, z8=z1loopb(cc,x7, y7, z7)
        x9, y9, z9=z2loopb((cc + dd), x8, y8, z8)
        x10, y10, z10=z1loopb(dd,x9, y9, z9)
        x11, y11, z11=z2loopb((cc + dd),x10, y10, z10)
        x12, y12, z12=z1loopb(cc,x11, y11, z11)
        x13, y13, z13=z2loopb((aa + cc),x12, y12, z12)
        x14, y14, z14=z1loopb(aa, x13, y13, z13)
        x15, y15, z15=z2loopb((aa + bb),x14, y14, z14)
        x16, y16, z16=z1loopb(bb,x15, y15, z15)
        x17, y17, z17=z2loopb((aa + bb),x16, y16, z16)
        x18, y18, z18=z1loopb(aa, x17, y17, z17)
        x19, y19, z19=z2loopb(aa,x18, y18, z18)
    
        return x19, y19, z19
    
    for k in range(n2-1):
        x2[k+1], y2[k+1],z2[k+1]= zaloop3b(x2[k], y2[k],z2[k])
       
    
else:
    raise ValueError("Invalid method. Choose one method.")

# -------------------------------
# Airy Method using Python Ai's and Bi's for comparison
# -------------------------------

t_vals=np.concatenate((x[::-1], x2))
step_size =(t_vals[1]) -(t_vals[0])
# Evaluate Airy functions for entire time range
Ai, Aip, Bi, Bip = airy(t_vals)

# Evaluate Ai, Ai', Bi, Bi' at ss to solve for constants
Ai_ss, Aip_ss, Bi_ss, Bip_ss = airy(ss)

# Solve the system: c1 * Ai(ss) + c2 * Bi(ss) = xx
#                   c1 * Ai'(ss) + c2 * Bi'(ss) = zz
A = np.array([[Ai_ss, Bi_ss], [Aip_ss, Bip_ss]])
b = np.array([xx, zz])
c1, c2 = np.linalg.solve(A, b)

# Compute full Airy solution
y_airy = c1 * Ai + c2 * Bi
dy_airy = c1 * Aip + c2 * Bip


# -------------------------------
# Plot 
# -------------------------------
os.makedirs("plots", exist_ok=True)
print("Current working directory is:")
print(os.getcwd())

plt.figure(figsize=(10, 5))
plt.plot(x2, y2, 'g-', label="")
plt.plot(x2, z2, 'm--', label="")
plt.plot(x, y, 'g-', label='y(t) Sixth')
plt.plot(x, z, 'm--', label="y'(t) Sixth")

plt.plot(t_vals, y_airy, 'r-', label=r'$A_i(t)$ Python')
plt.plot(t_vals, dy_airy, 'b--', label=r"$B_i(t)$ Python")

plt.axvline(ss, color='k', linestyle=':', label='s')
plt.title("")
plt.xlabel("t")
#plt.ylabel("Value")
plt.legend(fontsize='large') 
plt.grid(True)
plt.xlim(initial, final)
plt.savefig("plots/airy_plot.pdf", format='pdf', bbox_inches='tight')
plt.show()

# -------------------------------
# RMSE Calculation
# -------------------------------
yc=np.concatenate((y[::-1], y2))
zc=np.concatenate((z[::-1], z2))

rmse_y = np.sqrt(np.mean((yc- y_airy)**2 ))
rmse_dy= np.sqrt(np.mean((zc- dy_airy)**2 ))
print(f"\nRMSE of y between {method} and Python: {rmse_y:.4f}")
print(f"\nRMSE of dy between {method} and Python: {rmse_dy:.4f}")



# -------------------------------
# Pointwise Error Visualization
# -------------------------------

error = yc - y_airy
errorZ = zc - dy_airy
plt.plot(t_vals, error, label="Error y(t)")
plt.plot(t_vals, errorZ, label="Error y'(t)")
plt.title("Pointwise Error in y(t)")
plt.grid(True)
plt.legend()
plt.show()
