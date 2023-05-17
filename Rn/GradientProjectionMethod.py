#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from numpy.linalg import norm


# In[2]:


Z = np.array([1, -4])
C = (1,1)
gamma = 2
k1=-1
k2 = 1
b1 = -2
b2 = 2

def direct0(x1, k1, gamma):
    return k1*x1+gamma

def direct1(x, k2, b1):
    return k2*x+b1

def direct2(x, k2, b2):
    return k2*x+b2

def Proj(Z, C, gamma, k2, b1, b2):
    if (Z[0]<0)&(Z[1]<0):
        print(0)
        return np.array([0,0])
    if (Z[1] <= direct1(Z[0], k2, b1))&(Z[0] >= 2):
        print(1)
        return np.array([2,0])
    if (Z[1] >= direct2(Z[0], k2, b2))&(Z[1]>=2):
        print(2)
        return np.array([0,2])
    if (Z[1]>direct1(Z[0], k2, b1))&(Z[1]<direct2(Z[0], k2, b2)):
        if Z[1] > direct0(Z[0],k1,gamma):
            print(3)
            return np.array(Z + np.dot((gamma - np.dot(C, Z)) / norm(C)**2, C))
        if (Z[1] <= direct0(Z[0],k1,gamma))&(Z[1] >= 0)&(Z[0] >= 0):
            print(3)
            return Z
    if (Z[1]<0)&(Z[0]>=0)&(Z[0]<2):
        return np.array([Z[0],0])
    if (Z[1]<2)&(Z[1]>=0)&(Z[0]<0):
        return np.array([0,Z[1]])
    
def f(Z):
    return (Z[0]**2 + Z[1]**2 - 6*Z[0] - 4*Z[1]) 

def gradient(Z):
    return np.array([2*Z[0] - 6,
                     2*Z[1] - 4])


# In[3]:


Proj(Z, C, gamma, k2, b1, b2)


# In[4]:


lyambda = (1 + 5**(1/2))/2

def GoldenRatio(Z, lyambda, A=0, H=1, eps1=0.0000001, count=0):
    if (f(Z - A*gradient(Z)/norm(gradient(Z))) < f(Z - H*gradient(Z)/norm(gradient(Z)))):
        B = H
    else:
        while (f(Z - A*gradient(Z)/norm(gradient(Z))) >= f(Z - H*gradient(Z)/norm(gradient(Z)))):
            H *=2
            count+=1
        B = H
    delta = (B-A)/lyambda**2
    X = A + delta
    Y = B - delta
    if (abs(B - A) <= 2 * eps1):
        S = (A+B)/2
    else:
        while (abs(B - A) > 2 * eps1):
            if (f(Z - X*gradient(Z)/norm(gradient(Z))) > f(Z - Y*gradient(Z)/norm(gradient(Z)))):
                A = X;
                X = Y;
                Y = A + B - X;
            else:
                B = Y;
                Y = X;
                X = A + B - Y;
        S = (A+B)/2
    return S


# In[5]:


eps2 = 0.000001
count = 0
S = GoldenRatio(Z,lyambda)
Pz = Proj(Z-S*gradient(Z)/norm(gradient(Z)), C, gamma, k2, b1, b2)
print(f"Дано:  Pz = {Pz}  S = {S:.4f}  |f'(Z)| = {norm(gradient(Z)):.4f}  |Z-Pz| = {norm(Z - Pz)}")
while (norm(gradient(Z)) >= eps2) & (norm(Z - Pz) >= eps2):
    Z = Pz
    S = GoldenRatio(Z,lyambda)
    Pz = Proj(Z-S*gradient(Z)/norm(gradient(Z)), C, gamma, k2, b1, b2)
    print(Pz)
    count += 1
    print(f"Нормер итерации: {count}  Pz = [{Pz[0]:.4f} {Pz[1]:.4f}]  S = {S:.4f}  f(Z) = {f(Z):.4f}  |f'(Z)| = {norm(gradient(Z)):.4f}  |Z-Pz| = {norm(Z - Pz):.4f}")


# In[ ]:




