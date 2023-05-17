#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from numpy.linalg import norm
from math import e


# In[2]:


def f(Z):
    return (Z[0]**2 + Z[1]**2 - 6*Z[0] - 4*Z[1])


def gradient(Z):
    return np.array([2*Z[0] - 6,
                     2*Z[1] - 4])


def fi(Z, i):
    if i == 0:
        return Z[0] + Z[1] - 2
    elif i == 1:
        return -Z[0]
    elif i == 2:
        return -Z[1]
    

def H(Z, k=3, p=2):
    sigma=0
    for i in range(k):
        f_i = fi(Z, i)
        if f_i > 0:
            sigma += f_i**p
    return sigma


def dH(Z, k=3, p=2):
    sigma = np.array([0, 0])
    for i in range(k):
        if fi(Z, i) > 0:
            if i == 0:
                sigma[0] += p*(Z[0] + Z[1] - 2)**(p-1)
                sigma[1] += p*(Z[0] + Z[1] - 2)**(p-1)
            elif i == 1:
                sigma[0] += -p*(-Z[0])**(p-1)
            elif i == 2:
                sigma[1] += -p*(-Z[1])**(p-1)
    return sigma


def fn(Z, a):
    return f(Z) + a*H(Z)


def dfn(Z, a):
    return np.array([2*Z[0] - 6 + a*dH(Z)[0], 
                     2*Z[1] - 4 + a*dH(Z)[1]])


# проба
def fn(Z, a, k=3, p=2):
    sigma = 0
    for i in range(k):
        if fi(Z, i) > 0:
            if i == 0:
                sigma += (Z[0] + Z[1] - 2)**(p)
            elif i == 1:
                sigma += (-Z[0])**(p)
            elif i == 2:
                sigma += (-Z[1])**(p)
    return (Z[0]**2 + Z[1]**2 - 6*Z[0] - 4*Z[1] + a*sigma)


# проба
def dfn(Z, a, k=3, p=2):
    sigma = np.array([2*Z[0] - 6, 2*Z[1] - 4])
    for i in range(k):
        if fi(Z, i) > 0:
            if i == 0:
                sigma[0] += p*a*(Z[0] + Z[1] - 2)**(p-1)
                sigma[1] += p*a*(Z[0] + Z[1] - 2)**(p-1)
            elif i == 1:
                sigma[0] += -p*a*(-Z[0])**(p-1)
            elif i == 2:
                sigma[1] += -p*a*(-Z[1])**(p-1)
    return np.array([sigma[0], sigma[1]])
        
    
lyambda = (1 + 5**(1/2))/2
def GoldenRatio(Z, a, lyambda, A=0, H=0.00000001, eps1=0.000000001, count=0):
    if (fn((Z - A*dfn(Z, a)/norm(dfn(Z, a))), a) < fn((Z - H*dfn(Z, a)/norm(dfn(Z, a))), a)):
        B = H
    else:
        while (fn((Z - A*dfn(Z, a)/norm(dfn(Z, a))), a) >= fn((Z - H*dfn(Z, a)/norm(dfn(Z,a))), a)):
            H *= 2
            count += 1
        B = H
    delta = (B-A)/lyambda**2
    X = A + delta
    Y = B - delta
    if (abs(B - A) <= 2 * eps1):
        S = (A+B)/2
    else:
        while (abs(B - A) > 2 * eps1):
            if (fn((Z - X*dfn(Z, a)/norm(dfn(Z, a))), a) > fn((Z - Y*dfn(Z, a)/norm(dfn(Z, a))), a)):
                A = X;
                X = Y;
                Y = A + B - X;
            else:
                B = Y;
                Y = X;
                X = A + B - Y;
        S = (A+B)/2
    return S


# In[3]:


def SDescent(Z, a):
    eps2 = 0.00001
    count = 0
    if (norm(dfn(Z, a)) <= eps2):
        print(f"Нормер итерации: {count}    Z = [{Z[0]:.5f} {Z[1]:.5f}]   |fn'(Z)| = {norm(dfn(Z, a)):.5f}   fn(Z) = {fn(Z, a):.5f}   f(Z) = {f(Z):.5f}")
    else:
        while (norm(dfn(Z, a)) > eps2):
            S = GoldenRatio(Z, a, lyambda)
            Z = Z - S*dfn(Z, a)/norm(dfn(Z, a))
            count+=1
            print(f"Нормер итерации: {count}    Z = [{Z[0]:.5f} {Z[1]:.5f}]   |fn'(Z)| = {norm(dfn(Z, a)):.5f}   fn(Z) = {fn(Z, a):.5f}   f(Z) = {f(Z):.5f}")
    return Z


# In[4]:


Z = np.array([1, 2])
a = e
count = 1
eps3 = 0.001
if H(Z)**(1/2) >= eps3:
    print(f"\nШаг по штрафу {count}:   a = {a:.3f}     x = [{Z[0]:.3f} {Z[1]:.3f}]     H(x)^1/2 = {H(Z)**(1/2):.3f}\n")
Z = SDescent(Z, a)
while H(Z)**(1/2) >= eps3:
    count += 1
    a *= e
    print(f"\nШаг по штрафу {count}:   a = {a:.3f}     x = [{Z[0]:.3f} {Z[1]:.3f}]     H(x)^1/2 = {H(Z)**(1/2):.3f}\n")
    Z = SDescent(Z, a)
print(f"\nРезультат:   x = [{Z[0]:.3f} {Z[1]:.3f}]     H(x)^1/2 = {H(Z)**(1/2):.3f}\n")


# In[ ]:




