#!/usr/bin/env python
# coding: utf-8

# In[58]:


from math import e
import numpy as np
from numpy.linalg import norm


# In[59]:


def f(Z):
    return (Z[0]**2 + 2*Z[1]**2 + e**(Z[0]**2 + Z[1]**2) - Z[0] + 2*Z[1])


# In[60]:


def gradient(Z):
    return np.array([2*Z[0] * (1 + e**(Z[0]**2 + Z[1]**2)) - 1,
                     2*Z[1] * (2 + e**(Z[0]**2 + Z[1]**2)) + 2])


# In[61]:


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


# In[127]:


# без обновления
Z = (0, 0)
eps2 = 0.000001
count = 0
h = gradient(Z)
n = norm(h)
gk_1 = h
if (n <= eps2):
    print(f"Нормер итерации: {count}  Z = [{Z[0]:.6f} {Z[1]:.6f}]  f(Z) = {f(Z):.6f}  |f'(Z)| = {n:.6f}")
else:
    while (n > eps2):
        S = GoldenRatio(Z,lyambda)
        Z = Z - S*h/norm(h)
        gk = gradient(Z)
        n = norm(gk)
        beta = np.dot(gk_1 - gk, gk) / norm(gk_1)**2
        h = gk - beta*h
        gk_1 = gk
        count+=1
        print(f"Нормер итерации: {count} Z = [{Z[0]:.6f} {Z[1]:.6f}]  f(Z) = {f(Z):.6f}  |f'(Z)| = {n:.6f}")


# In[128]:


#  с обновлением
Z = (0, 0)
eps2 = 0.000001
count = 0
k = 2
h = gradient(Z)
n = norm(h)
gk_1 = h
while (n > eps2):
    S = GoldenRatio(Z,lyambda)
    Z = Z - S*h/norm(h)
    count += 1
    gk = gradient(Z)
    n = norm(gk)
    if count%k != 0:
        beta = np.dot(gk_1 - gk, gk) / norm(gk_1)**2
    else:
        beta = 0
    h = gk - beta*h
    gk_1 = gk
    print(f"Нормер итерации: {count} Z = [{Z[0]:.6f} {Z[1]:.6f}]  f(Z) = {f(Z):.6f}  |f'(Z)| = {n:.6f}")


# In[ ]:




