#!/usr/bin/env python
# coding: utf-8

# In[1]:


from math import e
import numpy as np
from numpy.linalg import norm


# In[2]:


Z = np.array([2,4])


# In[3]:


def f(Z):
    return (Z[0]**2 + 2*Z[1]**2 + e**(Z[0]**2 + Z[1]**2) - Z[0] + 2*Z[1])


# In[4]:


def gradient(Z):
    return np.array([2*Z[0] * (1 + e**(Z[0]**2 + Z[1]**2)) - 1,
                     2*Z[1] * (2 + e**(Z[0]**2 + Z[1]**2)) + 2])


# In[4]:


lyambda = (1 + 5**(1/2))/2


# In[5]:


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


# In[6]:


eps2 = 0.000001
count = 0
print(f"Дано:   Z = {Z}  f(Z) = {f(Z):.6f}  |f'(Z)| = {norm(gradient(Z)):.6f}")
if (norm(gradient(Z)) <= eps2):
    print(f"Нормер итерации: {count}  Z = [{Z[0]:.6f} {Z[1]:.6f}]  f(Z) = {f(Z):.6f}  |f'(Z)| = {norm(gradient(Z)):.6f}")
else:
    while (norm(gradient(Z)) > eps2):
        S = GoldenRatio(Z,lyambda)
        Z = Z - S*gradient(Z)/norm(gradient(Z))
        count+=1
        print(f"Нормер итерации: {count} Z = [{Z[0]:.6f} {Z[1]:.6f}]  f(Z) = {f(Z):.6f}  |f'(Z)| = {norm(gradient(Z)):.6f}")


# In[ ]:




