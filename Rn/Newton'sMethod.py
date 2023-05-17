#!/usr/bin/env python
# coding: utf-8

# In[3]:


from math import e
import numpy as np
from numpy.linalg import norm


# In[4]:


Z = np.array([1,2])


# In[83]:


def f(Z):
    return (Z[0]**2 + 2*Z[1]**2 + e**(Z[0]**2 + Z[1]**2) - Z[0] + 2*Z[1])

def grad1(Z):
    return np.array([2*Z[0] * (1 + e**(Z[0]**2 + Z[1]**2)) - 1,
                     2*Z[1] * (2 + e**(Z[0]**2 + Z[1]**2)) + 2])

def grad2(Z):
    return np.array([[2*(1 + e**(Z[0]**2 + Z[1]**2)) + 2*Z[0] * 2*Z[0] * e**(Z[0]**2 + Z[1]**2),
                     2*Z[0] * 2*Z[1] * (1 + e**(Z[0]**2 + Z[1]**2))],
                     [2*Z[0] * 2*Z[1] * (1 + e**(Z[0]**2 + Z[1]**2)),
                     2*(2 + e**(Z[0]**2 + Z[1]**2)) + 2*Z[1] * 2*Z[1] * e**(Z[0]**2 + Z[1]**2)]])


# In[84]:


def arg(X):
    return np.linalg.solve(grad2(Z), np.dot(grad2(Z), Z) - X * grad1(Z))

lyambda = (1 + 5**(1/2))/2

def GoldenRatio(Z, lyambda, A=0, B=1, eps1=0.00000001, count=0):

    delta = (B-A)/lyambda**2
    X = A + delta
    Y = B - delta
    if (abs(B - A) <= 2 * eps1):
        S = (A+B)/2
    else:
        while (abs(B - A) > 2 * eps1):
            if (f(arg(X)) > f(arg(Y))):
                A = X;
                X = Y;
                Y = A + B - X;
            else:
                B = Y;
                Y = X;
                X = A + B - Y;
        S = (A+B)/2
    return S


# In[85]:


eps = 0.0000000001
count = 0
while (norm(grad1(Z)) >= eps):
    #Z = Z - np.dot(np.linalg.inv(grad2(Z)), grad1(Z))
    #S = GoldenRatio(Z, lyambda)
    Z = arg(1)
    count += 1
    print(f"Нормер итерации: {count} Z = [{Z[0]:.10f} {Z[1]:.10f}]  f(Z) = {f(Z):.10f}  |f'(Z)| = {norm(grad1(Z)):.10f}")


# In[ ]:




