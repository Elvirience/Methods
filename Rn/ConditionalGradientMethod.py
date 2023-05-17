#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
from numpy.linalg import norm

X = np.array([1, 1])
D = np.array([(0,0),(2,0),(0,2)])

def f(X):
    return (X[0]**2 + X[1]**2 - 6*X[0] - 4*X[1]) 

def grad(X):
    return np.array([2*X[0] - 6,
                     2*X[1] - 4])

def fk(x, xk):
    return np.dot(grad(xk), (x-xk))

def y(X,a=0):
    for i in range(2):
        if fk(D[i], X) > fk(D[i+1], X):
            a = D[i+1]
    return a

y(X)


# In[15]:


fk(y(X), X)


# In[16]:


lyambda = (1 + 5**(1/2))/2

def GoldenRatio(Z, lyambda, A=0, B=1, eps1=0.000001, count=0):
    delta = (B-A)/lyambda**2
    X = A + delta
    Y = B - delta
    if (abs(B - A) <= 2 * eps1):
        S = (A+B)/2
    else:
        while (abs(B - A) > 2 * eps1):
            if (f(Z + X*(y(Z)-Z))) > f(Z + Y*(y(Z)-Z)):
                A = X;
                X = Y;
                Y = A + B - X;
            else:
                B = Y;
                Y = X;
                X = A + B - Y;
        S = (A+B)/2
    return S


# In[17]:


X


# In[18]:


eps = 0.00001
count=1
S = GoldenRatio(X,lyambda)
XX = X + S*(y(X)-X)
print(f"Нормер итерации: {count}  X = [{X[0]:.5f} {X[1]:.5f}]  y = [{y(X)[0]:.5f} {y(X)[1]:.5f}]  S = {S:.5f}  f(X) = {f(X):.5f}   fk = {fk(y(X), X)}")
while(fk(y(X), X)<-eps):
    S = GoldenRatio(X,lyambda)
    X = X + S*(y(X)-X)
    count+=1
    print(f"Нормер итерации: {count}  X = [{X[0]:.5f} {X[1]:.5f}]  y = [{y(X)[0]:.5f} {y(X)[1]:.5f}]  S = {S:.5f}  f(X) = {f(X):.5f}   fk = {fk(y(X), X)}")


# In[ ]:




