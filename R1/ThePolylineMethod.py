#!/usr/bin/env python
# coding: utf-8

# In[196]:


import numpy as np
from math import cos, sin


# In[248]:


def f(x):
    return cos(x)/x**2

def df(x):
    return -sin(x)/x**2 - 2*cos(x)/x**3

def get_R(a, b, h):
    X = np.arange(a, b+h, h)
    R = abs(df(X[0]))
    for i in range(X.size-1):
        if abs(df(X[i])) >= R:
            R = df(X[i])
    return R + h

def x(a, b, R):
    return (a+b)/2 + (f(a)-f(b))/(2*R)

def m(a, b, R):
    return R*(a-b)/2 + (f(a)+f(b))/2


# In[249]:


def Method(a, b, h, eps, get_R, x, m, count=1, q=0, k=0):
    R = get_R(a, b, h)
    num = 0
    digits = eps
    while digits < 1:
        digits *= 10
        num += 1
    X = np.array([])
    M = np.array([])
    x = x(a, b, R)
    m = m(a, b, R)
    while f(x) - m >= eps:
        delta = (f(x)-m)/2/R
        print(f"Номер итерации: {count}  x = {x:.{num}f}  f(x) = {f(x):.{num}f}  m = {m:.{num}f}  delta = {delta:.{num}f}")
        X = np.append(X, [x-delta, x+delta])
        M = np.append(M, [(f(x)+m)/2, (f(x)+m)/2])
        q = M[0]
        for i in range(M.size-1):
            if M[i] <= q:
                q = M[i]
                k = i
        M = np.delete(M, k)
        x = X[k]
        X = np.delete(X, k)
        m = q
        count += 1
    print(f"Номер итерации: {count}  x = {x:.{num}f}  f(x) = {f(x):.{num}f}  R = {R:.{num}f}")


# In[252]:


Method(7, 11, 0.1, 0.0000001, get_R, x, m)


# In[ ]:




