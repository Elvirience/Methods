#!/usr/bin/env python
# coding: utf-8

# In[25]:


from numpy import sign
from math import cos


# In[26]:


def f(x):
    return cos(x)/x**2

def find_convex_triple(x0, a, b, eps, i=1):
    h = (b-a)/10
    while (f(x0+h)>f(x0))&(f(x0-h)>f(x0)):
        h /= 2
        if (f(x0+h)-f(x0)<=eps)|(f(x0-h)-f(x0)<=eps):
            return x0
    if f(x0+h)<=f(x0):
        x1 = x0 + h
        x2 = x0 + h*2**(i)
        f0 = f(x0)
        while (f0<f(x1))|(f(x2)<f(x1))|(f0+f(x2)-2*f(x1)<=0):
            i+=1
            x1 = x0 + h*2**(i-1)
            x2 = x0 + h*2**(i)
    else:
        x1 = x0 - h
        x2 = x0 - h*2**(i)
        f0 = f0
        while (f0<f(x1))|(f(x2)<f(x1))|(f0+f(x2)-2*f(x1)<=0):
            i+=1
            x1 = x0 - h*2**(i-1)
            x2 = x0 - h*2**(i)       
    return (x0, x1, x2)
            
def get_top(x1, x2, x3, f1, f2, f3):
    return (x1+x2)/2 + (f2-f1)*(x3-x2)*(x3-x1)/(2*(f1*(x2-x3)+f2*(x3-x1)+f3*(x1-x2)))

def sort(points):
    for i in range(3,0,-1):
        for j in range(i):
            if points[j][1] > points[j+1][1]:
                points[j], points[j+1] = points[j+1], points[j]
    return points


# In[27]:


def method(x0, eps, a, b, count=1):
    try:
        x1, x2, x3 = find_convex_triple(x0, a, b, eps)
    except:
        x0 = find_convex_triple(x0, a, b, eps)
    x4 = get_top(x1, x2, x3, f(x1), f(x2), f(x3))
    f1 = f(x1)
    f2 = f(x2)
    f3 = f(x3)
    f4 = f(x4)
    points = sort([(x1, f1), (x2, f2), (x3, f3), (x4, f4)])
    print(f"Номер итерации: {count}   x1 = {points[0][0]:.4f} f1 = {points[0][1]:.4f}\n\t\t    x2 = {points[1][0]:.4f} f2 = {points[1][1]:.4f}\n\t\t    x3 = {points[2][0]:.4f} f3 = {points[2][1]:.4f}\n\t\t    x4 = {points[3][0]:.4f} f4 = {points[3][1]:.4f}\n")
    while points[1][1] - points[0][1] >= eps:
        if sign(points[1][0] - points[0][0]) == sign(points[2][0] - points[0][0]) == -sign(points[3][0] - points[0][0]):
            count += 1
            x3 = points[3][0]
            x2 = points[1][0]
            x1 = points[0][0]
            x4 = get_top(x1, x2, x3, f(x1), f(x2), f(x3))
            f1 = f(x1)
            f2 = f(x2)
            f3 = f(x3)
            f4 = f(x4)
            points = sort([(x1, f1), (x2, f2), (x3, f3), (x4, f4)])
            print(f"Номер итерации: {count}   x1 = {x1:.4f} f1 = {f1:.4f}\n\t\t    x2 = {x2:.4f} f2 = {f2:.4f}\n\t\t    x3 = {x3:.4f} f3 = {f3:.4f}\n\t\t    x4 = {x4:.4f} f4 = {f4:.4f}\n")
        else:
            count += 1
            x3 = points[2][0]
            x2 = points[1][0]
            x1 = points[0][0]
            x4 = get_top(x1, x2, x3, f(x1), f(x2), f(x3))
            f1 = f(x1)
            f2 = f(x2)
            f3 = f(x3)
            f4 = f(x4)
            points = sort([(x1, f1), (x2, f2), (x3, f3), (x4, f4)])
            print(f"Номер итерации: {count}   x1 = {x1:.4f} f1 = {f1:.4f}\n\t\t    x2 = {x2:.4f} f2 = {f2:.4f}\n\t\t    x3 = {x3:.4f} f3 = {f3:.4f}\n\t\t    x4 = {x4:.4f} f4 = {f4:.4f}\n")
    print(f"argmin(f) = {points[0][0]:.4f}")


# In[28]:


method(9, 0.0001, 7, 11)


# In[ ]:




