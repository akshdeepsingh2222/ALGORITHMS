############ PART B QUES 2############
import matplotlib.pyplot as plt
from scipy import integrate
from sympy import *
import numpy as np
import math as mt
########### PART B QUES 2############
print("PART B QUES 2")
print("normalisation-question")
n=10
def f(sigma):
    a=2

    x=symbols('x')
    h=(x-a)**2/(2*sigma**2)
    integrand=(E**(-h))/np.sqrt(2*np.pi*sigma)
    h1=(x+a)**2/(2*sigma**2)
    integrand1=(E**(-h1))/np.sqrt(2*np.pi*sigma)
    return (integrand*integrand1)**2
print("NEW FUNCTION AFTER MULTIPLICATION",f(0.01))    # INTEGRAND 
def integration(n):
    x=symbols('x')
    y=0
    sigma=1
    for i in range (n):
        z=f(sigma)
        xi,wi=np.polynomial.hermite.hermgauss(n)
        y+=wi[i]*np.exp(xi[i]**2)*z.subs(x,xi[i])
    return y

print("VALUE AFTER NORMALISATION :",integration(10))

"""
output:
    PART B QUES 2
normalisation-question
NEW FUNCTION AFTER MULTIPLICATION 253.302959105844*exp(-10000.0*(x - 2)**2)*exp(-10000.0*(x + 2)**2)
VALUE AFTER NORMALISATION : 1.06496326962636e-5
"""