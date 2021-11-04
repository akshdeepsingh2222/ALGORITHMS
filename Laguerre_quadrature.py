############# PART A AND PART C QUESTION-1#####
import matplotlib.pyplot as plt
from scipy import integrate
from sympy import *
import numpy as np
############# PART A AND PART C QUESTION-1#####
print("PART A AND PART C QUESTION-1")
def final_func(s):
    n=10    # this n is no. of quadrature point
    function_n=1   # this n is the n given in the function
    def func():
        x=symbols('x')
        return x**function_n
    
    def f():
        x=symbols('x')
        z2=func()
        integrand=z2*E**(-s*x)
        return integrand 
    
    def integration(n):
        x=symbols('x')
        y=0
        # n=10
        xi,wi=np.polynomial.laguerre.laggauss(n)
        for i in range(n):
            z=f()
            y+=wi[i]*np.exp(xi[i])*z.subs(x,xi[i])
        return y
    return integration(n)

s=np.arange(0,11,1)
print("S- VALUES :",s)
integrand=[]
for j in range(0,11,1):
    y1=final_func(s[j])
    integrand.append(y1)
print("I(s) VALUES FOR DIFFERENT S:",integrand)
plt.plot(s,integrand,'r*-',label="n=1")
plt.xlabel("S")
plt.ylabel("I(S)")
plt.title("PART A AND PART C QUESTION-1--LAPLACE'S TRANSFORM")
plt.legend()
plt.grid(True)
plt.show()
"""
output
PART A AND PART C QUESTION-1
S- VALUES : [ 0  1  2  3  4  5  6  7  8  9 10]
I(s) VALUES FOR DIFFERENT S: [614.215464721220, 1.00000000000000, 0.250000013112609, 0.111121913632303, 0.0626559302557678, 0.0405926814046912, 0.0290120226219872, 0.0222773034098904, 0.0179736120021536, 0.0149692867182496, 0.0127096477679539]
"""

