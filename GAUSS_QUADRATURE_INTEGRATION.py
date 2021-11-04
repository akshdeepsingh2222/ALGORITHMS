import numpy as np
from scipy.special import legendre
import scipy as sp
import math as mt
from sympy import *
import matplotlib.pyplot as plt
import pandas as pd
n=10
choice=2
a=0
b=5
xi,wi=np.polynomial.legendre.leggauss(n)
def inbuiltlegendre(x,n) :
    leg=sp.special.legendre(n)
    p_x=leg(x)
    return p_x
def nodal_points(n):
    x=symbols('x')
    xi=solve(legendre(n,x))
    return xi
# xi=nodal_points(n)
def f4(choice):
    x=symbols('x')
    if choice==1:
        n=3
        m=2
        return legendre(n,x)*legendre(m,x)
    if choice==2:
        return x**1*E**(-2*x)   #s=2
# def f2(choice):
#     x = Symbol('x')
#     if choice==1:
#       z=f4(choice)
#       y=integrate(z, x)
#       absolute=y.subs(x,b)-y.subs(x,a)
#       print("absolute value :",absolute)
#       return absolute 
#     if choice==2:
#          z=f4(choice)
#          y=integrate(z, x)
#          absolute=y.subs(x,b)-y.subs(x,a)
#          print("absolute value :",absolute)
#          return absolute
def weights(n):
    xi,wi=np.polynomial.legendre.leggauss(n)
    wn=[]
    for i in range (n):
        wi=2*(1-xi[i]**2)/(n*inbuiltlegendre(xi[i],n-1))**2
        wn.append(wi)
    return wn
def transform(a,b,choice):
    x=symbols('x')
    t=symbols('t')
    func=f4(choice)
    gt=func.subs(x,((b-a)/2)*t+(b+a)/2)
    # dt=diff(((b-a)/2)*t+(b+a)/2)
    return gt*((b-a)/2)
  
def gauss_integration(n,choice):
    t=symbols('t')
    xi,wi=np.polynomial.legendre.leggauss(n)
    wn=weights(n)
    func=0
    for i in range(n):
        func+=wn[i]*transform(a,b,choice).subs(t,xi[i])
    return func
print("ORIGINAL FUNCTION  is : ",f4(choice) ,"with limits ("+str(a)+","+str(b)+")")
print("TRANSFORMED FUNCTION is :",transform(a,b,choice),"with limits (-1,1)")
print("gauss legendre integration with "+str(n)+" points system :",gauss_integration(n,choice))
# y=np.arange(2,8,1)  # array for different n-pts
# integration=[]
# for i in range(len(y)):
#     integrand=gauss_integration(y[i],choice)
#     integration.append(integrand)
# list_of_tuples = list(zip(y,integration))    
# df = pd.DataFrame(list_of_tuples, columns = ['N', 'gauss integration'])  
# print(df)
# e=[]
# wn=weights(n)
# x_inbuilt,wi=np.polynomial.legendre.leggauss(n)
# for i in range(len(wn)):
#     err=abs(wn-wi)
#     e.append(err)  # e is error in weights
# k=[]
# absolute=f2(choice)
# for j in range(len(y)):
#     k1=(abs(integration[j]-absolute))
#     k.append(k1)
# k=np.asarray(k)
# x=np.linspace(a,b,100)

# plt.plot(y,k,"ro-")
# plt.yscale("log")
# plt.xscale("log")
# plt.ylabel("log error")
# plt.xlabel("n-pts")
# plt.grid("log")
# plt.show()
# ax1=plt.subplot(1,1,1)
# ax1.bar(y,k,color='#609000')
# ax2=ax1.twinx()
# plt.plot(y,k,'ro-')
# # plt.yscale("log")   ##
# plt.ylabel("log error")
# plt.xlabel("n-pts")
# plt.title("log error V/S n-points")
# plt.grid()
# plt.show()

