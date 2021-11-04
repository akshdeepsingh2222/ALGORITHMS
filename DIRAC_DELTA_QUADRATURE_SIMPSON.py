import matplotlib.pyplot as plt
from scipy import integrate
from sympy import *
import numpy as np
from scipy import integrate
from ipywidgets import interactive

def f(x,sigma):
    h=(x-2)**2/(2*sigma**2)
    integrand=((x+3)*np.exp(-h))/np.sqrt(2*np.pi*sigma)
    return integrand 

def user_simpson(func,a,b,n,sigma):
    h = float((b-a)/n)
    result = (1)*(func(a,sigma)+func(b,sigma))
    for i in range(1,n,2):
        result+= 4*(func(a+i*h,sigma))
    for j in range(2,n-1,2):
        result+=2*(func(a+j*h,sigma))
    result*=h/3
    return result

def plotting(sigma=1):
    n=np.arange(2,10,2)
    x=np.linspace(-sigma,sigma,100)
    y=[]
    for i in range(len(x)):
        y1=f(x[i],sigma)
        y.append(y1)
    z=integrate.simps(y, x,dx=0.01)
    z1=user_simpson(f,-sigma,sigma,10,sigma)   
    error=np.abs(z-z1) 
    print("INTEGRATION BY INBUILT SIMPSON METHOD : " ,z)
    print("INTEGRATION BY USER-MADE SIMPSON METHOD : " ,z1)
    print("ERROR IN INTEGRATION : ",error)
    plt.plot(y,x,label="SIGMA="+str(sigma))
    plt.xlabel("X- values")
    plt.ylabel("F(x)")
    plt.title("DIRAC DELTA FUNCTION")
    plt.legend()
    plt.grid(True)
    plt.show()
    return 
plotting(0.1)
interactive_plot=interactive(plotting,sigma={1,0.1,0.01})
interactive_plot
