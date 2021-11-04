import numpy as np
from scipy import integrate
def user_simpson(func,a,b,n):
    h = float((b-a)/n)
    result = (1)*(func(a)+func(b))
    for i in range(1,n,2):
        result+= 4*(func(a+i*h))
    for j in range(2,n-1,2):
        result+=2*(func(a+j*h))
    result*=h/3
    return result
    
    
def f1(x):
    return x


def f2(x):
    return 1/(1-x)



#user_simpson(f2,0,1,2)
h=user_simpson(f1,0,2,4)

# print(f2(0.25))
# print(f2(0.5))
# print(f2(0.75))

x1=[0,0.25,0.50,0.75,1]
y=[1,0.8000,0.6667,0.5714,0.5000]

z=integrate.simps(y, x1)


print("THE INTEGRATED VALUE USING USER DEFINED FUNCTION :")
print(h)
print("THE INTEGRATED VALUE USING INBUILT FUCTION :")
print(z)
error=float(h)-float(z)
print("THE TRUNKATION ERROR IS  ")
print(error)

