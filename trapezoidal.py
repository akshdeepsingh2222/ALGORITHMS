import numpy as np
import scipy
def user_trapezoidal (f,a,b,n):
    h=(b-a)/n
    x=[a]
    y=[]
    total=0
    for i in range (1,n):
        step=x[0]+i*h
        x.append(step)
    x.append(b)
    q=len(x)
    u=q-1
    for j in range(1,q-1):
        z=f(x[j])
        y.append(z)
    for ele in range(0, len(y)): 
        total = total +y[ele]
    total=2*(total) + f(x[0]) + f(x[u])
    res=(h/2)*total
        
    #res=(h/2)*(f(b)+f(a))
       
    
   
    return res
    # print(y)
    # print(x)
    # print(total)
    # print(h)
   
x1=np.linspace(0,1,100)
y=[]        

def f1(x):
    return x*(2.718)**x

def f2(x):
    return 1/(x*x+6*x+10)

y.append(f1(x1))
z=np.trapz(y, x1, dx=1.0, axis=- 1)
print("value of integral using inbuilt function :")
print(z)

# error=user_trapezoidal(f1,0,3,8)-trapz(f1,0,3,8)
print("the computed value of integrated function:")
h=user_trapezoidal(f1,0.5,2.5,8)
print(h)
error=h-z
print("THE TRUNKATION ERROR IS  ")
print(error)




 



