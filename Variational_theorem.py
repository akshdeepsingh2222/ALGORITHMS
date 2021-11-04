import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from scipy.integrate import quad
lamb=1
a=0*lamb  
b=1*lamb
l=b-a
x=Symbol('x')
def user_simpson(func,a,b,n):
    h1 = float((b-a)/n)
    result = (1)*(func.subs(x,a)+func.subs(x,b))
    for i in range(1,n,2):
        result+= 4*(func.subs(x,(a+i*h1)))
    for j in range(2,n-1,2):
        result+=2*(func.subs(x,(a+j*h1)))
    result*=h1/3
    return result
   
def H(alpha):
    N=150
    h=(b-a)/(N+1)
    psi=x*((l-x)**alpha)
    # psi=(b-x)*pow((x),alpha)      # trial wave function
    def p(alpha):
        return (x)*pow((b-x),alpha) 
    psi=p(alpha)
    differential_2=-1*diff(diff(psi,x),x)
    f_x=differential_2*psi
    g_x=psi*psi
    N=user_simpson(f_x, a+a*h,b- b*h, 8)
    D=user_simpson(g_x, a+a*h, b-b*h, 8)
    H_star=N/D
    return H_star
print("hdcg\s",H(2))

alpha=np.around(np.linspace( 0.55,1.55,20),2)
h=[]
for i in range(len(alpha)):    
    h.append(float(H(alpha[i])))
# print(h)
plt.plot(alpha,h,"-*")
# plt.bar(alpha,h,width=0.3,alpha=0.75)
plt.grid()
plt.title("Ev/s alpha")
plt.xlabel("alpha")
plt.ylabel("H star")
plt.show()
# print(h,alpha)
E0=min(h)
pos=h.index(min(h))
print("index of min energy ::",pos)
epsilon=500   # KeV
print("min energy :",E0)
l=alpha[pos]
print("alpha of min. energy::",l)
h_real=[]
for i in range(len(h)):
    h_temp=(h[i]*epsilon)/(1.708*8*np.pi**2)   # in kev   #1.708
    h_real.append(h_temp)
print("pos of min energy :",h_real.index(min(h_real)))  
g=min(h_real)
print("min energy(in KeV):",g)   
alpha=np.around(np.linspace( 0.6,1.5,10),2)
x2=np.linspace(a,b,100)
all_psi=[]
for j in range(len(alpha)):
    psi=(b-x)*((x)**alpha[j])      # trial wave function
    psi1=(b-x)*((x)**l) 
    psi_temp=[]
    psi11=[]
    for i in range(len(x2)):
        psi_temp.append(psi.subs(x,x2[i]))
        psi11.append(psi1.subs(x,x2[i]))
    all_psi.append(psi_temp)
for k in range(len(alpha)):
    plt.plot(x2,all_psi[k],label=str(alpha[k]))
plt.plot(x2,psi11,'*',label="min energy plot")
plt.title("psi V/S x for diff alpha")
plt.xlabel("x")
plt.ylabel("psi")
# plt.legend()
plt.grid()
plt.show()
"""   OUTPUT
index of min energy :: 9
min energy : 9.980995245312103
alpha of min. energy:: 1.02
pos of min energy : 9
min energy(in KeV): 37.005498103906895
"""