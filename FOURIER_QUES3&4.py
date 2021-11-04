#2019phy1073
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp
import scipy.integrate as spi
import random
import math as mt
from scipy.integrate import quad 

def user_simpson(func2,a,b,n,n1):
    h = float((b-a)/n)
    result = (1/3)*(func2(a,n1)+func2(b,n1))
    for i in range(1,n,2):
        result+= 4*(func2(a+i*h,n1))
    for j in range(2,n-1,2):
        result+=2*(func2(a+j*h,n1))
    result*=h/3
    return result

def user_trapezoidal (f,a,b,n,n1):
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
        z=f(x[j],n1)
        y.append(z)
    for ele in range(0, len(y)): 
        total = total +y[ele]
    total=2*(total) + f(x[0],n1) + f(x[u],n1)
    res=(h/2)*total
    return res


l1=-2  # for rectifier is l1 -2
l2=+2  # for rectifier l2 is +2
hp=(abs(l2)-(l1))/2
#####################################################################
k=20    # NO. OF TERMS OF AN AND BN 
period=3 #no. of periods
###################         FUNCTIONS     ########################################

# def func(x):                  # SQUARE FUNCTION
#     if(x >= l1 and x< 0):
#         return -1
#     elif(x>0 and x<= l2):
#         return 1 
#     elif(x==0):
#         return 0
####################
def func(x):                     #HALF WAVE RECTIFIER
        if(x >= l1 and x< 0): 
            return 0
        elif(x>0 and x<= l2):
            return np.sin((np.pi*x)/hp)
        elif(x==0):
            return 0
# ####################
def a_n(x,n):
    f_of_x=func(x)
    return f_of_x * np.cos((x*n*mt.pi)/hp) 

def b_n(x,n):
    f_of_x=func(x)
    return f_of_x * np.sin((x*n*mt.pi)/hp)
    


choice=1
####################
x=np.linspace(l1,period*l2,100)
def inbuiltquadrature(a,b,f,n):
    res = quad(f,a,b,args=(n))
    return res[0]/hp

###################   
G=[]
C=np.linspace(l1,l2,100)
for i in range(0,len(C)):
    g=func(C[i])
    G.append(g)
z=spi.simps(G,C,dx=0.01)
####################


a0=z/hp
a_0= a0/2
b0=0

            #####################################################

an_inbuilt=[]
an_user_simpson=[]
an_user_trapezoidal=[]
y=[]
series_of_sum_simpson=[]
series_of_sum_trapezoidal=[]

f_of_x=[]
bn_inbuilt=[]
bn_user_simpson=[]
bn_user_trapezoidal=[]
x4=np.linspace(l1,l2,100)
u1=[]
u2=[]
for i in range (1,k):
    y7=inbuiltquadrature(l1,l2,a_n,i)   # AN USING INBUILT
    an_inbuilt.append(y7)
  #############################
    y6=inbuiltquadrature(l1,l2,b_n,i)   # BN USING INBUILT
    bn_inbuilt.append(y6)
  #############################
    y3=user_simpson(a_n,l1,l2,50,i)  # AN USING USER SIMPSON
    an_user_simpson.append(y3)
  ################################
    y4=user_trapezoidal(a_n,l1,l2,500,i)   # AN USING USER TRAPEZOIDAL
    an_user_trapezoidal.append(y4)
  ################################
    y9=user_simpson(b_n,l1,l2,10,i)      # BN USING USER SIMPSON
    bn_user_simpson.append(y9)
################################
    y13=user_trapezoidal(b_n,l1,l2,10,i)  # BN USING USER TRAPEZOIDAL
    bn_user_trapezoidal.append(y13)
################################
n=np.arange(1,k,1)

plt.plot(n,an_inbuilt,"r*-",label="an_inbuilt")
plt.plot(n,bn_inbuilt,"b*-",label="bn_inbuit")
plt.plot(n,an_user_simpson,label="an_user_simpson")
plt.plot(n,bn_user_simpson,label="bn_user_simpson")
plt.plot(n,an_user_trapezoidal,label="an_user_trapezoidal")
plt.plot(n,bn_user_trapezoidal,label="bn_user_trapezoidal")
plt.xlabel("n")
plt.ylabel("an and bn")
plt.title("an and bn V/S n")
plt.legend()
plt.grid(True)
plt.show()

##################################
x=np.linspace(l1,period*l2,100)
for j in range(0, len(x)):
    x1=x[j]
    series_sum=a_0
    series_sum_simpson=a_0
    series_sum_trapezoidal=a_0
    for q in range (1,k):
        series_sum+=an_inbuilt[q-1]*mt.cos((q*x1*mt.pi)/hp) + bn_inbuilt[q-1]*mt.sin((q*x1*mt.pi)/hp)
        series_sum_simpson+=an_user_simpson[q-1]*mt.cos((q*x1*mt.pi)/hp) +bn_user_simpson[q-1]*mt.sin((q*x1*mt.pi)/hp)
        series_sum_trapezoidal += an_user_trapezoidal[q-1]*mt.cos((q*x1*mt.pi)/hp) + bn_user_trapezoidal[q-1]*mt.sin((q*x1*mt.pi)/hp)
    y.append(series_sum)
    series_of_sum_simpson.append(series_sum_simpson)
    series_of_sum_trapezoidal.append(series_sum_trapezoidal)
# #####################################################################

list_of_tuples = list(zip(an_inbuilt,an_user_simpson,an_user_trapezoidal,bn_inbuilt))
list_of_tuples
df = pd.DataFrame(list_of_tuples,columns=["a_n_inbilt","an_user_simpson","an_user_trapezoidal","bn_inbuilt"])
print(df)
list_of_tuples = list(zip(x,y,series_of_sum_simpson,series_of_sum_trapezoidal))
list_of_tuples
df = pd.DataFrame(list_of_tuples,columns=["X","s_n_inbilt","sn_user_simpson","sn_user_trapezoidal"])
print(df)

# ################################################
plt.plot(x,y,'r--',label=" inbuilt simpson integration")
plt.plot(x,series_of_sum_trapezoidal,'b:',label="user trapezoidal integration")
plt.plot(x,series_of_sum_simpson,'g',label="user simpson integration")
plt.title("FUNCTION PLOT")
plt.xlabel("X")
plt.ylabel("sum of series - F(x)")
plt.grid(True)
plt.legend()
plt.show()

x=np.arange(1,k,1)

######################################
an_2=np.array(an_inbuilt)
bn_2=np.array(bn_inbuilt)
an_2_simps=np.array(an_user_simpson)
bn_2_simps=np.array(bn_user_simpson)
an_2_trap=np.array(an_user_trapezoidal)
bn_2_trap=np.array(bn_user_trapezoidal)
anbnsum=np.array(np.square(an_2)+np.square(bn_2))
anbnsimpssum=np.array(np.square(an_2_simps)+np.square(bn_2_simps))
anbntrapsum=np.array(np.square(an_2_trap)+np.square(bn_2_trap))
plt.plot(x,anbnsum,'r*-',label='inbuilt')
plt.plot(x,anbnsimpssum,'g:',label='user simpson')
plt.plot(x,anbntrapsum,'b',label='user trap')
plt.xlabel('n')
plt.ylabel('an^2+bn^2')
plt.title('an^2+bn^2 V/S n')
plt.legend()
plt.grid(True)
plt.show()

#####################################
an_inbuilt=[x/an_inbuilt[0] for x in an_inbuilt]
bn_inbuilt=[x/bn_inbuilt[0] for x in bn_inbuilt]
an_user_simpson=[x/an_inbuilt[0] for x in an_user_simpson]
bn_user_simpson=[x/bn_inbuilt[0] for x in bn_user_simpson]
an_user_trapezoidal=[x/an_inbuilt[0] for x in an_user_trapezoidal]
bn_user_trapezoidal=[x/bn_inbuilt[0] for x in bn_user_trapezoidal]

plt.plot(x,an_inbuilt,"r.-",label="an_inbuilt")
plt.plot(x,bn_inbuilt,"b.-",label="bn_inbuilt")
plt.plot(x,an_user_simpson,label="an_user_simpson")
plt.plot(x,bn_user_simpson,label="bn_user_simpson")
plt.plot(x,an_user_trapezoidal,label="an_user_trapezoidal")
plt.plot(x,bn_user_trapezoidal,label="bn_user_trapezoidal")
plt.xlabel("n")
plt.ylabel("an/a1 and bn/b1")
plt.title("an/a1 and bn/b1 V/S n")
plt.legend()
plt.grid()
plt.show()
###########################################
"""
                    OUTPUT OF QUES-4  WITH GRAPH 
      a_n_inbilt  an_user_simpson  an_user_trapezoidal    bn_inbuilt
0   1.184856e-18        -0.001684        -9.454910e-18  5.000000e-01
1  -2.122066e-01        -0.426094        -4.244299e-01  2.402212e-17
2  -2.041599e-17        -0.001739        -1.249415e-17 -6.938894e-17
3  -4.244132e-02        -0.086580        -8.489939e-02 -4.132689e-17
4   7.553815e-17        -0.001856        -5.025561e-17  2.081668e-17
5  -1.818914e-02        -0.038103        -3.639503e-02  3.365532e-17
6  -4.906237e-17        -0.002051        -6.711201e-16  3.816392e-17
7  -1.010508e-02        -0.021974        -2.022692e-02 -1.057859e-16
8   1.033075e-17        -0.002357         7.443588e-16  3.816392e-17
9  -6.430503e-03        -0.014678        -1.287777e-02 -1.278392e-16
10 -1.615699e-16        -0.002834         7.664289e-17  0.000000e+00
11 -4.451887e-03        -0.010789        -8.920547e-03  9.404361e-17
12  1.664470e-16        -0.003596        -1.577668e-16  6.418477e-17
13 -3.264717e-03        -0.008501        -6.546215e-03 -1.862473e-16
14 -2.969877e-16        -0.004893         9.275733e-16  3.348016e-16
15 -2.496548e-03        -0.007072        -5.009885e-03  1.162975e-16
16 -9.745616e-17        -0.007325        -6.400200e-16 -1.162265e-16
17 -1.970959e-03        -0.006153        -3.958716e-03 -2.394471e-16
18 -5.509334e-16        -0.012701        -2.640429e-16  3.642919e-17
           X  s_n_inbilt  sn_user_simpson  sn_user_trapezoidal
0  -2.000000    0.016780        -0.260572            -0.284928
1  -1.919192   -0.004997         0.292613             0.263227
2  -1.838384    0.003083        -0.979091            -0.892237
3  -1.757576   -0.001440         1.218997             1.077643
4  -1.676768    0.000106        -1.087516            -0.885458
..       ...         ...              ...                  ...
95  5.676768    0.486302         1.400691             1.221625
96  5.757576    0.370223        -1.090220            -0.976770
97  5.838384    0.254231         0.847216             0.770277
98  5.919192    0.121595        -0.718105            -0.666775
99  6.000000    0.016780        -0.260572            -0.284928

[100 rows x 4 columns]

"""

