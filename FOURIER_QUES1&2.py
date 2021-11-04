#2019phy1073
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy as sp
import scipy.integrate as spi
import random
import math as mt


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


l1=0  # for parabola 0
l2=np.pi  # for parabola pi
hp=(abs(l2)-abs(l1))/2
#####################################################################
k=20    # NO. OF TERMS OF AN AND BN 
period=2  #no. of periods
###################         FUNCTIONS     ########################################
# def func(x):                                  
#     return x
#########
def func(x):              #PARABOLA FUNCTION
    return x*(x-(np.pi))
#########

def a_n(x,n):
    f_of_x=func(x)
    return f_of_x * np.cos((x*n*mt.pi)/hp) 

def b_n(x,n):
    f_of_x=func(x)
    return f_of_x * np.sin((x*n*mt.pi)/hp)
    


####################
x=np.linspace(l1,period*l2,100)

####################
z=spi.simps(func(np.linspace(l1,l2)),np.linspace(l1,l2),dx=0.01)


a0=z/hp
a_0= a0/2
b0=0
print("a0 =" + str(a0))
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

for i in range (1,k):
    integrand=a_n(x,i) 
    y1=func(x)*np.cos((i*x*mt.pi)/hp)  # AN USING INBUILT
    y2=spi.simps(y1,x,dx=0.01)
    an_inbuilt.append(y2)
 #############################
    y5=func(x)*np.sin((i*x*mt.pi)/hp)  # BN USING INBUILT
    y6=spi.simps(y5,x,dx=0.01) 
    bn_inbuilt.append(y6)
  #############################
    y3=user_simpson(a_n,l1,l2,50,i)  # AN USING USER SIMPSONN
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
# print(n,bn_inbuilt)
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
plt.plot(x,y,'g--',label=" inbuilt simpson integration")
plt.plot(x,series_of_sum_trapezoidal,'b:',label="user trapezoidal integration")
plt.plot(x,series_of_sum_simpson,'r*',label="user simpson integration")
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
# print(an_inbuilt)
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
                OUTPUT   QUESTION 2 WITH GRAPH

a0 =-3.289840170278713
    a_n_inbilt  an_user_simpson  an_user_trapezoidal  bn_inbuilt
0     3.140832         1.570790             1.570817   -9.869298
1     0.782242         0.392673             0.392720   -4.933960
2     0.341960         0.174472             0.174554   -3.288034
3     0.183804         0.098065             0.098195   -2.463898
4     0.106276         0.062655             0.062853   -1.967862
5     0.059746         0.043369             0.043654   -1.635233
6     0.027302         0.031682             0.032078   -1.395330
7     0.001974         0.024028             0.024564   -1.212738
8    -0.019471         0.018699             0.019413   -1.067726
9    -0.038643         0.014791             0.015729   -0.948420
10   -0.056368         0.011781             0.013002   -0.847237
11   -0.073056         0.009344             0.010929   -0.759112
12   -0.088891         0.007259             0.009315   -0.680535
13   -0.103918         0.005358             0.008035   -0.609009
14   -0.118107         0.003492             0.007002   -0.542713
15   -0.131378         0.001499             0.006157   -0.480304
16   -0.143627        -0.000831             0.005456   -0.420774
17   -0.154737        -0.003828             0.004869   -0.363367
18   -0.164589        -0.008080             0.004372   -0.307512
           X  s_n_inbilt  sn_user_simpson  sn_user_trapezoidal
0   0.000000    1.906432         0.812296             0.858793
1   0.063467  -12.585381         0.675634             0.655270
2   0.126933  -11.752875         0.339585             0.323217
3   0.190400  -11.688503         0.037431             0.064411
4   0.253866  -12.173663        -0.194952            -0.214909
..       ...         ...              ...                  ...
95  6.029319   14.992808        -0.194952            -0.214909
96  6.092786   14.286733         0.037431             0.064411
97  6.156252   17.824973         0.339585             0.323217
98  6.219719   19.083370         0.675634             0.655270
99  6.283185    1.906432         0.812296             0.858793

[100 rows x 4 columns]
"""