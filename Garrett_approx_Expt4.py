import numpy as np
import matplotlib.pyplot as plt
import math as mt
import pandas as pd
#line no 94,93
def garrett_approx(state,precision,Mass,Length,v):
    h_bar=1.054571e-34       ####J s
    def Energy(state,l,V):
        energy=(state*np.pi*h_bar)**2/(2*Mass*l**2)
        return energy  
    def New_len(en,v):
        delta=(h_bar/(2*np.pi))/(abs((2*Mass*(v-en))**0.5))
        return delta    # L=l+2*delta
    E=[]
    Len=[]
    deltaa=0
    count=0
    run=True
    # while run:
    #     Len.append(k)
    #     j=Energy(state,Length+2*deltaa, v)  #E_(i-1)||lower iteration energy value
    #     k=New_len(j,v)
    #     j1=Energy(state, Length+ 2*k, v)      #E_(i)|| higher iteration energy value
    #     E.append(j)
    #     deltaa=k
    #     count+=1
    #     # E=np.around(E,precision)
    #     # Len=np.around(Len,precision)
    
    #     if j==j1:
    #         break 
    for i in range(1,200):
        j=Energy(state,Length+2*deltaa, v)  #E_(i-1)||lower iteration energy value
        k=New_len(j,v)
        E.append(j)
        Len.append(k)
        deltaa=k
        count+=1
    E = [element* 6.242e+18  for element in E]     # converting from joule to eV
    Len=[(element+Length)/1e-10 for element in Len] # converting from mtr to angstrom
    # E=np.around(E,precision)
    Len=np.around(Len,precision)

    count_l=0
    while True:
        X=E[count_l+1]
        Y=E[count_l]
        count_l+=1
        if X==Y:
            break
    # # count_l=30
    return E,Len,count_l
Length=1e-10          # in mtr
me=9.10938e-31 # in Kg
# mp=1.6726219e-27
v=150e-19                    #in V
state=1
precision=5
E,Len,iteration=garrett_approx(state, precision, me, Length,v)
E_o=E[:iteration+1]
# print(iteration , E_o)
# print(Len_0)
df = pd.DataFrame(list(zip(E_o,Len)),columns =['Energy(in eV)','Length(in A^o)'])
print(df)
plt.plot(np.linspace(0,len(E_o),len(E_o)),E_o,'-*',label="STATE="+str(state)) ,plt.xlabel("iterations"),plt.ylabel("Energy(in eV)"), plt.legend()
plt.minorticks_on() ,plt.grid(b=True,which="both",axis="both")
plt.show()
dell=[(element-Length/1e-10) for element in Len]
dell=dell[:iteration]
plt.plot(np.linspace(0,len(dell),len(dell)),dell,'g-*')
plt.title("delta v/s iteration"), plt.grid()
plt.show()
state_list=[]
En=[]
iteration_list=[]
converg=[]
converg_l=[]
K=8
for i in range(1,K):
    E,Len,iteration=garrett_approx(i, precision, me, Length,v)
    # plt.plot()
    iteration_list.append(iteration)
    converg.append(E[iteration-1])
    converg_l.append(Len[iteration])
    En.append(E[0])
    state_list.append(i)
df = pd.DataFrame(list(zip(state_list,En,converg,converg_l,iteration_list)),index=None,columns =['States','Energy(in eV)','E_convergence',"L_convergence",'Iteration'])
print(df)
#################################################
V=np.arange(2,75,8)*(1e-19)*1e2
# V=V*(1e-19)         # V
E_new=[]
L_new=[]
Del_new=[]
for i in range(len(V)):
    E,L,C=garrett_approx(1, precision, me, Length, V[i])
    delta=[(element-Length/1e-10) for element in L]
    E_new.append(E[C])
    L_new.append(L[C])
    Del_new.append(delta[C])
V=[(element/1e-19) for element in V]  
plt.plot(V,E_new) , plt.title("Energy V/s V ")
plt.grid()
plt.show()
###################################################
potential=100e-19
L=np.linspace(1,50,50)
L1=[(element*1e-10) for element in L]
energy=[]
for i in range(len(L)):
    E,Len,C=garrett_approx(2,5, me, L[i], potential)
    E_new=E[C]
    energy.append(E_new)
plt.plot(L,energy)
plt.title('energy V/S Len')
plt.grid()
plt.show()