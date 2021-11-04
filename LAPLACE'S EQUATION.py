import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
# from ipywidgets import interactive

N=5
# a=10
# b=5
# hx=N/a
# hy=N/b
def boundary_north(choice):
    if choice==1:
        y=[]
        for i in range(N):
            y.append(100)
        return y
    if choice==2:
        y=[]
        for i in range(N):
            y.append(100)
        return y
    if choice==3:
        a=20
        x=np.arange(N)
        t0=150
        y=t0*np.cos((np.pi*x)/a)
        return y
        
    if choice==4:
        x=np.arange(N+1)
        a=20
        y=np.sin((np.pi*x)/a)
        return y  
    if choice==5:
        u0=1
        a=2
        # x=np.arange(N)
        x=(a/N)*np.arange(N)
        y=u0*x*(a-x)
        return y
def boundary_west(choice):
    if choice==1:
        return 75
    if choice==2:
        return 0  
    if choice==3:
        return 0
    if choice==4:
        return 0
    if choice==5:
        return 0
def boundary_east(choice):
    if choice==1:
        return 50
    if choice==2:
        return 0 
    if choice==3:
        return 0
    if choice==4:
        return 0
    if choice==5:
        return 0
def boundary_south(choice):
    if choice==1:
        return 0
    if choice==2:
        return 0
    if choice==3:
        return 0
    if choice==4:
        return 0
    if choice==5:
        return 0
wn=[]
choice=int(input(print("enter your choice :")))
y=boundary_north(choice)
for row in range(N):
    for col in range(N):
        if row==0 and col<=N-1:
            wn.append(y[col])    
        elif row==N-1 and 1<=col<=N-2: #previous (and col<=1 ); then (1<=col<=N-1)
            wn.append(boundary_south(choice))
        elif col==0 and row <=N-1:
            wn.append(boundary_west(choice))
        elif col==N-1 and row<=N-1:
            wn.append(boundary_east(choice))
        else :
            wn.append(0)
wn=np.array(wn)
wn=np.reshape(wn,(N,N))
wn=wn.astype('float64')
wm=np.array(wn)
wm=wm.astype('float64')
print("0'th iterated matrix")
print(wn)
count=0
run=True
while True:
    wn=np.array(wn)
    for i in range(1,N-1):
          for j in range(1,N-1):
              wn[i][j]=0.25*(wn[i+1][j]+wn[i-1][j]+wn[i][j+1]+wn[i][j-1])   #wn is k_lower 
              wm[i][j]=0.25*(wm[i+1][j]+wm[i-1][j]+wm[i][j+1]+wm[i][j-1])    #wm is K_higher
    k=1
    for k in range(k): 
        for i in range(1,N-1):
            for j in range(1,N-1):
                wm[i][j]=0.25*(wm[i+1][j]+wm[i-1][j]+wm[i][j+1]+wm[i][j-1])    #wm is K_higher
    W_diff=np.array(wm-wn)
    W_diff=np.abs(W_diff)
    W_max=float(np.max(W_diff))
    count+=1 
    if W_max<0.001:
        print("no. of iterations is",count)
        break          
print("K times iterated matrix")
print(wn)
print("K+1 times iterated matrix")
print(wm)
print("difference matrix")
print(W_diff)
print("max diff :",W_max)    
x,y=np.meshgrid(np.arange(N),np.arange(N))
ax=plt.axes(projection="3d")
plt.contourf(x,y,wm)
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Temprature(C)')
ax.plot_surface(x,y,wm)
plt.colorbar()
plt.show()

# ax=plt.axes(projection="3d")
# plt.contourf(x,y,wm)
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.set_zlabel('Temprature(C)')
# plt.colorbar()
# plt.show()

plt.matshow(wm)
plt.colorbar()
plt.show()