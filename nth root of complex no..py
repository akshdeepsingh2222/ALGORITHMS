import math as mt 


x=3 
y=-1
n=4
r=mt.sqrt(x**2+y**2)
theta=mt.acos(x/r)
w=[]
for i in range(n):
    w1=r**(1/n)*complex(mt.cos((theta+2*i*mt.pi)/n),mt.sin((theta+2*i*mt.pi)/n))
    w.append(w1)
print(r**(1/n))
print(w)
print(theta)
