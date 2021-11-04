#2019PHY1073
import sympy as sm
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import legendre
def fact(n): 
    if n < 0: 
        return 0
    elif n == 0 or n == 1: 
        return 1
    else: 
        fact = 1
        while(n > 1): 
            fact *= n 
            n -= 1
        
        return fact


def user_legendre(k,x):
    y=[]
    values=0
    if (k%2==0):
        M=int((k/2))
    else:
           M=int((k-1)/2)
    s=0
    p_coeff=[]
    
    
        
        
    for i in range(M+1):
        s=((-1)**i*fact(2*k-2*i))/((2)**k*fact(i)*fact(k-i)*fact(k-2*i))
        p_coeff.append(s)
        val = p_coeff[i]*(x**(k-(2*i)))
        values += val
    return values
        
        
   

k=int(input(print("enter the order of legendre polynomial :")))
x=np.linspace(-1,1 ,100)   
x_plot=[]
y_plot=[]
x_plot.append(x)
p_1=[["x",x_plot]]
for j in range(0,k):
    y=user_legendre(j,x)
    y_plot.append(y)
    p_2=[['P('+str(j)+')',y_plot[j]]]
    p_1.append(p_2)
    plt.plot(x,y_plot[j],label='P('+str(j)+')')
plt.xlim(-1, 1)
plt.ylim(-1, 1.2)
plt.grid()
plt.legend()
plt.title("THE USER MADE LEGENDRE POLYNOMIAL")
plt.xlabel("X ")
plt.ylabel("P(x) ")
plt.show()
print(x)
print(y_plot)
# array = np.array([x, y_plot])
# np.savetxt("legendrepolynomial.txt",array)
def inbuiltlegendre(x,n) :
          
    leg=legendre(n)
    p_x=leg(x)
    return p_x
def inbuiltlegendreplot():
    for i in range (0,k):
        func=inbuiltlegendre(x,i)
        plt.plot(x,func,label="n="+str(i))
    plt.title('THE INBUILT LEGENDRE POLYNOMIAL')
    plt.grid(True)
    plt.legend()
    plt.show()     
    
    

inbuiltlegendreplot()



o=int(input("Enter the no. of terms upto which you want to verify the recurrance relation : "))
print("_______________________________First relation______________________________")    
if k>=2:    
   
    rhs_1=k*(user_legendre(k,x))
    lhs_1=((2*k-1)*x*user_legendre(k-1, x))-((k-1)*user_legendre(k-2,x))
  #  print(f'{rhs_1:20}',f'{lhs_1:20}')
    print("Recurrance relation first :-")
    print(f'{" ":10}',f'{"R.H.S.":20}',f'{"|":15}',f'{"L.H.S.":20}')
    for i in range(o+1):
        print(f'{rhs_1[i]:20}',f'{"           |":5}',f'{lhs_1[i]:15}')
else:
    print("Can't find the relation for value of K less than 2.")
print("_______________________________Second relation______________________________")
if k>=1:
    
    rhs_2=((k+1)*user_legendre(k+1,x))
    lhs_2=((2*k+1)*x*user_legendre(k,x))-(k*user_legendre(k-1,x))
    print("Recurrance relation second :-")
    print(f'{" ":10}',f'{"R.H.S.":20}',f'{"|":5}',f'{"L.H.S.":15}')
    for i in range(o+1):
        print(f'{rhs_2[i]:20}',f'{"           |":5}',f'{lhs_2[i]:20}')
else:
    print("Can't find the relation for value of K less than 1.")
print(x)