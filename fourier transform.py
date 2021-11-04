import numpy as np
import math as mt
import cmath 
import matplotlib.pyplot as plt


N=4
i=complex(0,1)

def f(t):
    return 5 + 2*np.cos(2*np.pi*t - np.pi/2) + 3*np.cos(4*np.pi*t - np.pi/2)
# def f(k):
#     return 5 + 2*np.cos((2*mt.pi*k)/N-mt.pi/2)+3*np.cos((4*mt.pi*k)/N) # k=t*N


t0=1     # t0 is fundamental time period
T=t0/N   # T is step size
func_values=[]
for k in range(N):
    func_values.append(f((t0/N)*k))
func_values=np.array(func_values,"float")
print(func_values)
func_transpose=np.transpose([func_values])
def twid_matrix(N):
    iota=complex(0,1)
    W=cmath.exp((iota*2*np.pi)/N)  # this W is TWIDDLE FACTOR
    twid_matrix = []
    for i in range(0,N,1):
        row_list = [] 
        for j in range(0,N,1):
            row_list.append(W**(i*j))
        twid_matrix.append(row_list)         
    twid_matrix = (np.matrix(twid_matrix))
    return twid_matrix

# i =  complex(0,1) 
# w = cmath.exp(i*2*(np.pi)/N) 
# print(twid_matrix(N))
# print(twiddle(w,N))
DFT=np.matmul(twid_matrix(N),func_transpose)
print(DFT)
# def inverse_DFT(N):    # it gives original data points
#     twid_inverse=np.transpose(twid_matrix(N))
#     original_data_pts=np.matmul(twid_inverse,DFT)
#     return original_data_pts
# print("original",func_values)
# print("transformed",inverse_DFT(N))

# print(twid_mat(4))
# k=np.arange(0,N,1)      # k is from 0 to N-1
# t=np.arange(0,t0,T)   # T is the step size
# fk=f(t)
# #ft=f1(t)
# w=(2*mt.pi/N*T)*k    # this is list of all omegas(w) with  N-1 terms {n*w_o}[2pi/NT]
# f_wn=[]
# W=cmath.exp((-i*2*cmath.pi)/N)  # this W is TWIDDLE FACTOR
# #W=1

# for j in range(0,N):       # for different values of N (F_k)
#     f_n=0
#     for n in range(0,N):   # sums upto N-1  for constant N
#         f_n+=fk[n]*(W**(j*n))
        
#     # f_n=T*f_n
#     f_wn.append(f_n)
# # f_wn=np.array(f_wn)           # converting f_wn list to npumpy array
# # f_wn=np.vstack(f_wn)          # transpose of f_wn

# ## twiddle matrix
# twid_row1=[]

# twid_column=[]

# twid_mat=np.empty([N,N],dtype=complex)
# for i in range(0,N):
#     k=0
#     val=W**k
#     twid_row1.append(val)
#     # np.append(twid_mat,val)
# twid_column.append(twid_row1)
# for j in range(1,N):               #column
#     twid_elements=[]    
#     for k in range(0,N):                  
#         val=W**k
#         twid_elements.append(val)
#     twid_column.append(twid_elements)   
    
# # np.append(twid_mat,twid_row1)
# # np.append(twid_mat,twid_column)
# print(twid_column)




# # conjugate_f_wn=np.conj(f_wn)  # finding conjugate of array of f_wn
# # F_N=[]
# # for i in range(0,len(f_wn)):  
# #     mul=f_wn[i]*conjugate_f_wn[i]
# #     mul=cmath.sqrt(mul)
# #     mul=abs(mul)              # absolute
# #     F_N.append(mul) 
# # print(F_N) 
# t=np.arange(0,t0,0.01)           
# plt.plot(t,f(t),'-*')  
# # plt.plot(k/N,f(k),'r')
# plt.grid(True)
# plt.title("f(t) v/s t")
# plt.show()  
# plt.plot(w,f_wn,'-*')
# plt.title("F(w) v/s w")
# plt.grid(True)
# plt.show()
    

# plt.plot(k,fk)
# plt.scatter(w/2*mt.pi,f_wn)
# plt.bar(k,F_N,0.01)
# plt.grid()
# plt.show()
