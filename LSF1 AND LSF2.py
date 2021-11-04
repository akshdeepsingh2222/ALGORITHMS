import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import math as mt

def lsf1(x,y):
    
    y_star=[]
    dev_sq=[]
    sigma_y_star=[]
    for i in range(0,len(x)):
        temp_y=sum(y[i])/len(y[i])
        y_star.append(temp_y) 
    for j in range(0,len(x)):
        m=[]
        for k in range(0,len(y[j])):
            s=(y[j][k]-y_star[j])**2
            m.append(s)
        sum_m=sum(m)    
        dev_sq.append(sum_m)
    
    for t in range(0,len(x)):
        if len(y[t])>2:
            l=2
        elif len(y[t])==2 :
            l=1
        elif len(y[t])<2 :
            l=0
        temp=(dev_sq[t]/(int(len(y[t])-l)))**0.5
        sigma_y_star.append(temp)

    n=len(x)
    x=np.array(x)
    y_star=np.array(y_star)
    x_sq=np.square(x)
    y_sq=np.square(y_star)
    x_ystar=x*y_star
    x_sum=np.sum(x)
    ystar_sum=np.sum(y_star)
    x_sq_sum=np.sum(x_sq)
    ystar_sq_sum=np.sum(y_sq)
    x_y_sum=np.sum(x_ystar)    
    delta_o=ystar_sum*x_sq_sum-x_sum*x_y_sum
    delta=n*x_sq_sum-(x_sum)**2
    delta1=n*x_y_sum-x_sum*ystar_sum
    A=delta_o/delta                              # Y-INTERCEPT
    B=delta1/delta                               # SLOPE
    sigma=1
    sigma_B=sigma*mt.sqrt(n/delta)              # SIGMA B
    sigma_A=sigma*mt.sqrt(x_sq_sum/delta)       # SIGMA A
    y_fitted=B*x+A                                # Y FITTED
    y_plus_sigma=(B+sigma_B)*x+(A+sigma_A)      # Y +SIGMA
    y_minus_sigma=(B-sigma_B)*x+(A-sigma_A)     # Y -SIGMA
    
    plt.minorticks_on()
    plt.scatter(x,y_star,label="mean")
    plt.errorbar(x,y_star,sigma_y_star,color="red",fmt="o")
    plt.plot(x,y_fitted,label="fitted line")
    plt.plot(x,y_plus_sigma,label="Y+sigma")
    plt.plot(x,y_minus_sigma,label="Y-sigma")
    plt.xlabel("wavelength(in nano metre)")
    plt.ylabel("photocurrent {in microampere}")
    plt.title("photocurrent V\s wavelength SECTION-B QUES-4")
    plt.grid(b=True,which="both",axis="both")
    plt.legend()
    plt.show()
    print("2019PHY1073" )
    print("SECTION-B QUES-4")
    print("X :",x)
    print("y_mean :" , y_star)
    print("dev_sq:" , dev_sq)
    print("slope:",B)
    print("y-intercept :",A )
    print("error in intercept:",sigma_A)
    print("error in slope :" , sigma_B)
    return 

def lsf2(x,y,pe):
    if len(x)!=len(y) or len(pe)!=len(x) or len(y)!=len(pe):
        print("unequal lengths")
    n=len(x)
    pe_max=max(pe)
    lamb=pe_max**2   # CHECK THE FORMULA (PEMAX/PE[I])^2
    weights=[]
    for i in range(0,len(pe)):
        temp_val=(lamb/pe[i])
        weights.append(temp_val)
    weights=np.array(weights)
    x=np.array(x)
    y=np.array(y)
    x_square=np.square(x)
    x_sq_weight=x_square*weights
    y_square=np.square(y)
    y_sq_weight=y_square*weights
    w_x_y=x*y*weights
    w_x=x*weights
    w_y=y*weights
    x_y=x*y
    ############ SUM  #################
    weight_sum=np.sum(weights)
    sum_x_sq_w=np.sum(x_sq_weight)
    sum_y_sq_w=np.sum(y_sq_weight)
    sum_w_x=np.sum(w_x)
    sum_x=np.sum(x)
    sum_w_y=np.sum(w_y)
    sum_w_x_y=np.sum(w_x_y)
    sum_y=np.sum(y)
    sum_x_sq=np.sum(x_square)
    sum_x_y=np.sum(x_y)
    
    #--------------DELTAS WITHOUT WEIGHTS------------------#
    delta_without_weights=n*sum_x_sq-sum_x*sum_x
    delta_o_without_weights=sum_x_sq*sum_y-sum_x*sum_x_y
    delta1_without_weights=n*sum_x_y-sum_x*sum_y
    A_without_weights=delta_o_without_weights/delta_without_weights
    B_without_weights=delta1_without_weights/delta_without_weights
    Y_fit_without_weights=B_without_weights*x+A_without_weights
    #------------- DELTAS WITH WEIGHTS --------------------#
    delta=weight_sum*sum_x_sq_w-sum_w_x*sum_w_x
    delta1=weight_sum*sum_w_x_y-sum_w_x*sum_w_y
    delta_o=sum_w_y*sum_x_sq_w-sum_w_x_y*sum_w_x
    A=delta_o/delta
    B=delta1/delta
    Y_fitted=B*x+A
    sigma=sum_y_sq_w/(n-1)
    # Y_plus_sigma=(B+sigma)*x+(A+sigma)
    # Y_minus_sigma=(B-sigma)*x+(A-sigma)
    #-------------     PLOTTING ----------------#
    plt.minorticks_on()
    plt.grid(b=True,which="both",axis="both")
    plt.title("The Least Square Fitting")
    plt.xlabel("X values")
    plt.ylabel("Y values")
    plt.plot(x,Y_fitted,label="with weights")
    plt.plot(x,Y_fit_without_weights,label="withoutweights")
    plt.scatter(x,y,color="black",linewidths=0.1,label="original pts")
    plt.errorbar(x,y,yerr=pe,color="red",fmt="o")
    # plt.plot(x,Y_plus_sigma)
    # plt.plot(x,Y_minus_sigma)
    plt.legend()
    plt.show()
    print("X VALUES :",x)
    print("Y VALUES :",y)
    print("SLOPE :",B)
    print("INTERCEPT : ",A)
    print("SIGMA :",sigma)
    return A 



x=[185,195,202,210,230,260,280,300,340,380,420]
y=[[52.4],[46.6],[46.7],[42.5],[38.7],[30.1],[25.9],[22.3],[16.4],[11.85],[8.2]]

x=list(map(float,x))
for i in range(0,len(x)):
    y[i]=list(map(float,y[i]))
print(lsf1(x,y))



"""
                    OUTPUT LSF WITHOUT WEIGHTS
                    
X : [1.2 2.5 6.3]
y_mean : [ 2.83333333  4.43       19.33333333]
dev_sq: [11.166666666666668, 11.4794, 190.60666666666668]
slope: 3.390088593576966
y-intercept : -2.434739756367666
sigma_A: 1.0603524948732654
sigma_B : 0.26681691701343874
"""
"""
                       OUTPUT LSF WITH WEIGHTS
X VALUES : [1874.8 1879.5 1882.8 1902.4 1906.  1926.5 1928.  1932.5 1936.8 1940. ]
Y VALUES : [299989 299911 299855 299902 299794 299778 299796 299775 299777 299776]
SLOPE : -1.693012417471396
INTERCEPT :  303051.3211558141
SIGMA : -1194863775021.164
"""