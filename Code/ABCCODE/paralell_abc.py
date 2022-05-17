# Create N random values for each of the 4 parameters

#Import libraries
import random
import math
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
from sympy import ordered

#Write a function that caluclates the absolute distances between two lists of the same length
def Distance(x,y,sd):
    # computes the Euclidean distance between two lists of the same length, normalised by values at sd
    if len(x) == len(y):
        return np.absolute(sum([(((np.log(x[i])-np.log(y[i]))/sd[i])) for i in range(len(x))]))
    else:
        return 'lists not the same length'

# ODE to calculate numerically

def ode_model(contamination,t,r,C,m,g):
	
	Contamination = contamination;
	return(r*(1-Contamination/C)*Contamination-m*math.exp(-g*t)*Contamination)

# Extract specific time-ppints from ODE
def deterministic_run(precision,initial_contamination,r,C,m,g):
    tmax = 24
            
    time_space = np.linspace(0,tmax,precision+1)
    
    sim=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g))
    
    num_at_0=initial_contamination
    num_at_1=sim[int(precision*1.0/tmax)]
    num_at_2=sim[int(precision*2.0/tmax)]
    num_at_4=sim[int(precision*4.0/tmax)]
    num_at_8=sim[int(precision*8.0/tmax)]
    num_at_24=sim[int(precision*24.0/tmax)]
    
    return([num_at_0,num_at_1,num_at_2,num_at_4,num_at_8,num_at_24])

#set seed for reproducibility
random.seed(234)

# Create N random values for each of the 4 parameters: r, C, m, g
def random_values(N,r_min,r_max,C_min,C_max,m_min,m_max,g_min,g_max):
    r_values = [random.uniform(r_min,r_max) for i in range(N)]
    C_values = [random.uniform(C_min,C_max) for i in range(N)]
    m_values = [random.uniform(m_min,m_max) for i in range(N)]
    g_values = [random.uniform(g_min,g_max) for i in range(N)]
    return(r_values,C_values,m_values,g_values)

# Store random values in a dataframe
def random_values_df(N,r_min,r_max,C_min,C_max,m_min,m_max,g_min,g_max):
    r_values,C_values,m_values,g_values = random_values(N,r_min,r_max,C_min,C_max,m_min,m_max,g_min,g_max)
    rv = pd.DataFrame({'r':r_values,'C':C_values,'m':m_values,'g':g_values})
    return(rv)
rv=random_values_df(100000,0.01,0.1,0.01,0.1,0.01,0.1,0.01,0.1)

# Experimental data
data = [221.6,94.3,56.25,1.75,1.6,8.5]
# Calculate deterministic_run with each parameter combination from rv and calculate the distance between the prediction and the experimental data
def distance_calc(rv,data):
    distance_list = []
    for i in range(len(rv)):
        distance_list.append(Distance(deterministic_run(5000,data[0],rv['r'][i],rv['C'][i],rv['m'][i],rv['g'][i]),data,data))
    return(distance_list)

# test the code
print(distance_calc(rv,data))

# Plot the distance between the prediction and the experimental data
sns.set_style("whitegrid")
sns.set_context("paper")
plt.figure(figsize=(10,10))
plt.plot(distance_calc(rv,data))
plt.xlabel('Parameter Combinations')
plt.ylabel('Distance')
plt.title('Distance between prediction and experimental data')
plt.show()

# Bind distance_calc with rv 
rv['distance']=distance_calc(rv,data)

# Sort the dataframe by distance
rv_sorted = rv.sort_values(by='distance')
