# -*- coding: utf-8 -*-
"""
Created on Tue April 23 15:31:53 2022
Updated on Monday May 09 2022	

@author: MFK
"""
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

def ode_model(contamination,t,r,C,m,g,die_off):
	
	Contamination = contamination;
	return(r*(1-Contamination/C)*Contamination-m*math.exp(-g*t)*Contamination-die_off*Contamination)

# Extract specific time-ppints from ODE
def deterministic_run(precision,initial_contamination,r,C,m,g,die_off):
    tmax = 24
            
    time_space = np.linspace(0,tmax,precision+1)
    
    sim=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,die_off))
    
    num_at_0=initial_contamination
    num_at_1=sim[int(precision*1.0/tmax)]
    num_at_2=sim[int(precision*2.0/tmax)]
    num_at_4=sim[int(precision*4.0/tmax)]
    num_at_8=sim[int(precision*8.0/tmax)]
    num_at_24=sim[int(precision*24.0/tmax)]
    
    return([num_at_0,num_at_1,num_at_2,num_at_4,num_at_8,num_at_24])

# Store Experimental data: Detergent, Disinfectant, Control and Distilled Water

Detergent_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)] #surface, phase
Detergent_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Disinfectant_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Disinfectant_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Control_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Control_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Distilled_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Distilled_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]


# We store Detergent[Surface][Phase]

# Detergent for the Locker. Phases 1,2,3. 
Detergent_Means[0][0] = [221.6,94.3,56.25,1.75,1.6,8.5]
Detergent_SD[0][0] = [76.4,86.9,89.4,0.5,2.3,4.04]

Disinfectant_Means[0][0] = [59.4,18.7,5.4,5.4,2.4,8.6]
Disinfectant_SD[0][0] = [91.8,26.2,2.30,4.67,4.34,4.28] 

Distilled_Means[0][0] = [261.0,175.5,47,18.6,128,56]
Distilled_SD[0][0] = [31.5,61.7,9.0,13.2,78.2,76.4]

Control_Means[0][0] = [59.4,41,52,18.6,21,16.5]
Control_SD[0][0] = [23.5,9.9,21.7,13.2,16.2,6.54]


##################################################################################################################
## Applying the ABC algorithm
sample_size = 1000
parameter_sample = []

total_trials=0.0
accepted_trials=0.0

# File with posterior
Posterior = open("Posterior_Beth_Expt.txt","w")

distances=[]

# Precision of the ode solver
precision=5000

#delta
delta = 20.0

#create function to test different parameters in deterministic_run
# def test_parameters(parameters):


while len(parameter_sample) < sample_size:
	# The prior distributions we use are m ~ U(10^(-5),1.0), C ~ U(2,15), r ~ U(10^(-5),1.0), g ~ U(10^(-5),1.0), l ~ U(10^(-5),1.0)
    # We begin by sampling from these distributions and simulating the process
	trial_r = random.uniform(0.001,1.0)
	trial_C = random.uniform(1.0,70.0)
	trial_die_off = random.uniform(0.0001,20.0)
	
	# m and g for detergent
	trial_m_de = random.uniform(0.01,1.0)
	trial_g_de = random.uniform(0.0001,1.0)
	
	# m and g for disinfectant
	trial_m_di = random.uniform(0.01,1.0)
	trial_g_di = random.uniform(0.0001,1.0)

	# m and g for distilled water
	trial_m_dw = random.uniform(0.01,1.0)
	trial_g_dw = random.uniform(0.0001,1.0)

	# m and g for control 0
	trial_m_c = 0.0
	trial_g_c = 0.0

	total_trials+=1.0
    
	euclidean_distance=0
    
    
	# Learning from data for detergent
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Detergent_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_de,trial_g_de,trial_die_off)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			#Calculate the absolute difference between one_run and Detergent_Means[surface][phase]
			#euclidean_distance += np.sum(np.abs(np.subtract(np.log(one_run),np.log(Detergent_Means))))
			euclidean_distance += Distance(one_run,Detergent_Means[surface][phase],[1,1,1,1,1,1])#)#Detergent_SD[surface][phase]

	# Learning from data for disinfectant
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Disinfectant_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_di,trial_g_di,trial_die_off)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			#euclidean_distance += np.sum(np.abs(np.subtract(np.log(one_run),np.log(Disinfectant_Means))))
			euclidean_distance += Distance(one_run,Disinfectant_Means[surface][phase],[1,1,1,1,1,1])#[1,1,1,1,1,1])#,Disinfectant_SD[surface][phase]
	
	# Learning from data for distilled water
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Distilled_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_dw,trial_g_dw,trial_die_off)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			#euclidean_distance += np.sum(np.abs(np.subtract(np.log(one_run),np.log(Distilled_Means))))
			euclidean_distance += Distance(one_run,Distilled_Means[surface][phase],[1,1,1,1,1,1])#1,1,1,1,1,1]) #Distilled_SD[surface][phase]

# Learning from data for control
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Control_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,0.0,0.0,trial_die_off)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			#euclidean_distance += np.sum(np.abs(np.subtract(np.log(one_run),np.log(Distilled_Means))))
			euclidean_distance += Distance(one_run,Control_Means[surface][phase],[1,1,1,1,1,1])#1,1,1,1,1,1]) #Distilled_SD[surface][phase]

	if euclidean_distance < delta:
		parameter_sample.append([trial_r,trial_C,trial_m_de,trial_g_de,trial_m_di,trial_g_di,trial_m_dw,trial_g_dw,trial_die_off])
		distances.append(euclidean_distance)
		accepted_trials+=1.0
		print(accepted_trials)
		Posterior.write(str(trial_r))
		Posterior.write(",")
		Posterior.write(str(trial_C))
		Posterior.write(",")
		Posterior.write(str(trial_m_de))
		Posterior.write(",")
		Posterior.write(str(trial_g_de))
		Posterior.write(",")
		Posterior.write(str(trial_m_di))
		Posterior.write(",")
		Posterior.write(str(trial_g_di))
		Posterior.write(",")
		Posterior.write(str(trial_m_dw))
		Posterior.write(",")
		Posterior.write(str(trial_g_dw))
		Posterior.write(",")
		Posterior.write(str(trial_die_off))
		Posterior.write("\n")

#print(parameter_sample)
print("Percentage of trials accepted: ",100*accepted_trials/total_trials)
#print(np.argsort(distances)[:10])

# posterior_r=[]
# posterior_C=[]
# posterior_m_de=[]
# posterior_g_de=[]
# posterior_m_di=[]
# posterior_g_di=[]
# posterior_m_dw=[]
# posterior_g_dw=[]
# # posterior_l=[]

# for i in range(len(parameter_sample)):
# 	posterior_r.append(parameter_sample[i][0])
# 	posterior_C.append(parameter_sample[i][1])
# 	posterior_m_de.append(parameter_sample[i][2])
# 	posterior_g_de.append(parameter_sample[i][3])
# 	posterior_m_di.append(parameter_sample[i][4])
# 	posterior_g_di.append(parameter_sample[i][5])
# 	posterior_m_dw.append(parameter_sample[i][6])
# 	posterior_g_dw.append(parameter_sample[i][7])
# 	# posterior_l.append(parameter_sample[i][8])
