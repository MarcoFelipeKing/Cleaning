# -*- coding: utf-8 -*-
"""
Created on Tue May 23 15:31:53 2017

@author: Jonty
"""
import random
import math
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

#Write a function that caluclates the absolute distances between two lists of the same length
def Distance(x,y,sd):
    # computes the Euclidean distance between two lists of the same length, normalised by values at sd
    if len(x) == len(y):
        return np.absolute(sum([(((x[i]-y[i])/sd[i])) for i in range(len(x))]))
    else:
        return 'lists not the same length'

def num_at_x(x,L,T):
    for i in range(len(T)-1):
        if T[i] < x < T[i+1]:
            return L[i]

def ode_model(contamination,t,r,C,m,g,l):
	
	Contamination = contamination;
	
	return(r*Contamination-m*math.exp(-g*t)*Contamination)

def deterministic_run(precision,initial_contamination,r,C,m,g,l):
    tmax = 24
            
    time_space = np.linspace(0,tmax,precision+1)
    
    sim=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
    
    num_at_0=initial_contamination
    num_at_1=sim[int(precision*1.0/tmax)]
    num_at_2=sim[int(precision*2.0/tmax)]
    num_at_4=sim[int(precision*4.0/tmax)]
    num_at_8=sim[int(precision*8.0/tmax)]
   #num_at_12=sim[int(precision*12.0/tmax)]
    num_at_24=sim[int(precision*24.0/tmax)]
    #num_at_48=sim[int(precision*48.0/tmax)]
    
    
    return([num_at_0,num_at_1,num_at_2,num_at_4,num_at_8,num_at_24])


Detergent_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)] #surface, phase
Detergent_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Disinfectant_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Disinfectant_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Control_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Control_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

Distilled_Means=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]
Distilled_SD=[[np.zeros((1,6)) for i in range(1)] for j in range(1)]

# Locker: 1
# Right BR: 2
# Table: 3
# Left BR: 4

# We store Detergent[Surface][Phase]

# Detergent for the Locker. Phases 1,2,3. 
Detergent_Means[0][0] = [59.4,18.7,5.4,5.4,2.4,8.6]
Detergent_SD[0][0] = [91.8,26.2,2.30,4.67,4.34,4.28]
# Detergent_Means[0][1] = [5.608333333,4.85,4.25,3.608333333,5.05,3.65,5.05,5.45]
# Detergent_SD[0][1] = [2.777436282,2.989637275,3.051285766,3.032640915,2.940795107,2.989637275,2.940795107,2.796549598]
# Detergent_Means[0][2] = [6.875,3.85,3.65,3.05,2.85,3.25,4.25,6.05]
# Detergent_SD[0][2] = [4.264301495,3.024041598,2.989637275,2.796549598,2.698658671,2.876779809,3.051285766,2.441028613]

# Detergent for Right BR
Disinfectant_Means[0][0] = [96.0,94.3,96.2,1.75,1.6,8.5]
Disinfectant_SD[0][0] = [87.0,87.0,89,0.5,2.3,4.04]
# Detergent_Means[1][1]  = [5.641666667,5.033333333,3.208333333,3.166666667,3.408333333,3.975,5.075,4.991666667]
# Detergent_SD[1][1]  = [7.479267321,5.002499375,3.102071055,2.9871492,3.417311405,7.372391415,4.960751125,5.058238415]
# Detergent_Means[1][2]  = [6.875,3.85,3.65,3.05,2.85,3.25,4.25,6.05]
# Detergent_SD[1][2]  = [4.264301495,3.024041598,2.989637275,2.796549598,2.698658671,2.876779809,3.051285766,2.441028613]

Distilled_Means[0][0] = [109.0,159.3,47,47.2,128,56]
Distilled_SD[0][0] = [107,76.2,9,9,78,76]

Control_Means[0][0] = [59.4,41,52,18.6,21,16.5]
Control_SD[0][0] = [23.5,9.9,21.7,13.2,16.2,6.54]

# # Detergent for Table
# Detergent_Means[2][0] = [13.60833333,10.625,8.391666667,6.433333333,6.9,6.033333333,7.683333333,7.725]
# Detergent_SD[2][0] = [11.7817985,9.931467753,8.608210796,4.754195062,5.89526037,4.792684799,6.743882924,6.697496959]
# Detergent_Means[2][1] = [8.325,8.525,5.45,4.65,7.75,5.7,7.683333333,5.275]
# Detergent_SD[2][1] = [6.911826654,7.23756365,3.051285766,3.204264848,8.958621353,7.453592052,7.192303118,4.925172848]
# Detergent_Means[2][2] = [8.191666667,5.275,5.675,5.475,6.075,4.45,6.075,8.766666667]
# Detergent_SD[2][2] = [8.414055871,5.130915425,5.23964496,5.46292362,5.711943023,3.280033642,5.031152949,7.456829778]

# # Detergent for Left BR
# Detergent_Means[3][0] = [6.058333333,6.525,5.283333333,3.75,3.508333333,1.6,3.925,3.2]
# Detergent_SD[3][0] = [6.937746578,7.212148515,7.539074458,5.273285798,5.727862182,4.207444142,7.553473168,5.735551562]
# Detergent_Means[3][1] = [4.175,5.283333333,5.141666667,2.708333333,4.233333333,2.341666667,3.475,3.675]
# Detergent_SD[3][1] = [6.661202862,8.293848988,7.284992841,4.994170452,5.129781201,5.080939704,5.220855063,5.313166794]
# Detergent_Means[3][2] = [4.958333333,4.291666667,3.35,2.408333333,2.808333333,2.166666667,2.25,4.566666667]
# Detergent_SD[3][2] = [8.597501141,6.563318914,5.050264588,2.557359225,2.733965795,2.580597365,3.051285766,3.479620471]


# # We store Disinfectant[Surface][Phase]

# # Disinfectant for the Locker. Phases 1,2,3. 
# Disinfectant_Means[0][0] = [4.275,0.333333333,0.875,1.166666667,1.25,2.65,5.675,2.408333333]
# Disinfectant_SD[0][0] = [5.031088694,0.562220556,0.5826145,0.317135165,0.72,2.581098403,4.817868144,2.472681196]
# Disinfectant_Means[0][1] = [4.275,1.466666667,1.075,1.4,4.05,4.25,7.25,4.875]
# Disinfectant_SD[0][1] = [5.031088694,2.381441657,1.301772664,1.66816025,3.04449758,3.051285766,2.74,4.990832112]
# Disinfectant_Means[0][2] = [3.95,2.6,2.525,2.85,4.3,3.208333333,2.366666667,3.633333333]
# Disinfectant_SD[0][2] = [5.114381346,2.892916438,2.677613129,2.698658671,7.306232578,2.915537548,2.503216322,5.023107524]

# # Disinfectant for Right BR
# Disinfectant_Means[1][0] = [2.408333333,0.783333333,1.691666667,0.958333333,1.608333333,3.3,2.808333333,2.65]
# Disinfectant_SD[1][0] = [2.472681196,1.370460112,4.787143042,0.537728834,1.550421708,6.354607343,2.733637318,2.581098403]
# Disinfectant_Means[1][1]  = [4.725,2.058333333,1.475,1.641666667,2.525,2.45,7.25,4.458333333]
# Disinfectant_SD[1][1]  = [8.315790443,4.844425935,2.038794015,1.95827525,2.677613129,2.441028613,2.74,6.448430487]
# Disinfectant_Means[1][2]  = [2.441666667,1.716666667,2.083333333,1.883333333,2.483333333,2.05,1.2,1.641666667]
# Disinfectant_SD[1][2]  = [2.736789062,2.274041417,2.388309651,2.183138098,2.707694443,2.074475422,1.249827574,1.95827525]

# # Disinfectant for Table
# Disinfectant_Means[2][0] = [5.525,4.158333333,3.35,4.683333333,2.875,2.85,5.433333333,4.433333333]
# Disinfectant_SD[2][0] = [7.458895985,7.652327608,5.041038481,7.578592051,4.832875028,2.698658671,4.916954015,5.05867867]
# Disinfectant_Means[2][1] = [3.633333333,4.266666667,2.966666667,2.725,2.85,4.85,7.25,4.5]
# Disinfectant_SD[2][1] = [5.023107524,5.184082623,2.869078109,2.800361738,2.698658671,2.989637275,2.74,6.422616289]
# Disinfectant_Means[2][2] = [6.325,3.675,3.075,4.075,3.508333333,3.833333333,2.991666667,3.3]
# Disinfectant_SD[2][2] = [7.270675039,4.997046542,4.887153286,5.028003477,5.098378433,5.044343594,4.929452011,7.242237223]

# # Disinfectant for Left BR
# Disinfectant_Means[3][0] = [4.908333333,0.658333333,0.991666667,1.041666667,1.283333333,1.325,2.008333333,1.808333333]
# Disinfectant_SD[3][0] = [9.233205474,1.379472871,1.328538433,0.473811277,1.206543462,1.182031419,2.103448865,1.858906025]
# Disinfectant_Means[3][1] = [4.1,1.508333333,2.841666667,2.866666667,2.991666667,2.966666667,7.25,5.308333333]
# Disinfectant_SD[3][1] = [7.30481655,2.365779136,6.534551185,4.989529266,4.929452011,2.869078109,2.74,9.194658151]
# Disinfectant_Means[3][2] = [3.15,2,2.775,1.766666667,1.808333333,2.408333333,2.35,3.733333333]
# Disinfectant_SD[3][2] = [4.999827583,2.438908744,7.215987399,1.88566888,1.858906025,2.472681196,4.745051869,6.502961306]

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

while len(parameter_sample) < sample_size:
	# The prior distributions we use are m ~ U(10^(-5),1.0), C ~ U(2,15), r ~ U(10^(-5),1.0), g ~ U(10^(-5),1.0), l ~ U(10^(-5),1.0)
    # We begin by sampling from these distributions and simulating the process
	trial_r = random.uniform(0.0001,1.0)
	trial_C = random.uniform(2.0,120.0)
	trial_l = random.uniform(0.0001,1.0)
	
	# m and g for detergent
	trial_m_de = random.uniform(0.0001,1.0)
	trial_g_de = random.uniform(0.0001,1.0)
	
	# m and g for disinfectant
	trial_m_di = random.uniform(0.0001,1.0)
	trial_g_di = random.uniform(0.0001,1.0)
	
	total_trials+=1.0
    
	euclidean_distance=0
    
	delta = 5.0
    
	# Learning from data for detergent
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Detergent_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_de,trial_g_de,trial_l)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			
			euclidean_distance += Distance(one_run,Detergent_Means[surface][phase],Detergent_SD[surface][phase])

	# Learning from data for disinfectant
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Disinfectant_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_di,trial_g_di,trial_l)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			
			euclidean_distance += Distance(one_run,Disinfectant_Means[surface][phase],Disinfectant_SD[surface][phase])
			
	if euclidean_distance < delta:
		parameter_sample.append([trial_r,trial_C,trial_m_de,trial_g_de,trial_m_di,trial_g_di,trial_l])
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
		Posterior.write(str(trial_l))
		Posterior.write("\n")

#print(parameter_sample)
print("Percentage of trials accepted: ",100*accepted_trials/total_trials)
#print(np.argsort(distances)[:10])

posterior_r=[]
posterior_C=[]
posterior_m_de=[]
posterior_g_de=[]
posterior_m_di=[]
posterior_g_di=[]
posterior_l=[]

for i in range(len(parameter_sample)):
	posterior_r.append(parameter_sample[i][0])
	posterior_C.append(parameter_sample[i][1])
	posterior_m_de.append(parameter_sample[i][2])
	posterior_g_de.append(parameter_sample[i][3])
	posterior_m_di.append(parameter_sample[i][4])
	posterior_g_di.append(parameter_sample[i][5])
	posterior_l.append(parameter_sample[i][6])

# Plot histograms for each prarameter: r, C, m_de, g_de, m_di, g_di and l



# We plot the posteriors
f, ax = plt.subplots(3,3)

ax[0,0].set_title('Histogram for r')
ax[0,0].hist(posterior_r)

ax[0,1].set_title('Histogram for C')
ax[0,1].hist(posterior_C)

ax[1,0].set_title('Histogram for mu-detergent')
ax[1,0].hist(posterior_m_de)

ax[1,1].set_title('Histogram for gamma-detergent')
ax[1,1].hist(posterior_g_de)

ax[1,0].set_title('Histogram for mu-disinfectant')
ax[1,0].hist(posterior_m_di)

ax[1,1].set_title('Histogram for gamma-disinfectant')
ax[1,1].hist(posterior_g_di)

ax[2,0].set_title('Histogram for lambda')
ax[2,0].hist(posterior_l)

plt.show()


df = pd.DataFrame(parameter_sample)
sns.pairplot(df, diag_kind = 'kde',
             plot_kws = {'alpha': 0.6, 's': 80, 'edgecolor': 'k'},
             size = 4)
plt.show()
#posterior_l=np.array(posterior_l)
#posterior_r=np.array(posterior_r)
#posterior_C=np.array(posterior_C)
#posterior_m=np.array(posterior_m)
#posterior_g=np.array(posterior_g)

## plotting a bivariate histogram for lambda and r
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_l,posterior_r,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_l),max(posterior_l)) 
#plt.ylim(min(posterior_r),max(posterior_r)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$\lambda$ ($hours^{-1}$)',ylabel='$r$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for lambda and C
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_l,posterior_C,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_l),max(posterior_l)) 
#plt.ylim(min(posterior_C),max(posterior_C)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$\lambda$ ($hours^{-1}$)',ylabel='$C$ ($bacteria$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for lambda and mu
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_l,posterior_m,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_l),max(posterior_l)) 
#plt.ylim(min(posterior_m),max(posterior_m)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$\lambda$ ($hours^{-1}$)',ylabel='$\mu$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for lambda and gamma
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_l,posterior_g,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_l),max(posterior_l)) 
#plt.ylim(min(posterior_g),max(posterior_g)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$\lambda$ ($hours^{-1}$)',ylabel='$\gamma$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for r and mu
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_r,posterior_m,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_r),max(posterior_r)) 
#plt.ylim(min(posterior_m),max(posterior_m)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$r$ ($hours^{-1}$)',ylabel='$\mu$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for r and C
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_r,posterior_C,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_r),max(posterior_r)) 
#plt.ylim(min(posterior_C),max(posterior_C)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$r$ ($hours^{-1}$)',ylabel='$C$ (bacteria)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for r and g
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_r,posterior_g,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_r),max(posterior_r)) 
#plt.ylim(min(posterior_g),max(posterior_g)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$r$ ($hours^{-1}$)',ylabel='$\gamma$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for C and mu
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_C,posterior_m,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_C),max(posterior_C)) 
#plt.ylim(min(posterior_m),max(posterior_m)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$C$ (bacteria)',ylabel='$\mu$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for C and gamma
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_C,posterior_g,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_C),max(posterior_C)) 
#plt.ylim(min(posterior_g),max(posterior_g)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$C$ (Bacteria)',ylabel='$\gamma$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## plotting a bivariate histogram for mu and gamma
#sns.set(style="ticks", color_codes=True) 
#g = sns.JointGrid(posterior_m,posterior_g,space=0) 
#p1 = g.plot_joint(plt.hexbin,gridsize=50,cmap='Greens',color='0', edgecolor='white') 
#plt.xlim(min(posterior_m),max(posterior_m)) 
#plt.ylim(min(posterior_g),max(posterior_g)) 
#p2 = g.plot_marginals(sns.distplot, kde=False, color='g') 
#g.set_axis_labels(xlabel='$\mu$ ($hours^{-1}$)',ylabel='$\gamma$ ($hours^{-1}$)',fontsize=20) 
#g.fig.set_size_inches(9,8) 
#plt.colorbar() 
#plt.show()
#plt.close('all')

## We plot our 'best' 10 simulations versus the scatter plots of the data
#tmax = 50
#time_space = np.linspace(0,tmax,precision+1)

#fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)

#ax_Locker=axs[0,0]
#ax_Locker.errorbar([0,1,2,4,8,12,24,48],Data_Means[0][0],yerr=Data_SD[0][0], fmt='o', label='Phase 1')
#ax_Locker.errorbar([0,1,2,4,8,12,24,48],Data_Means[0][1],yerr=Data_SD[0][1], fmt='o', label='Phase 2')
#ax_Locker.errorbar([0,1,2,4,8,12,24,48],Data_Means[0][2],yerr=Data_SD[0][2], fmt='o', label='Phase 3')
#ax_Locker.set_xlim(-1,50)
#ax_Locker.set_title('Locker')

#i=np.argsort(distances)[0]
#r=parameter_sample[i][0]
#C=parameter_sample[i][1]
#m=parameter_sample[i][2]
#g=parameter_sample[i][3]
#l=parameter_sample[i][4]

#for phase in range(3):
	#initial_contamination=Data_Means[0][phase][0]
	#my_simulation=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
	#ax_Locker.plot(time_space,my_simulation)
		
#ax_RightBR=axs[0,1]
#ax_RightBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[1][0],yerr=Data_SD[1][0], fmt='o', label='Phase 1')
#ax_RightBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[1][1],yerr=Data_SD[1][1], fmt='o', label='Phase 2')
#ax_RightBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[1][2],yerr=Data_SD[1][2], fmt='o', label='Phase 3')
#ax_RightBR.set_xlim(-1,50)
#ax_RightBR.set_title('Right Bedrail')

#for phase in range(3):
	#initial_contamination=Data_Means[1][phase][0]
	#my_simulation=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
	#ax_RightBR.plot(time_space,my_simulation)

#ax_Table=axs[1,0]
#ax_Table.errorbar([0,1,2,4,8,12,24,48],Data_Means[2][0],yerr=Data_SD[2][0], fmt='o', label='Phase 1')
#ax_Table.errorbar([0,1,2,4,8,12,24,48],Data_Means[2][1],yerr=Data_SD[2][1], fmt='o', label='Phase 2')
#ax_Table.errorbar([0,1,2,4,8,12,24,48],Data_Means[2][2],yerr=Data_SD[2][2], fmt='o', label='Phase 3')
#ax_Table.set_xlim(-1,50)
#ax_Table.set_title('Table')

#for phase in range(3):
	#initial_contamination=Data_Means[2][phase][0]
	#my_simulation=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
	#ax_Table.plot(time_space,my_simulation)

#ax_LeftBR=axs[1,1]
#ax_LeftBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[3][0],yerr=Data_SD[3][0], fmt='o', label='Phase 1')
#ax_LeftBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[3][1],yerr=Data_SD[3][1], fmt='o', label='Phase 2')
#ax_LeftBR.errorbar([0,1,2,4,8,12,24,48],Data_Means[3][2],yerr=Data_SD[3][2], fmt='o', label='Phase 3')
#ax_LeftBR.set_xlim(-1,50)
#ax_LeftBR.set_title('Left Bedrail')

#for phase in range(3):
	#initial_contamination=Data_Means[3][phase][0]
	#my_simulation=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
	#ax_LeftBR.plot(time_space,my_simulation)

#fig.suptitle('Detergent')
#plt.show()

