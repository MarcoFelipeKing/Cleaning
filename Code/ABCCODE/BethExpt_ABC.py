# -*- coding: utf-8 -*-
"""
Created on Tue May 23 15:31:53 2022
Updated on Thurs April 28 2020

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
        return np.absolute(sum([(((x[i]-y[i])/sd[i])) for i in range(len(x))]))
    else:
        return 'lists not the same length'

def num_at_x(x,L,T):
    for i in range(len(T)-1):
        if T[i] < x < T[i+1]:
            return L[i]

def ode_model(contamination,t,r,C,m,g,l):
	
	Contamination = contamination;
	
	return((Contamination*r+l)*(1-Contamination/C)-m*math.exp(-g*t)*Contamination)

def deterministic_run(precision,initial_contamination,r,C,m,g,l):
    tmax = 24
            
    time_space = np.linspace(0,tmax,precision+1)
    
    sim=odeint(ode_model,initial_contamination,time_space,args=(r,C,m,g,l))
    
    num_at_0=initial_contamination
    num_at_1=sim[int(precision*1.0/tmax)]
    num_at_2=sim[int(precision*2.0/tmax)]
    num_at_4=sim[int(precision*4.0/tmax)]
    num_at_8=sim[int(precision*8.0/tmax)]
    num_at_24=sim[int(precision*24.0/tmax)]
    
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
Detergent_Means[0][0] = [221.6,94.3,96.25,1.75,1.6,8.5]
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

while len(parameter_sample) < sample_size:
	# The prior distributions we use are m ~ U(10^(-5),1.0), C ~ U(2,15), r ~ U(10^(-5),1.0), g ~ U(10^(-5),1.0), l ~ U(10^(-5),1.0)
    # We begin by sampling from these distributions and simulating the process
	trial_r = random.uniform(0.0001,20.0)
	trial_C = random.uniform(100.0,420.0)
	trial_l = random.uniform(0.0001,20.0)
	
	# m and g for detergent
	trial_m_de = random.uniform(0.0001,20.0)
	trial_g_de = random.uniform(0.0001,1.0)
	
	# m and g for disinfectant
	trial_m_di = random.uniform(0.0001,20.0)
	trial_g_di = random.uniform(0.0001,20.0)

	# m and g for distilled water
	trial_m_dw = random.uniform(0.0001,20.0)
	trial_g_dw = random.uniform(0.0001,20.0)
	
	total_trials+=1.0
    
	euclidean_distance=0
    
	delta = 100.0
    
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
			
			euclidean_distance += Distance(one_run,Disinfectant_Means[surface][phase],Distilled_Means[surface][phase])
	
	# Learning from data for distilled water
	for surface in range(1):
		for phase in range(1):
			initial_contamination=Distilled_Means[surface][phase][0]
			one_run = deterministic_run(precision,initial_contamination,trial_r,trial_C,trial_m_dw,trial_g_dw,trial_l)
			# Now we find the Euclidean distance between the simulated output and the
			# experimental results, normalised by the sd of the data. delta is the threshold that the Euclidean distance
			# must be less than for us to accept the trial parameters into our sample.
			
			euclidean_distance += Distance(one_run,Distilled_Means[surface][phase],Distilled_Means[surface][phase])
			
	if euclidean_distance < delta:
		parameter_sample.append([trial_r,trial_C,trial_m_de,trial_g_de,trial_m_di,trial_g_di,trial_m_dw,trial_g_dw,trial_l])
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
posterior_m_dw=[]
posterior_g_dw=[]
posterior_l=[]

for i in range(len(parameter_sample)):
	posterior_r.append(parameter_sample[i][0])
	posterior_C.append(parameter_sample[i][1])
	posterior_m_de.append(parameter_sample[i][2])
	posterior_g_de.append(parameter_sample[i][3])
	posterior_m_di.append(parameter_sample[i][4])
	posterior_g_di.append(parameter_sample[i][5])
	posterior_m_dw.append(parameter_sample[i][6])
	posterior_g_dw.append(parameter_sample[i][7])
	posterior_l.append(parameter_sample[i][8])

# Plot histograms for each prarameter: r, C, m_de, g_de, m_di, g_di, m_ds, g_ds, and l



# We plot the posteriors
f, ax = plt.subplots(3,3)

ax[0,0].set_title('r')
ax[0,0].hist(posterior_r)

ax[0,1].set_title('C')
ax[0,1].hist(posterior_C)

ax[1,0].set_title('mu-detergent')
ax[1,0].hist(posterior_m_de)

ax[1,1].set_title('gamma-detergent')
ax[1,1].hist(posterior_g_de)

ax[1,0].set_title('mu-disinfectant')
ax[1,0].hist(posterior_m_di)

ax[1,1].set_title('gamma-disinfectant')
ax[1,1].hist(posterior_g_di)

ax[1,0].set_title('mu-distilled-water')
ax[1,0].hist(posterior_m_di)

ax[1,1].set_title('gamma-disinfectant')
ax[1,1].hist(posterior_g_di)

ax[2,0].set_title('Histogram for lambda')
ax[2,0].hist(posterior_l)

plt.show()



#arrange parameter_sample by distances using argsort
# arr1inds = distances.argsort()
# parameter_sample_sorted = [parameter_sample[i] for i in arr1inds]
 

# We plot the posteriors using seaborn
df = pd.DataFrame({'r': parameter_sample[0], 'C': parameter_sample[1], 'mu-detergent': parameter_sample[2], 'gamma-detergent': parameter_sample[3], 'mu-disinfectant': parameter_sample[4], 'gamma-disinfectant': parameter_sample[5], 'mu-distilled-water': parameter_sample[6], 'gamma-distilled-water': parameter_sample[7], 'lambda': parameter_sample[8]})
df.to_csv("abc_results_BethExpt.csv", encoding='utf-8', index=False)
df.describe()

# sns.pairplot(df, diag_kind = 'kde',
#              plot_kws = {'alpha': 0.4, 's': 80, 'edgecolor': 'k'},
#              size = 0.5)
# plt.show()
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
# tmax = 24
# time_space = np.linspace(0,tmax,5000+1)

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

# for phase in range(1):

#Plot errorbars of experiments with predictions ontop

precision=5000
tmax = 24
time_space = np.linspace(0,tmax,precision+1)

#1. Detergent
initial_contamination=Detergent_Means[0][0][0]
P=odeint(ode_model,initial_contamination,time_space,args=(df["r"].mean(),df["C"].mean(),df["mu-detergent"].mean(),df["gamma-detergent"].mean(),df["lambda"].mean()))
Pmin=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.05),df["C"].quantile(0.05),df["mu-detergent"].quantile(0.05),df["gamma-detergent"].quantile(0.05),df["lambda"].quantile(0.05)))
Pmax=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.95),df["C"].quantile(0.95),df["mu-detergent"].quantile(0.95),df["gamma-detergent"].quantile(0.95),df["lambda"].quantile(0.95)))

# create a vector of values between 0 and 6
x = np.array([0,2,4,8,12,24])

#Define new sd just for plotting to avoid SD value at 0
s = Detergent_SD[0][0]
measurement_data = Detergent_Means[0][0]#np.array([1200,134.0,202.0,294.0])

#Plot errobars of experimental data
plt.errorbar(x,measurement_data,yerr=s,fmt='o', color='black',label='Experimental data')

#Plot the model prediction
plt.plot(time_space,P,label="Model prediction",color='blue')

#Plot confidence intervals around the model prediction

plt.fill_between(time_space,np.concatenate(Pmin),np.concatenate(Pmax),alpha=0.2,color='blue')
#plt.plot(time_space,Pmin,label="Model prediction",color='red')
#plt.plot(time_space,Pmax,label="Model prediction",color='red')

#plt.fill_between(x, np.array(map(operator.sub, P["Contamination"], Pmin["Contamination"])), np.array(map(operator.add, P["Contamination"], Pmax["Contamination"])), color='b', alpha=.1)
plt.xlim(-1,25)
plt.ylabel("CFU recovered from coupon \n after t hours")
plt.yscale("log")
plt.xlabel("Hours after surface cleaning")
plt.legend(loc="upper right")
plt.title("Detergent")
#make y axis logarithmic


#save the plot
plt.savefig("abc_prediction_BethExpt_detergent.png", dpi=600)

plt.show()

#2. Disinfectant
initial_contamination=Disinfectant_Means[0][0][0]
P=odeint(ode_model,initial_contamination,time_space,args=(df["r"].mean(),df["C"].mean(),df["mu-disinfectant"].mean(),df["gamma-disinfectant"].mean(),df["lambda"].mean()))
Pmin=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.05),df["C"].quantile(0.05),df["mu-disinfectant"].quantile(0.05),df["gamma-disinfectant"].quantile(0.05),df["lambda"].quantile(0.05)))
Pmax=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.95),df["C"].quantile(0.95),df["mu-disinfectant"].quantile(0.95),df["gamma-disinfectant"].quantile(0.95),df["lambda"].quantile(0.95)))

# create a vector of values between 0 and 6
x = np.array([0,2,4,8,12,24])

#Define new sd just for plotting to avoid SD value at 0
s = Disinfectant_SD[0][0]
measurement_data = Disinfectant_Means[0][0]#np.array([1200,134.0,202.0,294.0])

#Plot errobars of experimental data
plt.errorbar(x,measurement_data,yerr=s,fmt='o', color='black',label='Experimental data')

#Plot the model prediction
plt.plot(time_space,P,label="Model prediction",color='blue')

#Plot confidence intervals around the model prediction

plt.fill_between(time_space,np.concatenate(Pmin),np.concatenate(Pmax),alpha=0.2,color='red')
#plt.plot(time_space,Pmin,label="Model prediction",color='red')
#plt.plot(time_space,Pmax,label="Model prediction",color='red')

#plt.fill_between(x, np.array(map(operator.sub, P["Contamination"], Pmin["Contamination"])), np.array(map(operator.add, P["Contamination"], Pmax["Contamination"])), color='b', alpha=.1)
plt.xlim(-1,25)
plt.ylabel("CFU recovered from coupon \n after t hours")
plt.yscale("log")
plt.xlabel("Hours after surface cleaning")
plt.legend(loc="upper right")
#Add title
plt.title("Disinfectant")
#make y axis logarithmic


#save the plot
plt.savefig("abc_prediction_BethExpt_disinfectant.png", dpi=600)

plt.show()

#3. Distilled water
initial_contamination=Distilled_Means[0][0][0]
P=odeint(ode_model,initial_contamination,time_space,args=(df["r"].mean(),df["C"].mean(),df["mu-distilled-water"].mean(),df["gamma-distilled-water"].mean(),df["lambda"].mean()))
Pmin=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.05),df["C"].quantile(0.05),df["mu-distilled-water"].quantile(0.05),df["gamma-distilled-water"].quantile(0.05),df["lambda"].quantile(0.05)))
Pmax=odeint(ode_model,initial_contamination,time_space,args=(df["r"].quantile(0.95),df["C"].quantile(0.95),df["mu-distilled-water"].quantile(0.95),df["gamma-distilled-water"].quantile(0.95),df["lambda"].quantile(0.95)))

# create a vector of values between 0 and 6
x = np.array([0,2,4,8,12,24])

#Define new sd just for plotting to avoid SD value at 0
s = Distilled_SD[0][0]
measurement_data = Distilled_Means[0][0]#np.array([1200,134.0,202.0,294.0])

#Plot errobars of experimental data
plt.errorbar(x,measurement_data,yerr=s,fmt='o', color='black',label='Experimental data')

#Plot the model prediction
plt.plot(time_space,P,label="Model prediction",color='blue')

#Plot confidence intervals around the model prediction

plt.fill_between(time_space,np.concatenate(Pmin),np.concatenate(Pmax),alpha=0.2,color='green')
#plt.plot(time_space,Pmin,label="Model prediction",color='red')
#plt.plot(time_space,Pmax,label="Model prediction",color='red')

#plt.fill_between(x, np.array(map(operator.sub, P["Contamination"], Pmin["Contamination"])), np.array(map(operator.add, P["Contamination"], Pmax["Contamination"])), color='b', alpha=.1)
plt.xlim(-1,25)
plt.ylabel("CFU recovered from coupon \n after t hours")
plt.yscale("log")
plt.xlabel("Hours after surface cleaning")
plt.legend(loc="upper right")
#Add title
plt.title("Distilled-water")
#make y axis logarithmic


#save the plot
plt.savefig("abc_prediction_BethExpt_distilled-water.png", dpi=600)

plt.show()

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